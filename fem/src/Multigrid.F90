!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2001
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!> Module containing various versions of multigrid solvers: geometric, algebraic,
!> clustering and p-multigrid.
!------------------------------------------------------------------------------


MODULE Multigrid

   USE CRSMatrix
   USE IterSolve
   USE DirectSolve
   USE Smoothers
   USE ClusteringMethods

   IMPLICIT NONE


CONTAINS

!------------------------------------------------------------------------------
!> Multigrid solution subroutine common for all different multilevel strategies.
!> Multigrid methods include own developments of geometric, algebraic, clustering
!> and p-element versions.
!------------------------------------------------------------------------------
    RECURSIVE SUBROUTINE MultiGridSolve( Matrix1, Solution, &
        ForceVector, DOFs, Solver, Level, NewSystem )
!------------------------------------------------------------------------------
       USE ModelDescription
       IMPLICIT NONE

       TYPE(Matrix_t), POINTER :: Matrix1
       INTEGER :: DOFs, Level
       LOGICAL, OPTIONAL :: NewSystem
       TYPE(Solver_t) :: Solver
       REAL(KIND=dp) CONTIG :: ForceVector(:), Solution(:)
!------------------------------------------------------------------------------
       LOGICAL :: Found, Algebraic, Cluster, Geometric, Pelement
       CHARACTER(LEN=MAX_NAME_LEN) :: MGMethod
       TYPE(ValueList_t), POINTER :: Params

       IF( Level == Solver % MultigridLevel ) THEN 
         CALL Info('MultiGridSolve','*********************************',Level=10)
         WRITE( Message,'(A,I0)') 'Performing multigrid solution cycle: ',&
             Matrix1 % NumberOfRows
         CALL Info('MultiGridSolve',Message,Level=10 )
       END IF

       Params => Solver % Values
       MGMethod = ListGetString( Params,'MG Method',Found) 
       IF( Found ) THEN
         Pelement = ( MGmethod == 'p' )  
         Cluster = ( MGmethod == 'cluster' ) 
         Algebraic = ( MGmethod == 'algebraic' ) 
         Geometric = ( MGmethod == 'geometric' )  
       ELSE
         Algebraic = ListGetLogical( Params, 'MG Algebraic', Found ) 
         Cluster = ListGetLogical( Params, 'MG Cluster', Found )
         PElement = ListGetLogical( Params, 'MG PElement', Found )
         Geometric = ListGetLogical( Params, 'MG Geometric', Found )
       END IF

       IF ( Algebraic ) THEN
         CALL AMGSolve( Matrix1, Solution, ForceVector, DOFs, Solver, Level, NewSystem )
       ELSE IF( Cluster ) THEN
         CALL CMGSolve( Matrix1, Solution, ForceVector, DOFs, Solver, Level, NewSystem )
       ELSE IF( Pelement ) THEN
         CALL PMGSolve( Matrix1, Solution, ForceVector, DOFs, Solver, Level, NewSystem )
       ELSE
         CALL GMGSolve( Matrix1, Solution, ForceVector, DOFs, Solver, Level, NewSystem )
       END IF
!------------------------------------------------------------------------------
    END SUBROUTINE MultiGridSolve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Geometric multigrid solution procedure.
!------------------------------------------------------------------------------
    RECURSIVE SUBROUTINE GMGSolve( Matrix1, Solution, &
        ForceVector, DOFs, Solver, Level, NewSystem )
!------------------------------------------------------------------------------
       USE ModelDescription
       IMPLICIT NONE

       TYPE(Matrix_t), POINTER :: Matrix1
       INTEGER :: DOFs, Level
       LOGICAL, OPTIONAL :: NewSystem
       TYPE(Solver_t), TARGET :: Solver       
       REAL(KIND=dp), TARGET CONTIG :: ForceVector(:), Solution(:)
!------------------------------------------------------------------------------
       TYPE(Variable_t), POINTER :: Variable1, TimeVar, SaveVariable
       TYPE(Mesh_t), POINTER   :: Mesh1, Mesh2, SaveMesh
       TYPE(Matrix_t), POINTER :: Matrix2, PMatrix, SaveMatrix
       TYPE(Solver_t), POINTER :: PSolver             


       INTEGER :: i,j,k,l,m,n,n2,k1,k2,iter,MaxIter = 100,ndofs
       CHARACTER(LEN=MAX_NAME_LEN) :: Path,str,mgname, LowestSolver
       LOGICAL :: Condition, Found, Parallel, Project,Transient, LIter

       TYPE(Matrix_t), POINTER :: ProjPN, ProjQT
       INTEGER, POINTER :: Permutation(:), Permutation2(:)

       REAL(KIND=dp), POINTER :: Residual2(:), Solution2(:)
       REAL(KIND=dp), ALLOCATABLE, TARGET :: Residual(:)
       REAL(KIND=dp) :: ResidualNorm, RHSNorm, Tolerance, ILUTOL, tmp
       TYPE(ValueList_t), POINTER :: Params

#ifdef USE_ISO_C_BINDINGS
       REAL(KIND=dp) :: tt
#else
       REAL(KIND=dp) :: CPUTime, tt
#endif


       LOGICAL :: NewLinearSystem
       SAVE NewLinearSystem

!------------------------------------------------------------------------------
       tt = CPUTime()
       Params => Solver % Values

!      Initialize:
!      -----------
       Parallel = ParEnv % PEs > 1

       n = Matrix1 % NumberOfRows
       ALLOCATE( Residual(n) )
       Residual = 0.0d0


       RHSNorm = ParallelReduction(SQRT(SUM(ForceVector**2)))
       Solution(1:n) = Solution(1:n) / RHSnorm
       ForceVector(1:n) = ForceVector(1:n) / RHSnorm
!
!      Check for top & bottom levels:
!      ------------------------------ 

       IF ( Level == Solver % MultiGridLevel ) THEN
          NewLinearSystem = .TRUE.
          IF ( PRESENT( NewSystem ) ) THEN
             NewLinearSystem = NewLinearSystem .AND. NewSystem
          END IF
          IF( NewLinearSystem ) THEN
            CALL Info('GMGSolve','Encountered a new linear system',Level=10)
          END IF
          WRITE( Message,'(A,I0)') 'Size of variable field (in partition 1):',&
              SIZE( Solver % Variable % Values )
          CALL Info('GMGSolve', Message, Level=12 )
       ELSE IF ( Level <= 1 ) THEN


          CALL ListPushNamespace('mglowest:')
          
          CALL ListAddLogical( Params,'mglowest: Linear System Free Factorization', .FALSE. )
          CALL ListAddLogical( Params,'mglowest: Linear System Refactorize', NewLinearSystem )

          LowestSolver = ListGetString(Params,'MG Lowest Linear Solver',Found)
          IF ( .NOT. Found ) THEN
            LIter=ListGetLogical(Params,'MG Lowest Linear Solver Iterative',Found)
            IF ( .NOT. Found .AND. Parallel ) LIter=.TRUE.
            LowestSolver='direct'
            IF ( LIter ) LowestSolver='iterative'
          END IF
          
          CALL Info('GMGSolve','Starting lowest linear solver using: '//TRIM(LowestSolver),Level=10 )

          SELECT CASE(LowestSolver)

          CASE('iterative')
            IF ( Parallel ) THEN
              CALL ParallelIter( Matrix1, Matrix1 % ParallelInfo, DOFs, &
                  Solution, ForceVector, Solver, Matrix1 % ParMatrix )
            ELSE
              CALL IterSolver( Matrix1, Solution, ForceVector, Solver )
            END IF

          CASE('direct')
            CALL DirectSolver( Matrix1, Solution, ForceVector, Solver )

          CASE('smoother')
            IF ( Parallel ) THEN
               CALL ParallelInitSolve( Matrix1, Solution, &
                       ForceVector, Residual, NewLinearSystem )
            END IF
            PSolver => Solver
            tmp = MGSmooth( PSolver, Matrix1, Mesh1, Solution, &
                 ForceVector, Residual, Level, DOFs, LowestSmooth = .TRUE. )

          CASE('none') 
            CALL Info('GMGSolve','Applying no solver for coarsest level')

          CASE DEFAULT
             CALL Warn( 'GMGSolve', 'Unknown solver selection for MG lowest level' )
             CALL Warn( 'GMGSolve', 'Using iterative solver' )

             IF ( Parallel ) THEN
               CALL ParallelIter( Matrix1, Matrix1 % ParallelInfo, DOFs, &
                   Solution, ForceVector, Solver, Matrix1 % ParMatrix )
             ELSE
               CALL IterSolver( Matrix1, Solution, ForceVector, Solver )
             END IF
          END SELECT

          CALL ListPopNamespace('mglowest:')
          
          Solution(1:n) = Solution(1:n) * RHSNorm
          ForceVector(1:n) = ForceVector(1:n) * RHSnorm

          DEALLOCATE( Residual )
          NewLinearSystem = .FALSE.

          CALL Info('GMGSolve','Exiting lowest level',Level=12)

          RETURN
       END IF

!      Compute residual:
!      -----------------
       IF ( Parallel ) THEN
          CALL ParallelInitSolve( Matrix1, Solution, &
                  ForceVector, Residual, NewLinearSystem )
          PMatrix => ParallelMatrix( Matrix1 )
       END IF

       CALL MGmv( Matrix1, Solution, Residual, .TRUE. )
       Residual(1:n) = ForceVector(1:n) - Residual(1:n)
 
       Tolerance = ListGetConstReal( Params,'Linear System Convergence Tolerance' )
!---------------------------------------------------------------------
!
!      Initialize the multilevel solve:
!      --------------------------------
       SaveMesh     => Solver % Mesh
       SaveMatrix   => Solver % Matrix
       SaveVariable => Solver % Variable

       Mesh1 => Solver % Mesh
       Variable1   => Solver % Variable
       Permutation => Variable1 % Perm

       Mesh2 => Mesh1 % Parent
!---------------------------------------------------------------------
!
!      Allocate mesh and variable structures for the
!      next level mesh, if not already there:
!      ---------------------------------------------
       IF ( .NOT. ASSOCIATED( Mesh2 ) ) THEN

          mgname = ListGetString( Params, 'MG Mesh Name', Found )
          IF ( .NOT. Found ) mgname = 'mgrid'

          WRITE( Path,'(a,i1)' ) TRIM(OutputPath) // '/' // TRIM(mgname), Level - 1

          Mesh2 => LoadMesh2( CurrentModel, OutputPath, Path, &
               .FALSE., ParEnv % PEs, ParEnv % MyPE )

          CALL UpdateSolverMesh( Solver, Mesh2 )
          CALL ParallelInitMatrix( Solver, Solver % Matrix )

          TimeVar => VariableGet( Mesh1 % Variables, 'Time' )
          CALL VariableAdd( Mesh2 % Variables, Mesh2, Solver, &
                  'Time', 1, TimeVar % Values )

          TimeVar => VariableGet( Mesh1 % Variables, 'Timestep' )
          CALL VariableAdd( Mesh2 % Variables, Mesh2, Solver, &
                  'Timestep', 1, TimeVar % Values )

          TimeVar => VariableGet( Mesh1 % Variables, 'Timestep size' )
          CALL VariableAdd( Mesh2 % Variables, Mesh2, Solver, &
                  'Timestep size', 1, TimeVar % Values )

          TimeVar => VariableGet( Mesh1 % Variables, 'Timestep interval' )
          CALL VariableAdd( Mesh2 % Variables, Mesh2, Solver, &
                  'Timestep interval', 1, TimeVar % Values )

          TimeVar => VariableGet( Mesh1 % Variables, 'Coupled iter' )
          CALL VariableAdd( Mesh2 % Variables, Mesh2, Solver, &
                  'Coupled iter', 1, TimeVar % Values )

          TimeVar => VariableGet( Mesh1 % Variables, 'Nonlin iter' )
          CALL VariableAdd( Mesh2 % Variables, Mesh2, Solver, &
                  'Nonlin iter', 1, TimeVar % Values )

          CALL VariableAdd( Mesh2 % Variables, Mesh2, Solver, &
                'Coordinate 1', 1, Mesh2 % Nodes % x )

          CALL VariableAdd( Mesh2 % Variables, Mesh2, Solver, &
                'Coordinate 2', 1, Mesh2 % Nodes % y )

          CALL VariableAdd( Mesh2 % Variables,Mesh2, Solver, &
                'Coordinate 3', 1, Mesh2 % Nodes % z )

          Matrix2 => Solver % Matrix

          Mesh2 % Child    => Mesh1
          Mesh1 % Parent   => Mesh2

          Matrix2 % Child  => Matrix1
          Matrix1 % Parent => Matrix2

          Permutation2 => Solver % Variable % Perm
       ELSE
         Matrix2 => Matrix1 % Parent
         IF(.NOT. ASSOCIATED( Matrix2 ) ) THEN
           CALL Fatal('GMGSolve','Matrix2 not associated!')
         END IF

         Solver % Mesh => Mesh2
         CALL SetCurrentMesh( CurrentModel, Mesh2 )
         Solver % Variable => VariableGet( Mesh2 % Variables, &
             Variable1 % Name, ThisOnly = .TRUE. )
       END IF
       
 
!------------------------------------------------------------------------------
!
!      Some more initializations:
!      --------------------------
       n  = Matrix1 % NumberOfRows
       n2 = Matrix2 % NumberOfRows

       Residual2 => Matrix2 % RHS
       Solution2 => Solver % Variable % Values
       Permutation2 => Solver % Variable % Perm

!------------------------------------------------------------------------------
!
!      Mesh projector matrices from the higher level mesh
!      to the  lower and  transpose of the projector from
!      the lower level mesh to the higher:
!      ---------------------------------------------------
       ProjPN => MeshProjector( Mesh2, Mesh1, Trans = .TRUE. )
       ProjQT => MeshProjector( Mesh2, Mesh1, Trans = .TRUE. )

       Project = ListGetLogical( Params, 'MG Project Matrix', Found )
       IF ( .NOT. Found ) Project = .TRUE.

       IF ( Project ) THEN
          IF ( NewLinearSystem ) THEN
!            Project higher  level coefficient matrix to the
!            lower level: A_low = ProjPN * A_high * ProjQT^T
!            -----------------------------------------------
             CALL ProjectMatrix( Matrix1, Permutation, ProjPN, ProjQT, &
                         Matrix2, Permutation2, DOFs )
         END IF
       ELSE
          IF ( NewLinearSystem ) THEN
             Transient = ListGetString( CurrentModel % Simulation, &
                     'Simulation Type') == 'transient'

             Solver % Matrix => Matrix2
             i = Solver % MultigridLevel

             k  = ListGetInteger( Params,'Nonlinear System Max Iterations', Found )
             CALL ListAddInteger( Params,'Nonlinear System Max Iterations', 1 )

             Solver % MultigridLevel = -1

             CALL MSolverActivate( CurrentModel, Solver, Solver % dt, Transient )
             Solver % MultigridLevel = i
             CALL ListAddInteger( Params,'Nonlinear System Max Iterations',MAX(1, k) )
          END IF
       END IF
!------------------------------------------------------------------------------
!
!      Global iteration parameters:
!      ----------------------------
       MaxIter = 1

       IF ( Level == Solver % MultiGridTotal ) THEN
          MaxIter = ListGetInteger( Params,'MG Max Iterations', Found )

          IF ( .NOT. Found ) THEN
            IF( ListGetString( Params, 'Linear System Solver', Found ) == 'multigrid') THEN
              MaxIter = ListGetInteger( Params,'Linear System Max Iterations' )
            ELSE
              MaxIter = 1
            END IF
          END IF

          Tolerance = ListGetConstReal( Params,'MG Convergence Tolerance', Found )
          IF ( .NOT. Found ) THEN
             Tolerance = ListGetConstReal( Params,'Linear System Convergence Tolerance' )
          END IF
       ELSE
          MaxIter = ListGetInteger( Params,'MG Level Max Iterations', Found )
          IF ( .NOT. Found ) MaxIter = 1

          Tolerance = ListGetConstReal( Params,'MG Level Convergence Tolerance', Found )
          IF ( .NOT. Found ) Tolerance = HUGE(Tolerance)
       END IF

!------------------------------------------------------------------------------
!      Smoothing preconditiong, if not given
!      diagonal preconditioning is used:
!      -------------------------------------
       str = ListGetString( Params, 'MG Preconditioning', Found )
       IF ( .NOT. Found ) THEN
         str = ListGetString( Params,'Linear System Preconditioning', Found )
       END IF


       IF ( str == 'ilut' )  THEN

          IF ( NewLinearSystem ) THEN
             ILUTOL = ListGetConstReal( Params,'MG ILUT Tolerance', Found )
             IF ( .NOT. Found ) THEN
               ILUTOL = ListGetConstReal( Params,'Linear System ILUT Tolerance' )
             END IF
             
             IF ( Parallel ) THEN
               Condition = CRS_ILUT( PMatrix, ILUTOL )
             ELSE
               Condition = CRS_ILUT( Matrix1, ILUTOL )
             END IF
           END IF

       ELSE IF ( SEQL(str, 'ilu') ) THEN
          IF ( NewLinearSystem ) THEN
             k = ICHAR(str(4:4)) - ICHAR('0')
             IF ( k < 0 .OR. k > 9 ) k = 0
             IF ( Parallel ) THEN
                PMatrix % Cholesky = ListGetLogical( Params, &
                    'Linear System Symmetric ILU', Found )
                Condition = CRS_IncompleteLU( PMatrix, k )
             ELSE
                Matrix1 % Cholesky = ListGetLogical( Params, &
                  'Linear System Symmetric ILU', Found )
                Condition = CRS_IncompleteLU( Matrix1, k )
             END IF
          END IF

       END IF

!------------------------------------------------------------------------------
!
!      Ok, lets go:
!      ------------
       DO iter = 1,MaxIter
          ResidualNorm = GMGSweep()

          WRITE(Message,'(A,I0,A,I0,A,2E20.12E3)') 'MG Residual at level: ', &
                 Level, ' iter: ', iter,' is:', ResidualNorm/RHSNorm, ResidualNorm
          CALL Info( 'GMGSolve', Message, Level=5 )

          IF( ResidualNorm /= ResidualNorm .OR. ResidualNorm > 1.0d50 ) THEN
             CALL Fatal('GMGSolve','We seem to have diverged')
          END IF
          
          IF( Level == Solver % MultiGridTotal ) THEN
            IF ( ResidualNorm/RHSNorm < Tolerance ) EXIT
          ELSE
            IF ( ResidualNorm < Tolerance ) EXIT
          END IF
       END DO


!------------------------------------------------------------------------------
!
!      Finalize:
!      ---------
       IF ( Parallel ) THEN 
          CALL ParallelUpdateResult( Matrix1, Solution, Residual )
       END IF
       Solution(1:n) = Solution(1:n) * RHSNorm
       ForceVector(1:n) = ForceVector(1:n) * RHSNorm

       Solver % Variable => SaveVariable
       Solver % Mesh     => SaveMesh
       Solver % Matrix   => SaveMatrix

       CALL SetCurrentMesh( CurrentModel, Solver % Mesh )
       CurrentModel % Meshes => Solver % Mesh

       IF ( Level == Solver % MultiGridTotal ) THEN
          WRITE( Message,'(A,F8.2)') 'MG sweep performed in time: ', CPUTime() - tt
          CALL Info( 'GMGSolve', Message, Level=5 )
       END IF


       RETURN
!------------------------------------------------------------------------------

  CONTAINS

  
!------------------------------------------------------------------------------
    RECURSIVE FUNCTION GMGSweep() RESULT(RNorm)
!------------------------------------------------------------------------------
       INTEGER :: i,j,Rounds
       LOGICAL :: Found
       REAL(KIND=dp) :: RNorm
!------------------------------------------------------------------------------
       INTEGER :: Sweeps
       REAL(KIND=dp) :: Bnorm 
       INTEGER, POINTER :: Iters(:)
       REAL(KIND=dp), POINTER :: R1(:),R2(:)
!------------------------------------------------------------------------------

!      Presmoothing:
!      -------------

       Iters => ListGetIntegerArray( Params,'MG Sweeps',Found)
       IF(Found) THEN
         Sweeps = Iters(MIN(Level,SIZE(Iters)))
       ELSE        
         Sweeps = 1
       END IF

       WRITE( Message,'(A,I0,A,I0)') 'Performing ',Sweeps,' sweeps on Level ',Level
       CALL Info('GMGSweep',Message,Level=10)

       PSolver => Solver

       RNorm = MGSmooth( PSolver, Matrix1, Mesh1, Solution, &
            ForceVector, Residual, Level, DOFs, PreSmooth = .TRUE. )

!------------------------------------------------------------------------------
!
!     Solve (PAQ)z = Pr, x = x + Qz:
!     ==============================
!
!     Project current residual to the lower level mesh:
!     -------------------------------------------------
      CALL Info('GMGSweep','Projecting residual to lower level mesh',Level=12)
      DO i=1,DOFs
         R2 => Residual2(i:n2:DOFs)
         R1 => Residual (i:n :DOFs)

         CALL CRS_ApplyProjector( ProjPN, R1, Permutation, &
                R2, Permutation2, Trans = .FALSE. )
      END DO
      CALL Info('GMGSweep','Projector applied',Level=20)

!
!     Recursively solve (PAQ)z = Pr:
!     ------------------------------
      CALL Info('GMGSolve','Calling recursively MG solver',Level=12)
      DO i=1,Sweeps
         CALL MultigridSolve( Matrix2, Solution2, Residual2, &
             DOFs, Solver, Level-1, NewLinearSystem )
      END DO
      CALL Info('GMGSolve','Returning from MG solver call',Level=12)
 
!
!     Compute x = x + Qz:
!     -------------------

      CALL Info('GMGSweep','Projecting solution to higher level mesh',Level=12)
      DO i=1,DOFs
         R1 => Residual (i:n :DOFs)
         R2 => Solution2(i:n2:DOFs)

         CALL CRS_ApplyProjector( ProjQT, R2, Permutation2, &
                R1, Permutation, Trans = .TRUE.)
      END DO
      CALL Info('GMGSweep','Projector applied',Level=20)

      Solution(1:n) = Solution(1:n) + Residual(1:n)
!
!     Post smoothing:
!     ---------------

      RNorm = MGSmooth( PSolver, Matrix1, Mesh1, Solution, ForceVector, &
                   Residual, Level, DOFs )

!------------------------------------------------------------------------------
    END FUNCTION GMGSweep
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>     Project matrix A to B: B = PAQ^T.
!
!>     The code is a little complicated by the fact that there might be
!>     three different  numbering schemes for  the DOFs: one  for both
!>     matrices and the nodal numbering of the projectors.
!------------------------------------------------------------------------------
    SUBROUTINE ProjectMatrix( A, PermA, P, Q, B, PermB, DOFs )
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A,P,Q,B
      INTEGER :: DOFs
      INTEGER, POINTER :: PermA(:), PermB(:)
!------------------------------------------------------------------------------
      INTEGER, POINTER :: L1(:), L2(:), InvPermB(:)
      REAL(KIND=dp) :: s
      REAL(KIND=dp), POINTER :: R1(:),R2(:), bvals(:)
      INTEGER :: i,j,k,l,m,NA,NB,N,k1,k2,k3,ni,nj,RDOF,CDOF,k_a,m_q,dof
!------------------------------------------------------------------------------

      NA = A % NumberOfRows
      NB = B % NumberOfRows

!
!     Compute size for the work arrays:
!     ---------------------------------
      N = 0
      DO i=1, P % NumberOfRows
         N = MAX( N, P % Rows(i+1) - P % Rows(i) )
      END DO

      DO i=1, Q % NumberOfRows
         N = MAX( N, Q % Rows(i+1) - Q % Rows(i) )
      END DO
!
!     Allocate temporary workspace:
!     -----------------------------
      ALLOCATE( R1(N), R2(NA), L1(N), L2(N), InvPermB(SIZE(PermB)) )

!
!     Initialize:
!     -----------
      InvPermB = 0
      DO i = 1, SIZE(PermB)
         IF ( PermB(i) > 0 ) InvPermB( PermB(i) ) = i
      END DO

      R1 = 0.0d0  ! P_i:
      R2 = 0.0d0  ! Q_j: ( = Q^T_:j )

      L1 = 0      ! holds column indices of P_i in A numbering
      L2 = 0      ! holds column indices of Q_j in A numbering

      B % Values = 0.0d0

!------------------------------------------------------------------------------
!
!     Compute the projection:
!     =======================
!
!     The code below duplicated for DOFs==1, and otherwise:
!     -----------------------------------------------------
      IF ( DOFs == 1 ) THEN
!
!        Loop over rows of B:
!        --------------------
         DO i=1,NB
!
!           Get i:th row of projector P: R1=P_i
!           -----------------------------------
            ni = 0   ! number of nonzeros in P_i
            k1 = InvPermB(i)
            DO k = P % Rows(k1), P % Rows(k1+1) - 1
               l = PermA( P % Cols(k) )
               IF ( l > 0 ) THEN
                  ni = ni + 1
                  L1(ni) = l
                  R1(ni) = P % Values(k)
               END IF
            END DO

            IF ( ni <= 0 ) CYCLE
!
!           Loop over columns of row i of B:
!           --------------------------------
            DO j = B % Rows(i), B % Rows(i+1)-1
!
!              Get j:th row of projector Q: R2=Q_j
!              -----------------------------------
               nj = 0 ! number of nonzeros in Q_j
               k2 = InvPermB( B % Cols(j) )
               DO k = Q % Rows(k2), Q % Rows(k2+1)-1
                  l = PermA( Q % Cols(k) )
                  IF ( l > 0 ) THEN
                     nj = nj + 1
                     L2(nj) = l
                     R2(l)  = Q % Values(k)
                  END IF
               END DO

               IF ( nj <= 0 ) CYCLE
!
!              s=A(Q_j)^T, only entries correspoding to
!              nonzeros in P_i actually computed, then
!              B_ij = DOT( P_i, A(Q_j)^T ):
!              ------------------------------------------
               DO k=1,ni
                  k2 = L1(k)
                  s = 0.0d0
                  DO l = A % Rows(k2), A % Rows(k2+1)-1
                     s = s + R2(A % Cols(l)) * A % Values(l)
                  END DO
                  B % Values(j) = B % Values(j) + s * R1(k)
               END DO
               R2(L2(1:nj)) = 0.0_dp
            END DO
         END DO

      ELSE ! DOFs /= 1
!
!        Loop over rows of B:
!        --------------------
         DO i=1,NB/DOFs
!
!           Get i:th row of projector P: R1=P_i
!           -----------------------------------
            ni = 0   ! number of nonzeros in P_i
            k1 = InvPermB(i)
            DO k = P % Rows(k1), P % Rows(k1+1) - 1
               l = PermA( P % Cols(k) )
               IF ( l > 0 ) THEN
                  ni = ni + 1
                  L1(ni) = l
                  R1(ni) = P % Values(k)
               END IF
            END DO

            IF ( ni <= 0 ) CYCLE

            DO RDOF = 1,DOFs
!
!              Loop over columns of row i of B:
!              --------------------------------
               k1 = DOFs*(i-1) + RDOF
               DO j = B % Rows(k1), B % Rows(k1+1)-1, DOFs
!
!                 Get j:th row of projector Q: R2=Q_j
!                 -----------------------------------
                  nj = 0 ! number of nonzeros in Q_j
                  k2 = InvPermB( (B % Cols(j)-1) / DOFs + 1 )
                  DO k = Q % Rows(k2), Q % Rows(k2+1)-1
                     l = PermA( Q % Cols(k) )
                     IF ( l > 0 ) THEN
                        nj = nj + 1
                        L2(nj)  = l
                        DO CDOF=1,DOFs
                           R2(DOFs*(l-1)+CDOF) = Q % Values(k)
                        END DO
                     END IF
                  END DO

                  IF ( nj <= 0 ) CYCLE

                  DO CDOF=0,DOFs-1
!
!                    s = A(Q_j)^T, only entries correspoding to
!                    nonzeros in P_i actually  computed, then
!                    B_ij = DOT( P_i, A(Q_j)^T ):
!                    ------------------------------------------
                     DO k=1,ni
                        k2 = DOFs * (L1(k)-1) + RDOF
                        s = 0.0d0
                        DO l = A % Rows(k2)+CDOF, A % Rows(k2+1)-1, DOFs
                           s = s + R2(A % Cols(l)) * A % Values(l)
                        END DO
                        B % Values(j+CDOF) = B % Values(j+CDOF) + s * R1(k)
                     END DO
                     R2(DOFs*(L2(1:nj)-1)+CDOF+1) = 0.0d0
                  END DO
               END DO
            END DO
         END DO
      END IF

      DEALLOCATE( R1, R2, L1, L2, InvPermB )
!------------------------------------------------------------------------------
    END SUBROUTINE ProjectMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION MGnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       REAL(KIND=dp) :: s
       REAL(KIND=dp) CONTIG :: x(:)
!------------------------------------------------------------------------------
       IF ( .NOT. Parallel ) THEN
          s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
       ELSE
          s = ParallelNorm( n, x )
       END IF
!------------------------------------------------------------------------------
    END FUNCTION MGnorm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE MGmv( A, x, b, Update )
!------------------------------------------------------------------------------
       REAL(KIND=dp) CONTIG :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       IF ( .NOT. Parallel ) THEN
         CALL CRS_MatrixVectorMultiply( A, x, b )
       ELSE
         IF ( PRESENT( Update ) ) THEN
           CALL ParallelMatrixVector( A,x,b,Update,ZeroNotOwned=.TRUE. )
         ELSE
           CALL ParallelMatrixVector( A,x,b,ZeroNotOwned=.TRUE. )
         END IF
       END IF
!------------------------------------------------------------------------------
    END SUBROUTINE MGmv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION MGCnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       REAL(KIND=dp) :: s
       COMPLEX(KIND=dp) CONTIG :: x(:)
!------------------------------------------------------------------------------
       s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
!------------------------------------------------------------------------------
    END FUNCTION MGCnorm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE MGCmv( A, x, b, Update )
!------------------------------------------------------------------------------
       COMPLEX(KIND=dp) CONTIG :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       CALL CRS_ComplexMatrixVectorMultiply( A, x, b )
!------------------------------------------------------------------------------
    END SUBROUTINE MGCmv
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE GMGSolve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Multigrid solution in the case when different levels are different 
!> power of element basis functions. 
!------------------------------------------------------------------------------
    RECURSIVE SUBROUTINE PMGSolve( Matrix1, Solution, &
        ForceVector, DOFs, Solver, Level, NewSystem )
!------------------------------------------------------------------------------
       USE ModelDescription
       IMPLICIT NONE

       TYPE(Matrix_t), POINTER :: Matrix1
       INTEGER :: DOFs, Level
       LOGICAL, OPTIONAL :: NewSystem
       TYPE(Solver_t), TARGET :: Solver       
       REAL(KIND=dp), TARGET CONTIG :: ForceVector(:), Solution(:)
!------------------------------------------------------------------------------
       TYPE(Variable_t), POINTER :: Variable1, TimeVar, SaveVariable
       TYPE(Mesh_t), POINTER   :: Mesh1, Mesh2, SaveMesh
       TYPE(Matrix_t), POINTER :: Matrix2, PMatrix, SaveMatrix
       TYPE(Solver_t), POINTER :: PSolver             

       INTEGER :: i,j,k,l,m,n,n2,k1,k2,iter,MaxIter = 100, RDOF, CDOF,ndofs
       LOGICAL :: Condition, Found, Parallel, Project,Transient
       CHARACTER(LEN=MAX_NAME_LEN) :: Path,str,mgname, LowestSolver

       TYPE(Matrix_t), POINTER :: ProjPN, ProjQT
       INTEGER, POINTER :: Permutation(:), Permutation2(:), Indexes(:), Deg(:), &
                PatLevel(:), iPerm(:),dum(:)

       REAL(KIND=dp) :: ResidualNorm, RHSNorm, Tolerance, ILUTOL, tmp
       REAL(KIND=dp), POINTER CONTIG :: Residual(:), Residual2(:), Solution2(:), &
               SaveVals(:), Basis(:)

       REAL(KIND=dp), POINTER, SAVE :: SolutionStore(:,:), Degree(:)

#ifdef USE_ISO_C_BINDINGS
       REAL(KIND=dp) :: tt, detJ
#else
       REAL(KIND=dp) :: CPUTime, tt, detJ
#endif

       LOGICAL :: NewLinearSystem, stat, LIter

       SAVE NewLinearSystem

       TYPE(Nodes_t), SAVE :: Nodes
       TYPE(Element_t), POINTER :: Element
       TYPE(ValueList_t), POINTER :: Params
!------------------------------------------------------------------------------
       tt = CPUTime()

!
!      Initialize:
!      -----------
       Parallel = ParEnv % PEs > 1
       Params => Solver % Values

       IF ( Level == Solver % MultiGridLevel ) THEN
          NewLinearSystem = .TRUE.
          IF ( PRESENT( NewSystem ) ) THEN
             NewLinearSystem = NewLinearSystem .AND. NewSystem
          END IF
       END IF

!---------------------------------------------------------------------
!
!      If at lowest level, solve directly:
!      -----------------------------------
       IF ( Level <= 1 ) THEN

          CALL ListPushNamespace('mglowest:')

          CALL ListAddLogical( Params,'mglowest: Linear System Free Factorization', .FALSE. )
          IF ( NewLinearSystem ) THEN
            CALL ListAddLogical( Params,'mglowest: Linear System Refactorize', .TRUE. )
          ELSE
            CALL ListAddLogical( Params,'mglowest: Linear System Refactorize', .FALSE. )
          END IF
          NewLinearSystem = .FALSE.

          LowestSolver = ListGetString(Params,'MG Lowest Linear Solver',Found)
          IF ( .NOT. Found ) THEN
            LIter=ListGetLogical(Params,'MG Lowest Linear Solver Iterative',Found)
            IF ( .NOT. Found .AND. Parallel ) LIter=.TRUE.
            LowestSolver='direct'
            IF ( LIter ) LowestSolver='iterative'
          END IF

          CALL Info('PMGSolve','Starting lowest linear solver using: '//TRIM(LowestSolver),Level=10 )

          SELECT CASE(LowestSolver)

          CASE('iterative')
            IF ( Parallel ) THEN
              CALL ParallelIter( Matrix1, Matrix1 % ParallelInfo, DOFs, &
                  Solution, ForceVector, Solver, Matrix1 % ParMatrix )
            ELSE
              CALL IterSolver( Matrix1, Solution, ForceVector, Solver )
            END IF

          CASE('direct')
            CALL DirectSolver( Matrix1, Solution, ForceVector, Solver )

          CASE('smoother')
            IF ( Parallel ) THEN
               CALL ParallelInitSolve( Matrix1, Solution, &
                       ForceVector, Residual, NewLinearSystem )
            END IF
            PSolver => Solver
            tmp = MGSmooth( PSolver, Matrix1, Mesh1, Solution, &
                 ForceVector, Residual, Level, DOFs, LowestSmooth = .TRUE. )

          CASE('none') 
            CALL Info('GMGSolve','Applying no solver for coarsest level')

          CASE DEFAULT
             CALL Warn( 'GMGSolve', 'Unknown solver selection for MG lowest level' )
             CALL Warn( 'GMGSolve', 'Using iterative solver' )

             IF ( Parallel ) THEN
               CALL ParallelIter( Matrix1, Matrix1 % ParallelInfo, DOFs, &
                   Solution, ForceVector, Solver, Matrix1 % ParMatrix )
             ELSE
               CALL IterSolver( Matrix1, Solution, ForceVector, Solver )
             END IF
          END SELECT

          CALL ListPopNamespace('mglowest:')

          RETURN
       END IF
!
       n = Matrix1 % NumberOfRows
       ALLOCATE( Residual(n) )
       Residual = 0.0_dp
!
!      Parallel initializations:
!      -------------------------
       IF ( Parallel ) THEN
          CALL ParallelInitSolve( Matrix1, Solution, &
              ForceVector, Residual, NewLinearSystem )
          PMatrix => ParallelMatrix(Matrix1) 
       END IF
!
!      Compute residual:
!      -----------------
       CALL MGmv(Matrix1,Solution, Residual,.TRUE.)
       Residual(1:n) = ForceVector(1:n) - Residual(1:n)
 
       RHSNorm = MGnorm(n, ForceVector)
       ResidualNorm = MGnorm( n, Residual ) / RHSNorm
 
       Tolerance = ListGetConstReal( Params,'Linear System Convergence Tolerance' )

!---------------------------------------------------------------------
!
!      Initialize the multilevel solve:
!      --------------------------------
       SaveMatrix   => Solver % Matrix
       SaveVariable => Solver % Variable

       Variable1   => Solver % Variable
       Permutation => Variable1 % Perm

       Matrix2 => Matrix1 % Child
      
       IF ( NewLinearSystem ) THEN
         IF ( .NOT. ASSOCIATED(Matrix2) ) THEN

           n2 = Solver % Mesh % MaxElementDOFs
           ALLOCATE(Degree(n), Indexes(n2), Deg(n2))

           DO i=1,Solver % NumberOfActiveElements
             Element => Solver % Mesh % Elements(Solver % ActiveElements(i))
             n = PMGGetElementDOFs( Indexes, Element )
             CALL ElementBasisDegree(Element, Deg)

             DO j=1,n
               DO k=1,DOFs
                 Degree(DOFs*(Permutation(Indexes(j))-1)+k) = Deg(j)
               END DO
             END DO
           END DO
           DEALLOCATE(Indexes,Deg)
            
           PatLevel => ListGetIntegerArray( Params,'MG P at Level',Found)
           IF(.NOT.Found) THEN
             ALLOCATE(PatLevel(Level-1))
             PatLevel = (/(i,i=1,level-1)/)
           END IF

           n = Matrix1 % NumberOfRows
           ALLOCATE(Matrix1 % Eperm(n))
           Matrix1 % Eperm=(/(i,i=1,n)/)

           PMatrix => Matrix1
           DO l=level-1,1,-1
             Matrix2 => AllocateMatrix()
             Matrix2 % Parent => Pmatrix
             Pmatrix % Child  => Matrix2
             n2 = COUNT(Degree<=PatLevel(l))
             Matrix2 % NumberOfRows = n2
             ALLOCATE( Matrix2 % Rows(n2+1), Matrix2 % Diag(n2), Matrix2 % Eperm(n) )

             Matrix2 % Eperm = 0
             j = 0
             DO i=1,n
               IF ( Degree(i) <= PatLevel(l) ) THEN
                 j = j + 1
                 Matrix2 % Eperm(i) = j
               END IF
             END DO

             Matrix2 % Rows(1) = 1
             j = 0
             DO i=1,n
               IF ( Matrix2 % Eperm(i) == 0 ) CYCLE
               j = j + 1
               Matrix2 % Rows(j+1) = Matrix2 % Rows(j)
               DO k=Matrix1 % Rows(i),Matrix1 % Rows(i+1)-1
                 IF ( Matrix2 % Eperm(Matrix1 % Cols(k)) /= 0 ) &
                   Matrix2 % Rows(j+1) = Matrix2 % Rows(j+1)+1
               END DO
             END DO

             ALLOCATE( Matrix2 % Cols(Matrix2 % Rows(n2+1)), &
                       Matrix2 % Values(Matrix2 % Rows(n2+1)) )

             j = 0
             m = 0
             DO i=1,n
               IF ( Matrix2 % Eperm(i) == 0 ) CYCLE
               j = j + 1
               DO k=Matrix1 % Rows(i),Matrix1 % Rows(i+1)-1
                 IF ( Matrix2 % Eperm(Matrix1 % Cols(k)) /= 0 ) THEN
                   m = m + 1
                   Matrix2 % Cols(m) = Matrix2 % Eperm(Matrix1 % Cols(k))
                   IF ( j==Matrix2 % Cols(m)) Matrix2 % Diag(j)=m
                 END IF
               END DO
             END DO
             Pmatrix => Matrix2
           END DO
           Matrix2 => Matrix1 % Child
           DEALLOCATE(Degree)
         END IF

         IF (Level==Solver % MultigridLevel) THEN
           Pmatrix => Matrix1 % Child
           DO l=Level-1,1,-1
             j = 0
             m = 0
             DO i=1,n
               IF ( Pmatrix % Eperm(i) == 0 ) CYCLE
               j = j + 1
               DO k=Matrix1 % Rows(i),Matrix1 % Rows(i+1)-1
                 IF ( Pmatrix % Eperm(Matrix1 % Cols(k)) /= 0 ) THEN
                   m = m + 1
                   Pmatrix % Values(m) = Matrix1 % Values(k)
                 END IF
               END DO
             END DO
             Pmatrix => Pmatrix % Child
           END DO
         END IF
       END IF


!------------------------------------------------------------------------------
!
!      Some more initializations:
!      --------------------------
       n  = Matrix1 % NumberOfRows
       n2 = Matrix2 % NumberOfRows

       ALLOCATE( Residual2(n2), Solution2(n2) )
       Residual2 = 0._dp
       Solution2 = 0._dp

       IF ( Parallel ) THEN
         IF(.NOT.ASSOCIATED(Matrix2 % ParMatrix)) THEN
           k = SIZE(Solver % Variable % Perm)
           ALLOCATE(Permutation(k),iPerm(k))
           Permutation=0; iPerm=0
           DO i=1,k
             j=Solver % Variable % Perm(i)
             IF(j<=0) CYCLE
             iPerm(j) = i
           END DO

           j=0
           DO i=1,n
             IF(Matrix2 % Eperm(i)>0) THEN
               j = j + 1
               k=(i-1)/DOFs+1
               Permutation(iPerm(k))=(j-1)/DOFs+1
             END IF
           END DO
           Matrix2 % Comm = Matrix1 % Comm
           CALL ParallelInitMatrix(Solver,Matrix2,Permutation)
           DEALLOCATE(Permutation, iPerm)
         END IF
       END IF

!------------------------------------------------------------------------------
!
!      Global iteration parameters:
!      ----------------------------
       MaxIter = 1

       IF ( Level == Solver % MultiGridTotal ) THEN
          MaxIter = ListGetInteger(Params,'MG Max Iterations',Found)

          IF ( .NOT. Found ) THEN
            IF( ListGetString( Params,'Linear System Solver', Found ) == 'multigrid') THEN
              MaxIter = ListGetInteger(Params,'Linear System Max Iterations')
            ELSE
              MaxIter = 1
            END IF
          END IF
          Tolerance = ListGetConstReal( Params,'MG Convergence Tolerance', Found )
          IF ( .NOT. Found ) THEN
            Tolerance = ListGetConstReal( Params,'Linear System Convergence Tolerance' )
          END IF
       ELSE
          MaxIter = ListGetInteger( Params,'MG Level Max Iterations', Found )
          IF ( .NOT. Found ) MaxIter = 1

          Tolerance = ListGetConstReal( Params,'MG Level Convergence Tolerance', Found )
          IF ( .NOT. Found ) Tolerance = HUGE(Tolerance)
       END IF

!------------------------------------------------------------------------------
!      Smoothing preconditiong, if not given
!      diagonal preconditioning is used:
!      -------------------------------------
       IF ( NewLinearSystem ) THEN
         str = ListGetString( Params, 'MG Preconditioning', Found )
         IF ( .NOT. Found ) THEN
           str = ListGetString( Params,'Linear System Preconditioning', Found )
         END IF
         IF ( str == 'ilut' )  THEN
           ILUTOL = ListGetConstReal( Params,'MG ILUT Tolerance', Found )
           IF ( .NOT. Found ) THEN
             ILUTOL = ListGetConstReal( Params,'Linear System ILUT Tolerance' )
           END IF
           IF ( Parallel ) THEN
             Condition = CRS_ILUT( PMatrix, ILUTOL )
           ELSE
             Condition = CRS_ILUT( Matrix1, ILUTOL )
           END IF
         ELSE IF ( SEQL(str, 'ilu') ) THEN
           k = ICHAR(str(4:4)) - ICHAR('0')
           IF ( k < 0 .OR. k > 9 ) k = 0
           IF ( Parallel ) THEN
              PMatrix % Cholesky = ListGetLogical( Params, &
                  'Linear System Symmetric ILU', Found )
              Condition = CRS_IncompleteLU( PMatrix, k )
            ELSE
              Matrix1 % Cholesky = ListGetLogical( Params, &
                  'Linear System Symmetric ILU', Found )
              Condition = CRS_IncompleteLU( Matrix1, k )
            END IF
          END IF
        END IF


!------------------------------------------------------------------------------
!
!      Ok, lets go:
!      ------------
       DO iter = 1,MaxIter
          ResidualNorm = PMGSweep()

          WRITE(Message,'(A,I0,A,I0,A,2E20.12E3)') 'MG Residual at level: ', &
                 Level, ' iter: ', iter,' is:', ResidualNorm/RHSNorm, ResidualNorm
          CALL Info( 'PMGSolve', Message, Level=5 )


          IF( ResidualNorm /= ResidualNorm .OR. ResidualNorm > 1.0d50 ) THEN
             CALL Fatal('PMGSolve','We seem to have diverged')
          END IF

          IF( Level == Solver % MultiGridTotal ) THEN
            IF ( ResidualNorm/RHSNorm < Tolerance ) EXIT
          ELSE
            IF ( ResidualNorm < Tolerance ) EXIT
          END IF
       END DO

!------------------------------------------------------------------------------
!
!      Finalize:
!      ---------
       IF ( Parallel ) THEN 
          CALL ParallelUpdateResult( Matrix1, Solution, Residual )
       END IF
       Solver % Matrix   => SaveMatrix
       Solver % Variable => SaveVariable

       DEALLOCATE( Residual, Residual2, Solution2 )

       IF ( Level == Solver % MultiGridTotal ) THEN
          WRITE( Message,'(A,F8.2)') 'MG iter time: ', CPUTime() - tt
          CALL Info( 'PMGSolve', Message, Level=5 )
       END IF


       RETURN
!------------------------------------------------------------------------------

  CONTAINS

!------------------------------------------------------------------------------
  FUNCTION PMGGetElementDOFs( Indexes, UElement, USolver )  RESULT(NB)
!------------------------------------------------------------------------------
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     INTEGER :: Indexes(:)

     TYPE(Solver_t),  POINTER :: Solver
     TYPE(Element_t), POINTER :: Element, Parent

     LOGICAL :: Found, GB
     INTEGER :: nb,i,j,EDOFs, FDOFs, BDOFs,FaceDOFs, EdgeDOFs, BubbleDOFs

     Element => UElement

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     NB = 0

     IF ( ListGetLogical( Params, 'Discontinuous Galerkin', Found ) ) THEN
        DO i=1,Element % DGDOFs
           NB = NB + 1
           Indexes(NB) = Element % DGIndexes(i)
        END DO

        IF ( ASSOCIATED( Element % BoundaryInfo ) ) THEN
           IF ( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
              DO i=1,Element % BoundaryInfo % Left % DGDOFs
                 NB = NB + 1
                 Indexes(NB) = Element % BoundaryInfo % Left % DGIndexes(i)
              END DO
           END IF
           IF ( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
              DO i=1,Element % BoundaryInfo % Right % DGDOFs
                 NB = NB + 1
                 Indexes(NB) = Element % BoundaryInfo % Right % DGIndexes(i)
              END DO
           END IF
        END IF

        IF ( NB > 0 ) RETURN
     END IF

     DO i=1,Element % NDOFs
        NB = NB + 1
        Indexes(NB) = Element % NodeIndexes(i)
     END DO

     FaceDOFs   = Solver % Mesh % MaxFaceDOFs
     EdgeDOFs   = Solver % Mesh % MaxEdgeDOFs
     BubbleDOFs = Solver % Mesh % MaxBDOFs

     IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
        DO j=1,Element % TYPE % NumberOFEdges
          EDOFs = Solver % Mesh % Edges( Element % EdgeIndexes(j) ) % BDOFs
          DO i=1,EDOFs
             NB = NB + 1
             Indexes(NB) = EdgeDOFs*(Element % EdgeIndexes(j)-1) + &
                      i + Solver % Mesh % NumberOfNodes
          END DO
        END DO
     END IF

     IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
        DO j=1,Element % TYPE % NumberOFFaces
           FDOFs = Solver % Mesh % Faces( Element % FaceIndexes(j) ) % BDOFs
           DO i=1,FDOFs
              NB = NB + 1
              Indexes(NB) = FaceDOFs*(Element % FaceIndexes(j)-1) + i + &
                 Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges
           END DO
        END DO
     END IF

     GB = ListGetLogical( Params, 'Bubbles in Global System', Found )
     IF (.NOT.Found) GB = .TRUE.

     IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
       IF (.NOT. isActivePElement(Element) ) RETURN

       Parent => Element % BoundaryInfo % Left
       IF (.NOT.ASSOCIATED(Parent) ) &
         Parent => Element % BoundaryInfo % Right
       IF (.NOT.ASSOCIATED(Parent) ) RETURN

       IF ( ASSOCIATED( Parent % EdgeIndexes ) ) THEN
         EDOFs = Element % BDOFs
         DO i=1,EDOFs
           NB = NB + 1
           Indexes(NB) = EdgeDOFs*(Parent % EdgeIndexes(Element % PDefs % LocalNumber)-1) + &
                    i + Solver % Mesh % NumberOfNodes
         END DO
       END IF

       IF ( ASSOCIATED( Parent % FaceIndexes ) ) THEN
         FDOFs = Element % BDOFs
         DO i=1,FDOFs
           NB = NB + 1
           Indexes(NB) = FaceDOFs*(Parent % FaceIndexes(Element % PDefs % LocalNumber)-1) + i + &
              Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges
         END DO
       END IF
     ELSE IF ( GB ) THEN
        IF ( ASSOCIATED( Element % BubbleIndexes ) ) THEN
           DO i=1,Element % BDOFs
              NB = NB + 1
              Indexes(NB) = FaceDOFs*Solver % Mesh % NumberOfFaces + &
                 Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges + &
                   Element % BubbleIndexes(i)
           END DO
        END IF
     END IF
!------------------------------------------------------------------------------
  END FUNCTION PMGGetElementDOFs
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
    RECURSIVE FUNCTION PMGSweep() RESULT(RNorm)
!------------------------------------------------------------------------------
       INTEGER :: i,j,Rounds
       LOGICAL :: Found
       REAL(KIND=dp) :: RNorm
!------------------------------------------------------------------------------

!      Presmoothing:
!      -------------

       PSolver => Solver
       RNorm = MGSmooth( PSolver, Matrix1, Solver % Mesh, Solution, &
            ForceVector, Residual, Level, DOFs, PreSmooth = .TRUE. )
!
!------------------------------------------------------------------------------

!
!     Recursively solve (PAQ)z = Pr:
!     ------------------------------
      DO i=1,SIZE(Matrix2 % Eperm)
        j = Matrix2 % Eperm(i)
        IF ( j>0 )  THEN
          k = Matrix1 % Eperm(i)
          Residual2(j) = Residual(k)
        END IF
      END DO

      CALL MultigridSolve( Matrix2, Solution2, Residual2, &
           DOFs, Solver, Level-1, NewLinearSystem )

      DO i=1,SIZE(Matrix2 % Eperm)
        j = Matrix2 % Eperm(i)
        IF (j>0)  THEN
          k = Matrix1 % Eperm(i)
          Solution(k) = Solution(k) + Solution2(j)
        END IF
      END DO
!
!     Post smoothing:
!     ---------------
      RNorm = MGSmooth( PSolver, Matrix1, Solver % Mesh, Solution, &
             ForceVector, Residual, Level, DOFs )
!------------------------------------------------------------------------------
    END FUNCTION PMGSweep
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      FUNCTION MGnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
        INTEGER :: n
        REAL(KIND=dp)  :: s
        REAL(KIND=dp) CONTIG :: x(:)
!------------------------------------------------------------------------------
        IF ( .NOT. Parallel ) THEN
          s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
        ELSE
          s = ParallelNorm(n, x)
        END IF
!------------------------------------------------------------------------------
      END FUNCTION MGnorm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE MGmv( A, x, b, Update )
!------------------------------------------------------------------------------
        REAL(KIND=dp) CONTIG :: x(:), b(:)
        TYPE(Matrix_t), POINTER :: A
        LOGICAL, OPTIONAL :: Update
!------------------------------------------------------------------------------
        IF ( .NOT. Parallel ) THEN
          CALL CRS_MatrixVectorMultiply( A, x, b )
        ELSE
          IF ( PRESENT( Update ) ) THEN
            CALL ParallelMatrixVector( A,x,b,Update,ZeroNotOwned=.TRUE. )
          ELSE
            CALL ParallelMatrixVector(A,x,b,ZeroNotOwned=.TRUE.)
          END IF
        END IF
!------------------------------------------------------------------------------
      END SUBROUTINE MGmv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  END SUBROUTINE PMGSolve
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Subroutine containing algebraic multigrid solver using roughly the standard Ruge-Stuben interpolation.
!> Also some ideas of compatible relaxation have been tested within the context but their usability has
!> so far been rather limited.
!
!        Author: Peter R�back
!        Modified by: 
!        Date of modification: 30.10.2003
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE AMGSolve( Matrix1, Solution, &
    ForceVector, DOFs, Solver, Level, NewSystem )
!------------------------------------------------------------------------------
    USE ModelDescription
    IMPLICIT NONE
    
    TYPE(Matrix_t), POINTER :: Matrix1
    INTEGER :: DOFs, Level
    LOGICAL, OPTIONAL :: NewSystem
    TYPE(Solver_t), TARGET :: Solver
    REAL(KIND=dp) CONTIG :: ForceVector(:), Solution(:)
!------------------------------------------------------------------------------
    TYPE AMG_t
      INTEGER, POINTER :: CF(:)
      INTEGER, POINTER :: InvCF(:) 
    END TYPE AMG_t

    TYPE(AMG_t), POINTER :: AMG(:)
    TYPE(Mesh_t), POINTER   :: Mesh
    TYPE(Matrix_t), POINTER :: Matrix2, Pmatrix
    TYPE(Matrix_t), POINTER :: ProjPN, ProjQT, ProjT 
    TYPE(Solver_t), POINTER :: PSolver
   
    INTEGER :: i,j,k,l,m,n,n2,k1,k2,iter,MaxIter = 100, DirectLimit, &
        MinLevel, InvLevel
    LOGICAL :: Condition, Found, Parallel, EliminateDir, CoarseSave, RecomputeProjector
    CHARACTER(LEN=MAX_NAME_LEN) :: str,IterMethod,FileName
    INTEGER, POINTER :: CF(:), InvCF(:)
    LOGICAL, POINTER :: Fixed(:)
    
    REAL(KIND=dp), ALLOCATABLE, TARGET :: Residual(:), Solution2(:), Work2(:)
    REAL(KIND=dp), POINTER CONTIG :: Residual2(:)
    REAL(KIND=dp) :: ResidualNorm, RHSNorm, Tolerance, ILUTOL
#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: tt
#else
    REAL(KIND=dp) :: CPUTime, tt
#endif
    TYPE(ValueList_t), POINTER :: Params

    LOGICAL :: NewLinearSystem, gotit

    SAVE NewLinearSystem, AMG, MinLevel
    
!------------------------------------------------------------------------------

    IF( ParEnv % PEs > 1 ) THEN
      CALL Fatal('AMGSolve','This algebraic multigrid is not parallel')
    END IF

    IF(.FALSE.) THEN 
      WRITE(Message,*) 'Starting level ',Level,NewLinearSystem,NewSystem
      CALL Info('AMGSolve',Message)
    END IF

    Mesh => Solver % Mesh    
    Params => Solver % Values

    tt = CPUTime()
!
!      Initialize:
!      -----------
    Parallel = ParEnv % PEs > 1
    InvLevel = 1 + Solver % MultiGridTotal - Level

    ! This is a counter that for the first full resursive round keeps the 
    ! flag NewLinearSystem true.
    IF ( Level == Solver % MultiGridLevel ) THEN
      IF( ListGetLogical(Params,'MG Recompute Projector',GotIt) ) THEN
        NewLinearSystem = .TRUE.
      ELSE 
        NewLinearSystem = .NOT. ASSOCIATED(Matrix1 % Parent)
      END IF
      MinLevel = Solver % MultiGridLevel
    END IF

!---------------------------------------------------------------------
!
!      If at lowest level, solve directly:
!      -----------------------------------
    n = Matrix1 % NumberOfRows

    DirectLimit = ListGetInteger(Params,'MG Lowest Linear Solver Limit',GotIt) 
    IF(.NOT. GotIt) DirectLimit = 20

    IF ( Level <= 1 .OR. n < DirectLimit) THEN
      CALL ListPushNamespace('mglowest:')

      NewLinearSystem = .FALSE.
      IF(PRESENT(NewSystem)) NewSystem = .FALSE.
      IF ( ListGetLogical( Params, 'MG Lowest Linear Solver Unsolve',gotit ) ) THEN
        CALL Info('AMGSolve','Leaving lowest level of AMG cycle unsolved',Level=12)
      ELSE IF ( .NOT. Parallel ) THEN
        IF ( ListGetLogical( Params, 'MG Lowest Linear Solver Iterative',gotit ) ) THEN
          CALL IterSolver( Matrix1, Solution, ForceVector, Solver )
        ELSE
          CALL DirectSolver( Matrix1, Solution, ForceVector, Solver )
        END IF
      ELSE
        CALL ParallelIter( Matrix1, Matrix1 % ParallelInfo, DOFs, &
            Solution, ForceVector, Solver, Matrix1 % ParMatrix )
      END IF

      CALL ListPopNamespace('mglowest:')

      RETURN
    END IF

    n = Matrix1 % NumberOfRows
    ALLOCATE( Residual(n) )
    Residual = 0.0d0

!      Parallel initializations:
!      -------------------------
    IF ( Parallel ) THEN
      CALL ParallelInitSolve( Matrix1, Solution, ForceVector, Residual )
      PMatrix => ParallelMatrix( Matrix1 ) 
    END IF

!      Compute residual:
!      -----------------
    IF(SIZE(Solution) /= Matrix1 % NumberOfRows) THEN
      CALL WARN('AMGSolve','Solution and matrix sizes differ')
    END IF

    CALL MGmv( Matrix1, Solution, Residual, .TRUE. )
    Residual(1:n) = ForceVector(1:n) - Residual(1:n)

    RHSNorm = MGnorm( n, ForceVector )
    ResidualNorm = MGnorm( n, Residual ) / RHSNorm

    Tolerance = ListGetConstReal( Params,'Linear System Convergence Tolerance' )

    CoarseSave = ListGetLogical(Params,'MG Coarse Level Save',GotIt)

!---------------------------------------------------------------------
!
!      Initialize the multilevel solve:
!      --------------------------------


!---------------------------------------------------------------------
!      Create the Projectors between different levels
!      ---------------------------------------------

    IF ( NewLinearSystem ) THEN

      ! If the projection matrix is made again deallocate the old projectors
      ! Note that this should be done recursively!
      IF(ASSOCIATED(Matrix1 % Parent)) CALL FreeMatrix(Matrix1 % Parent)
      IF(ASSOCIATED(Matrix1 % Ematrix)) CALL FreeMatrix(Matrix1 % Ematrix)

      IF( CoarseSave .AND. Level == Solver % MultiGridTotal) THEN
        ALLOCATE(AMG(Solver % MultiGridTotal))
      END IF

      CALL Info('AMGSolve','------------------------------------------------')
      WRITE( Message, '(A,I3)' ) 'Creating a new matrix and projector for level',Level
      CALL Info('AMGSolve', Message)
      MinLevel = MIN(MinLevel,Level)
      
      EliminateDir = ListGetLogical(Params,'MG Eliminate Dirichlet',GotIt) 
      IF(.NOT. GotIt) EliminateDir = .TRUE.
      IF(Level /= Solver % MultiGridTotal) EliminateDir = .FALSE.

      ! Determine projector with one dof and use it for others
      ! Determine C/F split and projector from one dof, then project everything with it
      CALL ChooseCoarseNodes(Matrix1, Solver, ProjT, DOFs, CF, InvCF)        
      
      IF( ListGetLogical(Params,'MG Projector Matrix Save',GotIt) ) THEN
        WRITE(Filename,'(A,I1,A)') 'P',Level,'.dat'
        CALL SaveMatrix(ProjT,TRIM(FileName))          
      END IF
      
      ! The initial projection matrix is a transpose of the wanted one. 
      ! It is needed only in determining the matrix structure.       
      ProjPN => CRS_Transpose( ProjT )
      !     CALL ParallelInitMatrix( Solver, Solver % Matrix )
      CALL CRS_ProjectMatrixCreate( Matrix1, ProjPN, ProjT, Matrix2, DOFs)
      CALL FreeMatrix(ProjT)             

      Matrix2 % Child  => Matrix1
      Matrix1 % Parent => Matrix2
      ProjQT => ProjPN      
      Matrix1 % Ematrix => ProjPN

      WRITE( Message, '(A,F8.2)' ) 'MG coarse matrix creation time: ', CPUTime() - tt
      CALL Info( 'AMGSolve', Message, Level=5 )

      ! Make the cnodes point to the set of original nodes
      IF(CoarseSave) THEN
        ProjPN % Perm => CF
        ! Make the cnodes point to the set of original nodes         
        IF(Level < Solver % MultiGridTotal) THEN
          IF( ASSOCIATED(InvCF)) THEN
            AMG(Level) % InvCF => InvCF
            DO i=1,SIZE(InvCF)
              IF(InvCF(i) > 0) InvCF(i) = AMG(Level+1) % InvCF(InvCF(i))
            END DO
          END IF
        END IF
      END IF

      IF( ListGetLogical(Params,'MG Projected Matrix Save', GotIt ) ) THEN
        WRITE(Filename,'(A,I1,A)') 'B',Level,'.dat'
        CALL SaveMatrix(Matrix2,TRIM(FileName))          
      END IF
            
      IF( CoarseSave ) THEN
        CALL AMGTest(0)
        !        CALL AMGTest(1)
        !        CALL AMGTest(2)         
      END IF

    ELSE
      ! .NOT. new linear system
      Matrix2 => Matrix1 % Parent      
      ProjPN => Matrix1 % Ematrix
      CF => ProjPN % Perm
      ProjQT => ProjPN      
    END IF  
 
    n  = Matrix1 % NumberOfRows
    n2 = Matrix2 % NumberOfRows
    Residual2 => Matrix2 % RHS
    ALLOCATE( Work2(n2), Solution2(n2) )

!------------------------------------------------------------------------------
!      Global iteration parameters:
!      ----------------------------

    MaxIter = 1
    IF ( Level == Solver % MultiGridTotal ) THEN
      MaxIter = ListGetInteger( Params,'MG Max Iterations', Found )
      IF ( .NOT. Found ) THEN
        IF( ListGetString( Params, 'Linear System Solver', Found ) == 'multigrid') THEN
          MaxIter = ListGetInteger( Params,'Linear System Max Iterations' )
        ELSE
          MaxIter = 1
        END IF
      END IF
      Tolerance = ListGetConstReal( Params,'MG Convergence Tolerance', Found )
      IF ( .NOT. Found ) THEN
        Tolerance = ListGetConstReal( Params,'Linear System Convergence Tolerance' )
      END IF
    ELSE
      MaxIter = ListGetInteger( Params,'MG Level Max Iterations', Found )
      IF ( .NOT. Found ) MaxIter = 1         
      Tolerance = ListGetConstReal( Params,'MG Level Convergence Tolerance', Found )
      IF ( .NOT. Found ) Tolerance = HUGE(Tolerance)
    END IF
   
!   Smoothing preconditiong, if not given diagonal preconditioning is used:
!   ----------------------------------------------------------------------
    str = ListGetString( Params, 'MG Preconditioning', Found )
    IF ( .NOT. Found ) THEN
      str = ListGetString( Params,'Linear System Preconditioning', Found )
    END IF
   
    IF ( str == 'ilut' )  THEN
      IF ( NewLinearSystem ) THEN
        ILUTOL = ListGetConstReal( Params,'MG ILUT Tolerance', GotIt )
        IF ( .NOT. GotIt ) THEN
          ILUTOL = ListGetConstReal( Params,'Linear System ILUT Tolerance' )
        END IF
        
        IF ( Parallel ) THEN
          Condition = CRS_ILUT( PMatrix, ILUTOL )
        ELSE
          Condition = CRS_ILUT( Matrix1, ILUTOL )
        END IF
      END IF
      
    ELSE IF ( SEQL(str, 'ilu') ) THEN      
      IF ( NewLinearSystem ) THEN
        k = ICHAR(str(4:4)) - ICHAR('0')
        IF ( k < 0 .OR. k > 9 ) k = 0
        IF ( Parallel ) THEN
          PMatrix % Cholesky = ListGetLogical( Params, &
                 'Linear System Symmetric ILU', Found )
          Condition = CRS_IncompleteLU( PMatrix, k )
        ELSE
          Matrix1 % Cholesky = ListGetLogical( Params, &
                 'Linear System Symmetric ILU', Found )
          Condition = CRS_IncompleteLU( Matrix1, k )
        END IF
      END IF      
    END IF

!------------------------------------------------------------------------------
!      Ok, lets go:
!      ------------
    DO iter = 1,MaxIter
      ResidualNorm = AMGSweep() / RHSNorm
    
      WRITE(Message,'(A,I0,A,I0,A,E20.12E3)') 'MG Residual at level: ', &
          Level, ' iter: ', iter,' is:', ResidualNorm
      CALL Info( 'AMGSolve', Message, Level=5 )

      IF( ResidualNorm /= ResidualNorm .OR. ResidualNorm > 1.0d50 ) THEN
         CALL Fatal('AMGSolve','We seem to have diverged')
      END IF
            
      IF ( ResidualNorm < Tolerance ) EXIT

    END DO
    
!------------------------------------------------------------------------------
!
!      Finalize:
!      ---------
    IF ( Parallel ) THEN 
      CALL ParallelUpdateResult( Matrix1, Solution, Residual )
    END IF
    
    DEALLOCATE( Residual, Solution2, Work2 )
    
    IF ( Level == Solver % MultiGridTotal ) THEN
      WRITE( Message, '(A,F8.2)' ) 'MG iter time: ', CPUTime() - tt
      CALL Info( 'AMGSolve', Message, Level=5 )

      IF(CoarseSave) THEN
        IF(ASSOCIATED(AMG)) THEN
          DO i = MinLevel,Solver % MultiGridTotal
            IF(ASSOCIATED(AMG(i) % InvCF)) DEALLOCATE(AMG(i) % InvCF)
          END DO
          DEALLOCATE(AMG)
        END IF
      END IF
    END IF
 
    RETURN
!------------------------------------------------------------------------------

  CONTAINS

  
!------------------------------------------------------------------------------
    RECURSIVE FUNCTION AMGSweep() RESULT(RNorm)
!------------------------------------------------------------------------------
      INTEGER :: i,j,Rounds
      LOGICAL :: GotIt
      REAL(KIND=dp) :: RNorm
!------------------------------------------------------------------------------
      INTEGER :: Sweeps
      INTEGER, POINTER :: Iters(:)
      REAL(KIND=dp), POINTER :: R1(:),R2(:)

!------------------------------------------------------------------------------

!      Presmoothing:
!      -------------
      
      Iters => ListGetIntegerArray( Params,'MG Sweeps',GotIt)
      IF(GotIt) THEN
        Sweeps = Iters(MIN(InvLevel,SIZE(Iters)))
      ELSE        
        Sweeps = 1
      END IF
      
      PSolver => Solver
      RNorm = MGSmooth( PSolver, Matrix1, Solver % Mesh, Solution, ForceVector, &
             Residual, Level, DOFs, PreSmooth = .TRUE. )
      
!------------------------------------------------------------------------------
!
!      Solve (PAQ)z = Pr, x = x + Qz:
!      ==============================
!
!      Project current residual to the lower level mesh:
!      -------------------------------------------------
      R1 => Residual(1:n)
      R2 => Residual2(1:n2)

      CALL CRS_ProjectVector( ProjPN, R1, R2, DOFs, Trans = .FALSE. )
 
!      Recursively solve (PAQ)z = Pr:
!      ------------------------------
      Solution2 = 0.0d0

!      numbers of W-cycles
      DO i=1,Sweeps
        Work2(1:n2) = Solution2(1:n2)
        
        CALL MultigridSolve( Matrix2, Work2, Residual2, DOFs, &
            Solver, Level-1,NewLinearSystem )

        Solution2(1:n2) = Solution2(1:n2) + Work2(1:n2)
      END DO

!      Compute x = x + Qz:
!      -------------------
      R1 => Residual (1:n)
      R2 => Solution2(1:n2)
      
      CALL CRS_ProjectVector( ProjQT, R2, R1, DOFs, Trans = .TRUE. )      

      Solution(1:n) = Solution(1:n) + Residual(1:n)


!      Post smoothing:
!      ---------------
      RNorm = MGSmooth( PSolver, Matrix1, Solver % Mesh, Solution, ForceVector, &
          Residual, Level, DOFs )
!------------------------------------------------------------------------------
    END FUNCTION AMGSweep
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    FUNCTION MGnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       REAL(KIND=dp) :: s
       REAL(KIND=dp) CONTIG :: x(:)
!------------------------------------------------------------------------------
       IF ( .NOT. Parallel ) THEN
          s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
       ELSE
          s = ParallelNorm( n, x )
       END IF
!------------------------------------------------------------------------------
    END FUNCTION MGnorm
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE MGmv( A, x, b, Update )
!------------------------------------------------------------------------------
       REAL(KIND=dp) CONTIG :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       IF ( .NOT. Parallel ) THEN
         CALL CRS_MatrixVectorMultiply( A, x, b )
       ELSE
         IF ( PRESENT( Update ) ) THEN
           CALL ParallelMatrixVector( A,x,b,Update,ZeroNotOwned=.TRUE. )
         ELSE
           CALL ParallelMatrixVector( A,x,b,ZeroNotOwned=.TRUE. )
         END IF
       END IF
!------------------------------------------------------------------------------
    END SUBROUTINE MGmv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION MGCnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       REAL(KIND=dp) :: s
       COMPLEX(KIND=dp) CONTIG :: x(:)
!------------------------------------------------------------------------------
       s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
!------------------------------------------------------------------------------
    END FUNCTION MGCnorm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE MGCmv( A, x, b, Update )
!------------------------------------------------------------------------------
       COMPLEX(KIND=dp) CONTIG :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       CALL CRS_ComplexMatrixVectorMultiply( A, x, b )
!------------------------------------------------------------------------------
     END SUBROUTINE MGCmv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Apply compatible relaxation for the current C/F split and based on the 
!> convergence behavior include some of the F nodes to the new candidate list.
!------------------------------------------------------------------------------

     SUBROUTINE CompatibleRelaxation(Amat, Solver, CF, CandList, nods, newcands)
       
       TYPE(Matrix_t), POINTER  :: Amat
       TYPE(solver_t), TARGET :: Solver
       INTEGER, POINTER :: CF(:)
       INTEGER, POINTER :: CandList(:)
       INTEGER :: nods, newcands
       
       INTEGER :: i, j, k, Rounds
       REAL(KIND=dp), POINTER :: Ones(:), Zeros(:)
       REAL(KIND=dp) :: Limit, MaxRatio, Ratio
       INTEGER :: cnods, RatioClasses(101)
       LOGICAL :: GotIt
       
       Rounds = ListGetInteger(Params,'MG Compatible Relax Rounds',GotIt)
       IF(.NOT. GotIt) Rounds = 3
       
       ALLOCATE(Ones(nods),Zeros(nods))
       
       k = 0
       DO i=1,nods
         IF(CF(i) <= 0) THEN
           Ones(i) = 1.0
           k = k+1
         ELSE
           Ones(i) = 0.0
         END IF
       END DO
       Zeros = 0.0
       
       IterMethod = ListGetString( Params, 'MG Smoother', Found )
       
       IF ( .NOT. Found ) THEN
         IterMethod = ListGetString( Params,'Linear System Iterative Method', Found )
       END IF
       IF ( .NOT. Found ) IterMethod = 'jacobi'
       
       SELECT CASE( IterMethod )
       CASE( 'gs' )                         
         CALL CR_GS( Amat, Ones, Zeros, CF, Rounds )
         
       CASE( 'sor','sgs','psgs')                                     
         CALL CR_SGS( Amat, Ones, Zeros, CF, Rounds)
         
       CASE( 'csgs')                                     
         CALL CR_CSGS( nods, Amat, Ones, Zeros, CF, Rounds)       
         
       CASE DEFAULT
         CALL CR_Jacobi( Amat, Ones, Zeros, CF, Rounds )
       END SELECT
       
       MaxRatio = 0.0d0
       RatioClasses = 0
       DO i=1,nods
         IF(CF(i) <= 0) THEN
           Ratio = ABS(Zeros(i))
           j = INT(10*Ratio)+1
           IF(j > 0 .AND. j <= 10) THEN
             RatioClasses(j) = RatioClasses(j) + 1
           ELSE IF(j > 0) THEN
             RatioClasses(11) = RatioClasses(11) + 1            
           ELSE IF(j < 0) THEN
             RatioClasses(12) = RatioClasses(12) + 1                        
           END IF
           MaxRatio = MAX(Ratio, MaxRatio) 
         END IF
       END DO
              
       WRITE( Message, '(A)' ) 'Compatible relaxation classes (interval, no and %)'
       CALL Info('CompatibleRelaxation',Message)
       
       DO i=1,13
         IF(RatioClasses(i) > 0) THEN
           IF(i==11) THEN
             WRITE( Message, '(F3.1,A,A,I9,F9.3)' ) 1.0,' - ','...',RatioClasses(i),100.0*RatioClasses(i)/k
           ELSE IF(i==12) THEN
             WRITE( Message, '(A,A,F3.1,I9,F9.3)' ) '...',' - ',0.0,RatioClasses(i),100.0*RatioClasses(i)/k
           ELSE
             WRITE( Message, '(F3.1,A,F3.1,I9,F9.3)' ) 0.1*(i-1),' - ',0.1*i,RatioClasses(i),100.0*RatioClasses(i)/k
           END IF
           CALL Info('CompatibleRelaxation',Message)
         END IF
       END DO
       
       WRITE( Message, '(A,ES15.5)' ) 'Compatible relaxation merit',MaxRatio
       CALL Info('CompatibleRelaxation',Message)
       
       newcands = 0
       Limit = ListGetConstReal(Params,'MG Compatible Relax Limit',GotIt)
       IF(GotIt) THEN
         DO i=1,nods         
           Ratio = ABS(Zeros(i))
           IF(Ratio <= Limit) CYCLE 
           
           IF(CF(i) > 0) THEN
             CALL Warn('CompatibleRelaxation','Coarse nodes should relax well!?')
           END IF
           
           ! Add the node to the list of candidates and eliminate the old C/F info
           newcands = newcands + 1
           CF(i) = 0
           CandList(i) = 1
         END DO
         
         IF(newcands > 0) THEN
           WRITE(Message,'(A,I8)') 'Number of new candidate nodes using CR ',newcands
           CALL Info('CompatibleRelaxation',Message)         
         END IF
       END IF

       DEALLOCATE(Ones, Zeros)
       
     END SUBROUTINE CompatibleRelaxation


!------------------------------------------------------------------------------
!> Create a coarse mesh given the fine mesh and the stiffness matrix. 
!> The coarse nodes may be selected in a number of ways and after the 
!> selection a new mesh is made of them. This mesh is then used to create
!> a projection between the coarse and fine degrees of freedom. 
!------------------------------------------------------------------------------

  SUBROUTINE ChooseCoarseNodes(Amat, Solver, Projector, Components, CF, InvCF) 
    
    TYPE(Matrix_t), POINTER  :: Amat
    TYPE(solver_t), TARGET :: Solver
    TYPE(Matrix_t), POINTER :: Projector
    INTEGER :: Components
    INTEGER, POINTER :: CF(:)
    INTEGER, POINTER, OPTIONAL :: InvCF(:)

    INTEGER :: nods, cnods, newcands, Rounds, RatioClasses(101), elimnods, &
        Component1, cj, CRiter
    INTEGER, POINTER :: CandList(:)
    LOGICAL :: CompMat, PseudoGeometric, UseCR
    LOGICAL, POINTER :: Bonds(:)
    REAL(KIND=dp), POINTER :: Ones(:), Zeros(:)
    INTEGER, POINTER :: Cols(:),Rows(:)

    Component1 = 1
    IF(Components > 1) THEN
      Component1 = ListGetInteger(Params,'MG Determining Component',&
          GotIt,minv=1,maxv=Components)
      IF(.NOT. GotIt) Component1 = 1
    END IF

    CompMat = ListGetLogical(Params,'MG Complex Matrix',GotIt)
    CompMat = (CompMat .OR. .NOT. GotIt) .AND. Amat % COMPLEX
    IF(CompMat) THEN
      CALL Info('ChooseCoarseNodes','Assuming complex valued matrix')
    END IF

    PseudoGeometric = ListGetLogical(Params,'MG Pseudo Geometric',GotIt)

    Rows => Amat % Rows
    Cols => Amat % Cols
    nods = Amat % NumberOfRows
    ALLOCATE( Bonds(SIZE(Amat % Cols)), CandList(nods), CF(nods), Fixed(nods) )

    ! Make the candidate and strong bond list for determining the coarse nodes    
    Fixed = .FALSE.
    CF = 0
    UseCR =  ListGetLogical(Params,'MG Compatible Relax Init',GotIt) 
    IF(UseCR .AND. Components > 1) THEN
      CALL Fatal('CompatibleRelaxation','CR may only be applied to cases with 1 DOFs!')
    END IF
    CRiter = 0

    IF( UseCR ) THEN
      CandList = 0
      Bonds = .TRUE.
      CALL AMGBondsDirichlet(Amat, Bonds, CandList)
      CALL CompatibleRelaxation(Amat, Solver, CF, CandList, nods, newcands)
    ELSE 
      IF( Components > 1) THEN
        CandList = 0
        CandList(Component1:nods:Components) = 1
       ELSE
        CandList = 1
      END IF
      
      IF(CompMat) THEN
        CALL AMGBondsComplex(Amat, Bonds, CandList)
      ELSE IF(PseudoGeometric) THEN
        CALL AMGBondsGeometric(Amat, Bonds, CandList, Components)      
      ELSE
        CALL AMGBonds(Amat, Bonds, CandList, Components)      
      END IF
    END IF


100 CALL AMGCoarse(Amat, CandList, Bonds, CF, CompMat)     

    IF(.NOT. (CompMat .OR. UseCR ) ) THEN
      IF( ListGetLogical(Params,'MG Positive Connection Eliminate',GotIt)) THEN
        CALL AMGPositiveBonds(Amat, Bonds, CandList, CF)
      END IF
    END IF


    IF( Components == 1 .AND. &
        ListGetLogical(Params,'MG Compatible Relax Merit',GotIt) ) THEN
      CRIter = CRiter + 1
      CandList = 0
      CALL CompatibleRelaxation(Amat, Solver, CF, CandList, nods, newcands)
      IF(newcands > 0) GOTO 100
    END IF

    ! Set the CF vector to be zero for fine nodes and an order number for coarse nodes
    cnods = 0
    DO i=1,nods
      IF(CF(i) < 0) THEN
        CF(i) = 0
      ELSE IF(CF(i) > 0) THEN
        cnods = cnods + 1
        CF(i) = cnods
      END IF
    END DO

    WRITE(Message,'(A,F8.3)') 'Coarse node ratio',1.0d0*nods/cnods
    CALL Info('ChooseCoarseNodes',Message)

    DEALLOCATE(Bonds, CandList)

    IF(CompMat) THEN
      Projector => ComplexInterpolateF2C( Amat, CF )      
    ELSE IF(PseudoGeometric) THEN
      Projector => InterpolateF2CDistance( Amat, CF, Components)
    ELSE
      Projector => InterpolateF2C( Amat, CF, Components)
    END IF

    DEALLOCATE(Fixed)

    IF(PRESENT(InvCF)) THEN
      ALLOCATE(InvCF(cnods))
      DO i=1,nods
        IF(CF(i) > 0) InvCF(CF(i)) = i
      END DO
    END IF
    
!    CALL Info('ChooseCoarseNodes','Coarse set chosen')

  END SUBROUTINE ChooseCoarseNodes


!------------------------------------------------------------------------------
!> Create the initial list for measure of importance and
!> make the inverse table for important connections.
!------------------------------------------------------------------------------

  SUBROUTINE AMGBonds(Amat, Bonds, Cands, Components)
    
    LOGICAL, POINTER :: Bonds(:)
    TYPE(Matrix_t), POINTER  :: Amat
    INTEGER, POINTER :: Cands(:)
    INTEGER :: Components

    REAL(KIND=dp) :: NegLim, PosLim
    INTEGER :: nods, cnods, diagsign, maxconn, posnew, negnew, MaxConns, MinConns
    INTEGER :: i,j,k,cj,ci,ind, elimnods,posbonds,negbonds,measind
    INTEGER, POINTER :: Cols(:),Rows(:)
    LOGICAL :: debug, ElimDir, minmaxset, AllowPosLim
    REAL(KIND=dp), POINTER :: Values(:), measures(:)
    REAL(KIND=dp) :: maxbond, minbond, dirlim, meas, measlim

    debug = .FALSE.
    IF(debug) CALL Info('AMGBonds','Making a list of strong matrix connections')

    NegLim = ListGetConstReal(Params,'MG Strong Connection Limit',GotIt)
    IF(.NOT. GotIt) NegLim = 0.06

    ! Negative connections are more useful for the interpolation, but also 
    ! positive strong connection may be taken into account
    AllowPosLim = ListGetLogical(Params,'MG Positive Connection Allow',GotIt)
    PosLim = ListGetConstReal(Params,'MG Positive Connection Limit',GotIt)
    IF(.NOT. GotIt) PosLim = 1.0

    ! In the first time deselect the Dirichlet nodes from the candidate list
    ! their value is determined at the finest level and need not to be recomputed
    ElimDir = EliminateDir
    DirLim = ListGetConstReal(Params,'MG Eliminate Dirichlet Limit',GotIt)
    IF(.NOT. GotIt) DirLim = 1.0d-8      

    nods = Amat % NumberOfRows
    Rows   => Amat % Rows
    Cols   => Amat % Cols
    Values => Amat % Values

    maxconn = 0
    DO ind=1,nods
      maxconn = MAX(maxconn,Rows(ind+1)-Rows(ind))
    END DO
    MaxConns = ListGetInteger(Params,'MG Strong Connection Maximum',GotIt)
    MinConns = ListGetInteger(Params,'MG Strong Connection Minimum',GotIt)


    ALLOCATE(measures(maxconn))

    Bonds = .FALSE.
    posbonds = 0
    negbonds = 0
    elimnods = 0

    DO ind=1,nods

      IF(Cands(ind) == 0) CYCLE

      ! Matrix entries will be treated differently depending if they have the same or
      ! different sign than the diagonal element
      diagsign = 1
      IF(Values (Amat % Diag(ind)) < 0.0) diagsign = -1

      minmaxset = .FALSE.
      DO j=Rows(ind),Rows(ind+1)-1
        cj = Cols(j)

        IF( MOD(ind, Components) /= MOD(cj, Components) ) CYCLE

        IF(Cands(cj) /= 0 .AND. cj /= ind ) THEN
          IF(minmaxset) THEN
            maxbond = MAX(maxbond, Values(j))
            minbond = MIN(minbond, Values(j))
          ELSE
            maxbond = Values(j)
            minbond = Values(j)
            minmaxset = .TRUE.
          END IF
        END IF
      END DO

      IF(.NOT. minmaxset) THEN
        maxbond = 0.0
      ELSE IF(AllowPosLim) THEN
        maxbond = MAX(ABS(maxbond),ABS(minbond))
      ELSE
        IF(minbond * diagsign < 0.0) maxbond = minbond
        maxbond = ABS(maxbond)
      END IF
   
      ! Mark Dirichlet nodes with negative sign in order to favour boundaries in future
      IF( maxbond <= DirLim * ABS(Values (Amat % Diag(ind)) ) ) THEN
        IF(ElimDir) THEN
          Fixed(ind) = .TRUE.
          Cands(ind) = -1
          elimnods = elimnods + 1
        END IF
        CYCLE
      END IF

      IF(.NOT. minmaxset) CYCLE

      ! Make the direct table of important bonds
      posnew = 0
      negnew = 0

      DO j=Rows(ind),Rows(ind+1)-1
        cj = Cols(j)

        ! For multi-component equation assume block structure when building the projector
        IF( MOD(ind, Components) /= MOD(cj, Components) ) CYCLE

        IF(Cands(cj) /= 0 .AND. cj /= ind ) THEN
          IF(diagsign * Values(j) <= 0.0) THEN
            meas = ABS(Values(j)) / (NegLim * maxbond)
            measures(j-Rows(ind)+1) = -meas
            IF( meas > 1.0) THEN
              Bonds(j) = .TRUE.
              negnew = negnew + 1
            END IF
          ELSE IF(AllowPosLim) THEN
            meas = ABS(Values(j)) / (PosLim * maxbond)
            measures(j-Rows(ind)+1) = meas
            IF( meas > 1.0) THEN
              Bonds(j) = .TRUE.
              posnew = posnew + 1
            END IF
          ELSE
            measures(j-Rows(ind)+1) = 0.0d0
          END IF
        END IF
      END DO


      IF(MaxConns > 0) THEN
        DO WHILE(posnew + negnew > MaxConns)
          
          ! Find the weakest used connection
          measlim = HUGE(measlim)
          DO j=Rows(ind),Rows(ind+1)-1
            IF(.NOT. Bonds(j)) CYCLE
            
            meas = ABS(measures(j-Rows(ind)+1))
            IF(meas < measlim) THEN
              measlim = meas
              measind = j
            END IF
          END DO

          IF(measures(measind-Rows(ind)+1) < 0.0) THEN
            negnew = negnew - 1
          ELSE
            posnew = posnew - 1
          END IF
          Bonds(measind) = .FALSE.          
        END DO
      END IF


      IF(MinConns > 0) THEN
        DO WHILE(posnew + negnew < MinConns)
          
          ! Find the strongest unused connection
          measlim = 0.0
          measind = 0
          DO j=Rows(ind),Rows(ind+1)-1
            IF(Bonds(j)) CYCLE
            
            cj = Cols(j)
            IF(Cands(cj) == 0 .OR. cj == ind ) CYCLE
            IF( MOD(ind, Components) /= MOD(cj, Components) ) CYCLE
 
            meas = ABS(measures(j-Rows(ind)+1))
            IF(meas > measlim) THEN
              measlim = meas
              measind = j
            END IF
          END DO

          ! Check if there exist possible new connections
          IF(measind == 0 .OR. measlim > 1.0d-50) EXIT

          IF(measures(measind-Rows(ind)+1) < 0.0) THEN
            negnew = negnew + 1
          ELSE
            posnew = posnew + 1
          END IF
          Bonds(measind) = .TRUE.          
        END DO
      END IF

      posbonds = posbonds + posnew
      negbonds = negbonds + negnew
    END DO


    IF(elimnods > 0) THEN
      WRITE(Message,'(A,I8)') 'Number of eliminated nodes',elimnods
      CALL Info('AMGBonds',Message)
    END IF
!    WRITE(Message,'(A,I8)') 'Number of possible connections',SIZE(Bonds)
!    CALL Info('AMGBonds',Message)
!    j = posbonds + negbonds
!    WRITE(Message,'(A,I8)') 'Number of strong connections',j
!    CALL Info('AMGBonds',Message)
    IF(posbonds > 0) THEN
      WRITE(Message,'(A,I8)') 'Number of positive connections',posbonds
      CALL Info('AMGBonds',Message)
    END IF
    WRITE(Message,'(A,F8.3)') 'Average number of strong bonds',1.0*Components*j/nods
    CALL Info('AMGBonds',Message)

  END SUBROUTINE AMGBonds


!------------------------------------------------------------------------------
!> Mark the Dirichlet bonds in the AMG target matrix.
!------------------------------------------------------------------------------
  SUBROUTINE AMGBondsDirichlet(Amat, Bonds, Cands)
    
    LOGICAL, POINTER :: Bonds(:)
    TYPE(Matrix_t), POINTER  :: Amat
    INTEGER, POINTER :: Cands(:)
    INTEGER :: Components, Component1    

    REAL(KIND=dp) :: NegLim, PosLim
    INTEGER :: nods, cnods, diagsign, maxconn, posnew, negnew, MaxConns, MinConns
    INTEGER :: i,j,k,cj,ci,ind, elimnods,posbonds,negbonds,measind
    INTEGER, POINTER :: Cols(:),Rows(:)
    LOGICAL :: debug, ElimDir, minmaxset, AllowPosLim
    REAL(KIND=dp), POINTER :: Values(:), measures(:)
    REAL(KIND=dp) :: maxbond, minbond, dirlim, meas, measlim

    debug = .FALSE.
    IF(debug) CALL Info('AMGBondsDirichlet','Marking the dirichlet conditions')

    ! In the first time deselect the Dirichlet nodes from the candidate list
    ! their value is determined at the finest level and need not to be recomputed
    ElimDir = EliminateDir
    IF(.NOT. ElimDir ) RETURN

    DirLim = ListGetConstReal(Params,'MG Eliminate Dirichlet Limit',GotIt)
    IF(.NOT. GotIt) DirLim = 1.0d-8      

    nods = Amat % NumberOfRows
    Rows   => Amat % Rows
    Cols   => Amat % Cols
    Values => Amat % Values

    elimnods = 0
    DO ind=1,nods

      IF(Cands(ind) == 0) CYCLE

      maxbond = 0.0d0
      DO j=Rows(ind),Rows(ind+1)-1
        cj = Cols(j)
        IF(Cands(cj) /= 0 .AND. cj /= ind ) THEN
          maxbond = ABS( Values(j) )
        END IF
      END DO

      ! Mark Dirichlet nodes with negative sign in order to favour boundaries in future
      IF( maxbond < DirLim * ABS(Values (Amat % Diag(ind)) ) ) THEN
        IF(ElimDir) THEN
          Fixed(ind) = .TRUE.
          Cands(ind) = -1
          elimnods = elimnods + 1
        END IF
      END IF
    END DO
    
    IF(elimnods > 0) THEN
      WRITE(Message,'(A,I8)') 'Number of eliminated nodes',elimnods
      CALL Info('AMGBondsDirichlet',Message)
    END IF

  END SUBROUTINE AMGBondsDirichlet

!------------------------------------------------------------------------------
!> Create the initial list for measure of importance using geometric distance 
!> information and topology information from the matrix..
!------------------------------------------------------------------------------

  SUBROUTINE AMGBondsGeometric(Amat, Bonds, Cands, Components)
    
    LOGICAL, POINTER :: Bonds(:)
    TYPE(Matrix_t), POINTER  :: Amat
    INTEGER, POINTER :: Cands(:)
    INTEGER :: Components

    INTEGER :: nods, cnods, diagsign, maxconn, posnew, negnew, MaxConns, MinConns
    INTEGER :: i,j,k,cj,ci,ind, elimnods,posbonds,nobonds,measind
    INTEGER, POINTER :: Cols(:),Rows(:)
    LOGICAL :: debug, ElimDir, minmaxset
    REAL(KIND=dp), POINTER :: Values(:), measures(:)
    REAL(KIND=dp) :: maxbond, minbond, dirlim, meas, measlim, Neglim, &
        x0, y0, z0, x1, y1, z1, s2, Pow

    debug = .FALSE.
    IF(debug) CALL Info('AMGBonds','Making a list of strong matrix connections')

    NegLim = ListGetConstReal(Params,'MG Strong Connection Limit',GotIt)
    IF(.NOT. GotIt) NegLim = 0.25
    
    Pow = ListGetConstReal(Params,'MG Geometric Power',GotIt)
    IF(.NOT. GotIt) Pow = 1.0d0

    ElimDir = EliminateDir
    DirLim = ListGetConstReal(Params,'MG Eliminate Dirichlet Limit',GotIt)
    IF(.NOT. GotIt) DirLim = 1.0d-8      

    nods   = Amat % NumberOfRows
    Rows   => Amat % Rows
    Cols   => Amat % Cols
    Values => Amat % Values

    maxconn = 0
    DO ind=1,nods
      maxconn = MAX(maxconn,Rows(ind+1)-Rows(ind))
    END DO
    MaxConns = ListGetInteger(Params,'MG Strong Connection Maximum',GotIt)
    MinConns = ListGetInteger(Params,'MG Strong Connection Minimum',GotIt)

    ALLOCATE(measures(maxconn))

    Bonds = .FALSE.
    nobonds = 0
    elimnods = 0

    DO ind=1,nods

      IF(Cands(ind) == 0) CYCLE

      minmaxset = .FALSE.
      maxbond = 0.0
       
      IF ( Level == Solver % MultiGridTotal ) THEN
        k = ind
      ELSE
        k = AMG(Level+1) % InvCF(ind)
      END IF
      k = (k-1) / Components + 1

      x0 = Mesh % Nodes % x(k)
      y0 = Mesh % Nodes % y(k)
      z0 = Mesh % Nodes % z(k)

      DO j=Rows(ind),Rows(ind+1)-1
        cj = Cols(j)
        measures(j-Rows(ind)+1) = 0.0d0   

        IF(Cands(cj) == 0 .OR. cj == ind ) CYCLE
        IF(ElimDir .AND. ABS(Values(j)) < Dirlim) CYCLE 

        minmaxset = .TRUE.          
        
        IF ( Level == Solver % MultiGridTotal ) THEN
          k = cj
        ELSE
          k = AMG(Level+1) % InvCF(cj)
        END IF
        k = (k-1) / Components + 1
        
        x1 = Mesh % Nodes % x(k)
        y1 = Mesh % Nodes % y(k)
        z1 = Mesh % Nodes % z(k)
        
        s2 = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0)
        meas = s2 ** (-Pow/2)
        
        measures(j-Rows(ind)+1) = meas
        maxbond = MAX(maxbond, meas)
      END DO

      IF(ElimDir .AND. maxbond < Dirlim) THEN
        Cands(ind) = -1
        elimnods = elimnods + 1
        CYCLE
      END IF

      IF(.NOT. minmaxset) CYCLE

      ! Make the direct table of important bonds
      negnew = 0

      DO j=Rows(ind),Rows(ind+1)-1
        cj = Cols(j)
        IF(Cands(cj) /= 0 .AND. cj /= ind ) THEN
          meas = measures(j-Rows(ind)+1) / (Neglim * maxbond)
          measures(j-Rows(ind)+1) = meas
          IF( meas > 1.0) THEN
            Bonds(j) = .TRUE.
            negnew = negnew + 1
          END IF
        END IF
      END DO


      IF(MaxConns > 0) THEN
        DO WHILE(negnew > MaxConns)
          
          ! Find the weakest used connection and drop it
          measlim = HUGE(measlim)
          DO j=Rows(ind),Rows(ind+1)-1
            IF(.NOT. Bonds(j)) CYCLE
            
            meas = measures(j-Rows(ind)+1)
            IF(meas < measlim) THEN
              measlim = meas
              measind = j
            END IF
          END DO

          negnew = negnew - 1
          Bonds(measind) = .FALSE.          
        END DO
      END IF


      IF(MinConns > 0) THEN
        DO WHILE(negnew < MinConns)
          
          ! Find the strongest unused connection
          measlim = 0.0
          measind = 0
          DO j=Rows(ind),Rows(ind+1)-1
            IF(Bonds(j)) CYCLE
            
            cj = Cols(j)
            IF(Cands(cj) == 0 .OR. cj == ind ) CYCLE
 
            meas = measures(j-Rows(ind)+1)
            IF(meas > measlim) THEN
              measlim = meas
              measind = j
            END IF
          END DO

          ! Check if there exist possible new connections
          IF(measind == 0 .OR. measlim < 1.0d-50) EXIT

          negnew = negnew + 1
          Bonds(measind) = .TRUE.          
        END DO
      END IF

      nobonds = nobonds + negnew

    END DO

    IF(elimnods > 0) THEN
      WRITE(Message,'(A,I8)') 'Number of eliminated nodes',elimnods
      CALL Info('AMGBondsGeometric',Message)
    END IF
    WRITE(Message,'(A,I8)') 'Number of possible connections',SIZE(Bonds)
    CALL Info('AMGBondsGeometric',Message)
    WRITE(Message,'(A,I8)') 'Number of strong connections',nobonds
    CALL Info('AMGBondsGeometric',Message)
    WRITE(Message,'(A,F8.3)') 'Strong bonds for each dof',1.0*nobonds/nods
    CALL Info('AMGBondsGeometric',Message)

  END SUBROUTINE AMGBondsGeometric


!------------------------------------------------------------------------------
!> Create the initial list for measure of importance for a complex matrix.
!------------------------------------------------------------------------------

  SUBROUTINE AMGBondsComplex(Amat, Bonds, Cands)
    
    LOGICAL, POINTER :: Bonds(:)
    TYPE(Matrix_t), POINTER  :: Amat
    INTEGER, POINTER :: Cands(:)

    REAL(KIND=dp) :: NegLim
    INTEGER :: nods, cnods, anods, maxconn, negnew, MaxConns, MinConns
    INTEGER :: i,j,k,cj,ci,ind, elimnods, negbonds, measind,j2
    INTEGER, POINTER :: Cols(:),Rows(:)
    LOGICAL :: ElimDir
    REAL(KIND=dp), POINTER :: Values(:), measures(:)
    REAL(KIND=dp) :: maxbond, dirlim, meas, measlim

    CALL Info('AMGBondsComplex','Making a list of strong matrix connections')

    NegLim = ListGetConstReal(Params,'MG Strong Connection Limit',GotIt)
    IF(.NOT. GotIt) NegLim = 0.25

    ! In the first time deselect the Dirichlet nodes from the candidate list
    ! their value is determined at the finest level and need not to be recomputed
    ElimDir = EliminateDir
    DirLim = ListGetConstReal(Params,'MG Eliminate Dirichlet Limit',GotIt)
    IF(.NOT. GotIt) DirLim = 1.0d-8      

    nods = Amat % NumberOfRows
    Rows   => Amat % Rows
    Cols   => Amat % Cols
    Values => Amat % Values

    maxconn = 0
    DO ind=1,nods
      maxconn = MAX(maxconn,Rows(ind+1)-Rows(ind))
    END DO
    MaxConns = ListGetInteger(Params,'MG Strong Connection Maximum',GotIt)
    MinConns = ListGetInteger(Params,'MG Strong Connection Minimum',GotIt)


    ALLOCATE(measures(maxconn))

    Bonds = .FALSE.
    negbonds = 0

    DO ind=1,nods

      IF(Cands(ind) == 0) CYCLE

      maxbond = 0.0d0
      DO j=Rows(ind),Rows(ind+1)-1
        cj = Cols(j)
        IF(Cands(cj) == 0 .OR. cj == ind ) CYCLE
        
        j2 = j + Rows(ind+1) - Rows(ind)
        meas = SQRT( Values(j)**2 + Values(j2)**2 )
        maxbond = MAX(meas,maxbond)
      END DO
     
      meas = SQRT(Values (Amat % Diag(ind)) ** 2 + Values (Amat % Diag(ind+1)) ** 2 )

      IF( maxbond < DirLim * meas ) THEN
        IF(ElimDir) THEN
          Fixed(ind) = .TRUE.
          Cands(ind) = 0
          elimnods = elimnods + 1
        END IF
        CYCLE
      END IF

      ! Make the direct table of important bonds
      negnew = 0

      DO j=Rows(ind),Rows(ind+1)-1
        cj = Cols(j)
        IF(Cands(cj) == 0 .OR. cj == ind ) CYCLE

        j2 = j + Rows(ind+1) - Rows(ind)
        meas = SQRT( Values(j)**2 + Values(j2)**2 ) / (NegLim * maxbond)
        measures(j-Rows(ind)+1) = meas
        IF(meas > 1.0) THEN
          Bonds(j) = .TRUE.
          negnew = negnew + 1
        END IF
      END DO


      IF(MaxConns > 0) THEN
        DO WHILE(negnew > MaxConns)
          
          ! Find the weakest used connection
          measlim = HUGE(measlim)
          DO j=Rows(ind),Rows(ind+1)-1
            IF(.NOT. Bonds(j)) CYCLE
            
            meas = measures(j-Rows(ind)+1)
            IF(meas < measlim) THEN
              measlim = meas
              measind = j
            END IF
          END DO

          negnew = negnew - 1
          Bonds(measind) = .FALSE.          
        END DO
      END IF


      IF(MinConns > 0) THEN
        DO WHILE(negnew < MinConns)
          
          ! Find the strongest unused connection
          measlim = 0.0
          measind = 0
          DO j=Rows(ind),Rows(ind+1)-1
            IF(Bonds(j)) CYCLE
            
            cj = Cols(j)
            IF(Cands(cj) == 0 .OR. cj == ind ) CYCLE
 
            meas = measures(j-Rows(ind)+1)
            IF(meas > measlim) THEN
              measlim = meas
              measind = j
            END IF
          END DO

          ! Check if there exist possible new connections
          IF(measind == 0 .OR. measlim > 1.0d-50) EXIT

          negnew = negnew + 1
          Bonds(measind) = .TRUE.          
        END DO
      END IF
      negbonds = negbonds + negnew
    END DO

    WRITE(Message,'(A,I8)') 'Number of eliminated nodes',elimnods
    CALL Info('AMGBondsComplex',Message)
    j = negbonds
    WRITE(Message,'(A,I8)') 'Number of strong connections',j
    CALL Info('AMGBondsComplex',Message)
    WRITE(Message,'(A,F8.3)') 'Average number of strong bonds for each dof',2.0*j/nods
    CALL Info('AMGBondsComplex',Message)
    WRITE(Message,'(A,F8.3)') 'Fraction of strong bonds in sparse matrix',4.0*j/SIZE(Bonds)
    CALL Info('AMGBondsComplex',Message)

  END SUBROUTINE AMGBondsComplex



!------------------------------------------------------------------------------
!> Add a posteriori some nodes to be coarse nodes which have 
!> strong positive-positive connections
!------------------------------------------------------------------------------
  SUBROUTINE AMGPositiveBonds(Amat, Bonds, Cands, CF)
    
    LOGICAL, POINTER :: Bonds(:)    
    TYPE(Matrix_t), POINTER  :: Amat
    INTEGER, POINTER :: Cands(:), CF(:)

    INTEGER :: i, j, k, cj, ci, ind, nods, posnods
    INTEGER, POINTER :: Cols(:),Rows(:)
    LOGICAL :: ElimDir, AllowPosLim
    REAL(KIND=dp) :: maxbond, minbond, minbond2, PosLim, DirLim, diagvalue
    REAL(KIND=dp), POINTER :: Values(:)

    CALL Info('AMGPositiveBonds','Adding some F-nodes with positive connections to C-nodes')

    ! Negative connections are more useful for the interpolation, but also 
    ! positive strong connection may be taken into account
    AllowPosLim = ListGetLogical(Params,'MG Positive Connection Allow',GotIt)
    PosLim = ListGetConstReal(Params,'MG Positive Connection Limit',GotIt)
    IF(.NOT. GotIt) PosLim = 0.5

    ! In the first time deselect the Dirichlet nodes from the candidate list
    ! their value is determined at the finest level and need not to be recomputed
    ElimDir = EliminateDir
    DirLim = ListGetConstReal(Params,'MG Eliminate Dirichlet Limit',GotIt)
    IF(.NOT. GotIt) DirLim = 1.0d-8      

    nods = Amat % NumberOfRows
    Rows   => Amat % Rows
    Cols   => Amat % Cols
    Values => Amat % Values

    posnods = 0

    DO ind=1,nods
      
      IF(Cands(ind) == 0 .OR. CF(ind) > 0) CYCLE

      ! Matrix entries will be treated differently depending if they have the same or
      ! different sign than the diagonal element

      diagvalue = Values (Amat % Diag(ind))
      minbond = 0.0
      minbond2 = 0.0
      maxbond = 0.0      

      DO j=Rows(ind),Rows(ind+1)-1
        cj = Cols(j)
        IF(Cands(cj) == 0) CYCLE
        IF(cj == ind) CYCLE

        maxbond = MAX(maxbond,ABS(Values(j)))
        IF(diagvalue * Values(j) > 0.0) THEN
          minbond = minbond + ABS(Values(j))
          IF(CF(cj) <= 0) minbond2 = MAX(minbond2,ABS(Values(j)))
        END IF
      END DO

      IF(maxbond < DirLim * ABS(diagvalue)) CYCLE
      IF(minbond2 < minbond) CYCLE 

      IF(minbond2 > PosLim * maxbond) THEN
        CF(ind) = 1
        posnods = posnods + 1
      END IF
        
    END DO

    WRITE(Message,'(A,I9)') 'Number of added positive connection nodes',posnods
    CALL Info('AMGPositiveBonds',Message)

  END SUBROUTINE AMGPositiveBonds


!------------------------------------------------------------------------------
!> Creates a coarse mesh using the given list of strong connections. 
!> Only nodes assigned by the Cands vector may be included in the 
!> coarse set - others are ignored. The subroutine returns the vector CF which is 
!> nonzero for coarse nodes. The nodes are chosen using a heuristics which
!> takes into account the strength of the coupling.
!
! There are three possible coarsening modes
! 0) find a new start only in dire need
! 1) find a new start if the suggested start is worse than previous new start
! 2) find a new start if the suggested start is worse that theoretical best start
! 3) keep list of all points and always choose the highest priority
!
! CF includes the classification of the nodes
! coarse nodes > 0, undecided = 0, fine nodes < 0
!------------------------------------------------------------------------------

  SUBROUTINE AMGCoarse(Amat, Cands, Bonds, CF, CompMat)
    
    TYPE(Matrix_t), POINTER :: Amat
    LOGICAL :: Bonds(:)
    INTEGER, POINTER :: Cands(:), CF(:)
    LOGICAL :: CompMat

    INTEGER :: nods, cnods    
    INTEGER :: i,j,k,cj,ci,ind,maxi,maxi2,newstarts, loops, minneigh, maxneigh, &
        MaxCon, MaxConInd, RefCon, RefCon2, prevstart, cind, find, ind0, &
        CoarseningMode, GatMax, anods, FollowBounds, points, &
        neworder, oldorder
    INTEGER, POINTER :: Con(:)
    LOGICAL :: debug
    INTEGER, POINTER :: Cols(:),Rows(:)
    INTEGER, ALLOCATABLE :: GatLims(:), GatNos(:), ConInd(:), RevConInd(:)
    REAL(KIND=dp), POINTER :: Values(:)
    
    nods = Amat % NumberOfRows
    Rows   => Amat % Rows
    Cols   => Amat % Cols
    Values => Amat % Values

    FollowBounds = ListGetInteger(Params,'MG Boundary Priority',GotIt)
    IF(.NOT. GotIt) FollowBounds = 1
    CoarseningMode = ListGetInteger(Params,'MG Coarsening Mode',&
        GotIt,minv=0,maxv=3)
    IF(.NOT. GotIt) CoarseningMode = 3

    debug = .FALSE.

    ALLOCATE(Con(nods))
    Con = 0
    newstarts = 0
    cnods = 0
    loops = 0
    cind = 0
    find = 0

    ! Make the tightly bonded neighbours of coarse nodes to be fine nodes 
    ! This has an effect only if the CF list has been initialized
    DO ind = 1, nods
      IF(CF(ind) == 0) CYCLE
      
      IF(CF(ind) < 0) THEN
        find = find + 1
        CF(ind) = -find
        CYCLE
      END IF

      cind = cind + 1
      CF(ind) = cind

      IF(Cands(ind) == 0) CYCLE

      DO i=Rows(ind),Rows(ind+1)-1
        IF(.NOT. Bonds(i)) CYCLE
        ci = Cols(i)

        IF(Cands(ci) == 0) CYCLE
        IF(CF(ci) == 0) THEN
          find = find + 1
          CF(ci) = -find
        END IF
      END DO
    END DO
    IF(find+cind > 0) PRINT *,'A priori set fine and coarse nodes',find,cind,nods-find-cind


    ! Calculate the initial measure of of importance
    ! unvisited neighbour get 1 point, and a decided coarse node 2 points   

    DO ind = 1, nods      
      ! Points are obtained when a neighbour is undicided or fine node
      IF(CF(ind) > 0) CYCLE
      IF(Cands(ind) <= 0) CYCLE
      
      DO i=Rows(ind),Rows(ind+1)-1
        IF(.NOT. Bonds(i)) CYCLE
        ci = Cols(i)

        IF(Cands(ci) > 0) THEN
          ! Only give points to undicided points
          IF(CF(ci) == 0) THEN
            IF(CF(ind) == 0) Con(ci) = Con(ci) + 1
            IF(CF(ind) < 0)  Con(ci) = Con(ci) + 2
          END IF
        END IF
      END DO
    END DO


    ! Add some weight if the node is connected to a Dirichlet node
    ! This way the coarse nodes selection will follow the natural boundaries
    IF(FollowBounds /= 0) THEN
      DO ind = 1, nods
        IF(Cands(ind) <= 0) CYCLE
        DO i=Rows(ind),Rows(ind+1)-1
          IF(.NOT. Bonds(i)) CYCLE
          ci = Cols(i)

          IF(Cands(ci) == -1) THEN
            Con(ind) = Con(ind) + FollowBounds
          END IF
        END DO
      END DO
      IF(FollowBounds < 0) THEN
        Con = Con - MINVAL(Con)
      END IF
    END IF

    ! Boundaries may have been marked with negative index, remove that
    DO ind = 1, nods
      Cands(ind) = MAX(Cands(ind),0)
    END DO


    IF(CoarseningMode == 3) THEN
      MaxCon = MAXVAL( Con ) 
      GatMax = 2 * MaxCon + 1
      
      ! Bookkeeping is kept on category boundaries of the points
      ALLOCATE( GatLims(GatMax), GatNos(GatMax) )
      GatLims = 0
      GatNos = 0
      
      ! Number of points in different categories 
      DO ind = 1, nods
        IF(Con(ind) == 0) CYCLE
        GatNos(Con(ind)) = GatNos(Con(ind)) + 1 
      END DO
      anods = SUM(GatNos)

      DO ind = 2, MaxCon
        GatLims(ind) = GatLims(ind-1) + GatNos(ind-1)
      END DO
    
      ALLOCATE( ConInd(nods), RevConInd(anods) )
      ConInd = 0
      RevConInd = 0
      
      GatNos = 0
      DO ind = 1, nods
        IF(Con(ind) == 0) CYCLE
        GatNos(Con(ind)) = GatNos(Con(ind)) + 1 
        ConInd(ind) = GatNos(Con(ind)) + GatLims(Con(ind))
      END DO

      RevConInd = 0
      DO ind = 1, nods
        IF(Con(ind) == 0) CYCLE
        RevConInd(ConInd(ind)) = ind
      END DO
      

      DO WHILE( MaxCon > 0 ) 
        ind = RevConInd(anods)
 
        cind = cind + 1
        CF(ind) = cind
        cnods = cnods + 1
        
        ! Go through all strongly bonded neighbours to 'ind'
        DO j=Rows(ind),Rows(ind+1)-1
          IF(.NOT. Bonds(j)) CYCLE
          cj = Cols(j)
          
          IF(CF(cj) == 0) THEN
            find = find + 1
            CF(cj) = -find
            IF(Cands(cj) == 0) CYCLE
            
            ! Recompute the measure of importance for the neighbours
            DO i=Rows(cj),Rows(cj+1)-1
              IF(Bonds(i)) THEN
                ci = Cols(i)          
                IF(Cands(ci) > 0 .AND. CF(ci) == 0) THEN
                  loops = loops + 1

                  points = Con(ci) + 1
                  
                  oldorder = ConInd(ci)
                  IF(GatLims(points) == 0) GatLims(points) = anods
                  neworder = GatLims(points)
                  Con(ci) = points
                  GatLims(points) = GatLims(points) - 1
                  
                  IF(neworder /= oldorder) THEN
                    k = RevConInd(neworder)
                    
                    ConInd(ci) = neworder
                    ConInd(k) = oldorder
                    
                    RevConInd(neworder) = ci
                    RevConInd(oldorder) = k
                  END IF

                END IF
              END IF
            END DO
          END IF
        END DO

        ! Shorten the list from top in case the node is already set
        Refcon = 0
        ! DO WHILE( anods > 0 .AND. CF(RevConInd(anods)) /= 0 )
        DO WHILE( anods > 0)
          IF (CF(RevConInd(anods)) == 0) EXIT
          MaxCon = Con(RevConInd(anods))          
          Refcon = MAX(MaxCon,RefCon)
          anods = anods - 1
          Refcon = MAX(MaxCon,RefCon)            
        END DO

        IF(anods > 0) THEN
          MaxCon = Con(RevConInd(anods))
          RefCon = MAX(MaxCon,RefCon)
          GatLims((MaxCon+1):(Refcon+1)) = 0
        ELSE        
          MaxCon = 0
        END IF
      END DO

      DEALLOCATE(GatLims, GatNos, ConInd, RevConInd)

    ELSE

      ! Find the point to start the coarse node selection from
      MaxCon = 1000 
      prevstart = 1
      
10    MaxConInd =  0
      RefCon = 0
      newstarts = newstarts + 1
      DO ind=prevstart, nods
        loops = loops + 1
        IF(Cands(ind) == 0) CYCLE
        IF(Con(ind) > RefCon) THEN
          RefCon = Con(ind)
          MaxConInd = ind
        END IF
        IF(RefCon >= MaxCon) EXIT
      END DO
      
      IF(RefCon < MaxCon .AND. prevstart > 1) THEN
        DO ind=1, prevstart-1        
          loops = loops + 1
          IF(Cands(ind) == 0) CYCLE
          IF(Con(ind) > RefCon) THEN
            RefCon = Con(ind)
            MaxConInd = ind
          END IF
          IF(RefCon >= MaxCon) EXIT
        END DO
      END IF
      
      MaxCon = RefCon
      ind = MaxConInd
      prevstart = MAX(1,ind)
      ind0 = prevstart
      
      
      DO WHILE( MaxCon > 0) 
        
        cind = cind + 1
        CF(ind) = cind
        Con(ind) = 0
        cnods = cnods + 1
        
        ! Go through all strongly bonded neighbours to 'ind'
        DO j=Rows(ind),Rows(ind+1)-1
          IF(.NOT. Bonds(j)) CYCLE
          cj = Cols(j)
          
          IF(CF(cj) == 0) THEN
            find = find + 1
            CF(cj) = -find
            Con(cj) = 0
            IF(Cands(cj) == 0) CYCLE
            
            ! Recompute the measure of importance for the neighbours
            DO i=Rows(cj),Rows(cj+1)-1
              IF(Bonds(i)) THEN
                ci = Cols(i)          
                IF(Cands(ci) > 0 .AND. CF(ci) == 0) THEN
                  Con(ci) = Con(ci) + 1
                END IF
              END IF
            END DO
          END IF
        END DO
        
        ! The next candidate is probably among the secondary neighbours
        RefCon = 0
        RefCon2 = 0
        maxi = 0
        maxi2 = 0
        
        
        IF(CoarseningMode == 2) THEN       
          DO j=Rows(ind),Rows(ind+1)-1
            IF(.NOT. Bonds(j)) CYCLE
            cj = Cols(j)
            IF(Cands(cj) == 0) CYCLE
            
            DO i=Rows(cj),Rows(cj+1)-1
              IF(.NOT. Bonds(i)) CYCLE
              ci = Cols(i)
              
              IF(Cands(ci) == 0 .OR. CF(ci) /= 0) CYCLE
              
              IF(Con(ci) > RefCon .AND. ci /= maxi) THEN
                RefCon2 = RefCon
                maxi2 = maxi
                RefCon = Con(ci)
                maxi = ci
              ELSE IF (Con(ci) > RefCon2 .AND. ci /= maxi) THEN
                RefCon2 = Con(ci)
                maxi2 = ci
              END IF
            END DO
          END DO
          
          ! If none of the neighbouring nodes is a potential node go and find a new candidate
          IF(RefCon < MaxCon) THEN
            GOTO 10
          END IF
          
          MaxCon = MAX(RefCon2,MaxCon)
          
        ELSE
          maxneigh = 0
          DO j=Rows(ind),Rows(ind+1)-1
            IF(.NOT. Bonds(j)) CYCLE
            cj = Cols(j)
            IF(Cands(cj) == 0) CYCLE
            
            DO i=Rows(cj),Rows(cj+1)-1
              IF(.NOT. Bonds(i)) CYCLE
              ci = Cols(i)
              
              IF(Cands(ci) == 0 .OR. CF(ci) /= 0) CYCLE
              
              IF((Con(ci) > RefCon) .OR. (Con(ci) == RefCon .AND. -CF(cj) > maxneigh)) THEN
                RefCon = Con(ci)
                maxi = ci
                maxneigh = MAX(maxneigh,-CF(cj))
              END IF
            END DO
          END DO
          
          IF(RefCon < MaxCon) THEN
            
            DO j=Rows(ind0),Rows(ind0+1)-1
              IF(.NOT. Bonds(j)) CYCLE
              cj = Cols(j)
              IF(Cands(cj) == 0) CYCLE
              
              DO i=Rows(cj),Rows(cj+1)-1
                IF(.NOT. Bonds(i)) CYCLE
                ci = Cols(i)
                
                IF(Cands(ci) == 0 .OR. CF(ci) /= 0) CYCLE
                
                IF(Con(ci) > RefCon2) THEN
                  RefCon2 = Con(ci)
                  maxi2 = ci
                  IF(CF(cj) < 0) minneigh = MIN(minneigh,-CF(cj))
                ELSE IF (Con(ci) == RefCon2) THEN
                  IF(CF(cj) < 0 .AND. -CF(cj) < minneigh) THEN
                    RefCon2 = Con(ci)
                    maxi2 = ci
                    minneigh = -CF(cj)
                  END IF
                END IF
              END DO
            END DO
            
            ! Favor the vicinity of the previous starting point
            IF(RefCon2 >= RefCon) THEN
              maxi = maxi2
              RefCon = RefCon2
            END IF
            ind0 = maxi
          END IF
          
          IF(CoarseningMode == 1 .AND. RefCon < MaxCon) GOTO 10
          
          IF(CoarseningMode == 0 .AND. RefCon < 1) GOTO 10
          
        END IF

        ind = maxi
      END DO

      WRITE(Message,'(A,I8)') 'Coarsening algorithm starts',newstarts
      CALL Info('AMGCoarse',Message)

    END IF


    ! Assign unassigned nodes in the candidate list to be coarse nodes
    k = 0
    DO ind = 1, nods
      IF(Cands(ind) > 0 .AND. CF(ind) == 0) THEN
        k = k + 1
        cind = cind + 1
        CF(ind) = cind
      END IF
    END DO

    IF(k > 0) THEN
      WRITE(Message,'(A,I8)') 'Enforced coarse nodes',k
      CALL Info('AMGCoarse',Message)      
    END IF

!    WRITE(Message,'(A,I8)') 'Coarsening algorithm tests',loops
!    CALL Info('AMGCoarse',Message)
    WRITE(Message,'(A,I8)') 'Number of coarse nodes',cind
    CALL Info('AMGCoarse',Message)
!    WRITE(Message,'(A,I8)') 'Number of fine nodes',find
!    CALL Info('AMGCoarse',Message)

  END SUBROUTINE AMGCoarse



!-----------------------------------------------------------------------------
!>     make a matrix projection such that fine nodes are expressed with coarse nodes
!>     Xf = P Xc, with given set of coarse nodes, CF, and the stiffness matrix
!>     used in the projection.
!------------------------------------------------------------------------------
     FUNCTION InterpolateF2C( Fmat, CF, DOFs) RESULT (Projector)
!------------------------------------------------------------------------------
       USE Interpolation
       USE CRSMatrix
       USE CoordinateSystems
!-------------------------------------------------------------------------------
       TYPE(Matrix_t), TARGET  :: Fmat
       INTEGER, POINTER :: CF(:)
       INTEGER :: DOFs
       TYPE(Matrix_t), POINTER :: Projector
!------------------------------------------------------------------------------
       INTEGER, PARAMETER :: FSIZE=1000, CSIZE=100
       INTEGER :: i, j, k, l, Fdofs, Cdofs, ind, ci, cj, Component1, Components, node
       REAL(KIND=dp), POINTER :: PValues(:), FValues(:)
       INTEGER, POINTER :: FRows(:), FCols(:), PRows(:), PCols(:), CoeffsInds(:)
       REAL(KIND=dp) :: bond, VALUE, possum, negsum, poscsum, negcsum, diagsum, &
           ProjLim, negbond, posbond, maxbond
       LOGICAL :: Debug, AllocationsDone, DirectInterpolate, Lumping
       INTEGER :: inds(FSIZE), posinds(CSIZE), neginds(CSIZE), no, diag, InfoNode, &
           posno, negno, negi, posi, projnodes, DirectLimit, entries
       REAL(KIND=dp) :: coeffs(FSIZE), poscoeffs(CSIZE), negcoeffs(CSIZE), wsum, &
           refbond, posmax, negmax, favorneg

       Debug = .FALSE.


       IF(Debug) CALL Info('InterpolateF2C','Starting interpolation')
       
       ProjLim = ListGetConstReal(Params,'MG Projection Limit',GotIt)
       IF(.NOT. GotIt) ProjLim = 0.6

       Lumping = ListGetLogical(Params,'MG Projection Lumping',GotIt)

       DirectInterpolate = ListGetLogical(Params,'MG Direct Interpolate',GotIt)
       DirectLimit = ListGetInteger(Params,'MG Direct Interpolate Limit',GotIt)      
       IF(.NOT. GotIt) THEN
         IF(DirectInterpolate) THEN
           DirectLimit = 0
         ELSE
           DirectLimit = HUGE(DirectLimit)
         END IF
       END IF

       FavorNeg = ListGetConstReal(Params,'MG Negative Connection Favor Factor',GotIt)
       IF(.NOT. GotIt) FavorNeg = 1.0d0

       Components = DOFs
       Component1 = 1
       IF(Components > 1) THEN
         Component1 = ListGetInteger(Params,'MG Determining Component',GotIt,minv=1,maxv=DOFs)
         IF(.NOT. GotIt) Component1 = 1
       END IF

       InfoNode = ListGetInteger(Params,'MG Info Node',GotIt)

       Fdofs = Fmat % NumberOfRows
       FRows   => Fmat % Rows
       FCols   => Fmat % Cols
       FValues => Fmat % Values
       
       ALLOCATE(CoeffsInds(Fdofs))
       CoeffsInds = 0

       ! Calculate the order of the new coarse nodes
       Cdofs = 0
       DO ind = 1,Fdofs
         IF(CF(ind) > 0) THEN
           Cdofs = Cdofs + 1
           CF(ind) = Cdofs
         END IF
       END DO
       
       ! Initialize stuff before allocations
       AllocationsDone = .FALSE.
       ALLOCATE( PRows(Fdofs/Components+1) )
       PRows(1) = 1

       ! Go through the fine dofs and make a projection based on the strongly coupled nodes
       ! The first time only compute the structure of matrix to be allocated
10     inds = 0
       coeffs = 0.0d0      
       posinds = 0
       neginds = 0
       poscoeffs = 0.0
       negcoeffs = 0.0
       entries = 0

       DO ind=Component1,Fdofs,Components

         ! In component mode this corresponds to a physical node with multiple dofs
         node = (ind-Component1) / Components + 1

         Debug = (node == InfoNode) 

         ! For C-nodes use 1-to-1 mapping
         IF(CF(ind) > 0) THEN
           projnodes = 1
           IF(AllocationsDone) THEN
             PCols(PRows(node)) = CF(ind)
             Pvalues(PRows(node)) = 1.0d0            
           END IF
         ELSE
           
           no = 0
           projnodes = 0
           j = 0
           
           IF(DirectInterpolate) THEN
             DO i=FRows(ind),FRows(ind+1)-1
               ci = Fcols(i)

               ! We desire a block-diagonal projector
               IF(MOD(ci,Components) /= MOD(ind,Components)) CYCLE

               IF(Fixed(ci)) CYCLE

               VALUE = FValues(i)
               IF(ABS(VALUE) < 1.0d-50) CYCLE
               no = no + 1
               inds(no) = ci
               coeffs(no) = VALUE
               IF(ci == ind) THEN
                 diag = no
               ELSE IF(CF(ci) > 0) THEN
                 j = j + 1
               END IF
             END DO
           END IF

           IF(j < DirectLimit) THEN

             IF(no > 0) THEN
               inds(1:no) = 0
               coeffs(1:no) = 0
               no = 0
             END IF

             ! First make the list of the C-neighbouts
             diag = 0
             DO i=FRows(ind),FRows(ind+1)-1
               ci = Fcols(i)

               IF(MOD(ci,Components) /= MOD(ind,Components)) CYCLE

               IF(Fixed(ci)) CYCLE

               IF(CF(ci) > 0 .OR. ci == ind) THEN
                 VALUE = FValues(i)
                 IF(ABS(VALUE) < 1.0d-50) CYCLE 
                 no = no + 1
                 inds(no) = ci
                 coeffs(no) = VALUE
                 CoeffsInds(ci) = no
                 IF(ci == ind) diag = no
               END IF
             END DO

             
             ! Then go though the F-neigbours and express them with linear combinations
             DO i=FRows(ind),FRows(ind+1)-1
               ci = Fcols(i)

               IF(MOD(ci,Components) /= MOD(ind,Components)) CYCLE

               IF(CF(ci) > 0 .OR. ci == ind) CYCLE
 
               IF(Fixed(ci)) CYCLE
             
               DO j=FRows(ci),FRows(ci+1)-1
                 cj = Fcols(j)
                 IF(ci == cj) CYCLE

                 IF(MOD(cj,Components) /= MOD(ind,Components)) CYCLE

                 IF(Fixed(cj)) CYCLE

                 VALUE = Fvalues(i) * Fvalues(j) / Fvalues(Fmat % diag(ci))
                 IF(ABS(VALUE) < 1.0d-50) CYCLE
                 
                 k = CoeffsInds(cj) 
                 IF(k == 0) THEN
                   no = no + 1
                   inds(no) = cj
                   k = no
                   CoeffsInds(cj) = no
                 END IF

                 IF(k > FSIZE) THEN
                   PRINT *,'k',k,l,cj,CoeffsInds(cj)
                   CALL Fatal('InterpolateFineToCoarse','There are more neighbours than expected')
                 END IF
                 coeffs(k) = coeffs(k) - VALUE
               END DO
             END DO

             CoeffsInds(inds(1:no)) = 0
           END IF

             
           k = MOD(inds(1),Components)
           IF( ANY(MOD(inds(1:no),Components) /= k)) THEN
             PRINT *,'no',no
             PRINT *,'inds',inds(1:no)
             PRINT *,'mod',MOD(inds(1:no),Components)
           END IF


           IF(Debug) THEN
             PRINT *,'ind no diag',ind,no,diag
             PRINT *,'coeffs',coeffs(1:no)
             PRINT *,'inds',inds(1:no)
           END IF

           ! Check for Dirichlet points which should not be projected
           IF(no <= 1) THEN
!             IF(diag == 0) THEN
!               CALL Fatal('InterpolateFineToCoarse','Diagonal seems to be zero!')
!             ELSE
               projnodes =  0
               inds(1:no) = 0
               coeffs(1:no) = 0.0
               GOTO 20 
!             END IF
           END IF

          
           ! Calculate the bond limit and the total sums 
           possum = 0.0
           negsum = 0.0
           negmax = 0.0
           posmax = 0.0
           posi = 0
           negi = 0

           IF(Lumping) THEN
             coeffs(diag) = coeffs(diag) - SUM(coeffs(1:no))
           END IF
          
           ! If the diagonal is negative then invert the selection 
           IF(coeffs(diag) < 0.0) THEN
             coeffs(1:no) = -coeffs(1:no)
           END IF

           DO i=1,no
             VALUE = coeffs(i)
             IF(i == diag) THEN
               diagsum = VALUE
             ELSE IF(VALUE > 0.0) THEN
               possum = possum + VALUE
               ci = inds(i)
               IF(CF(ci) > 0) THEN
                 posmax = MAX(posmax, VALUE)
                 posi = posi + 1
               END IF
             ELSE 
               negsum = negsum + VALUE
               ci = inds(i)
               IF(CF(ci) > 0) THEN
                 negmax = MIN(negmax, VALUE)
                 negi = negi + 1
               END IF
             END IF
           END DO

           IF(posi == 0 .AND. negi == 0) THEN
             IF(.TRUE.) PRINT *,'The node is not connected to c-neighbours',ind,inds(1:no)
             projnodes =  0
             inds(1:no) = 0
             coeffs(1:no) = 0.0
             GOTO 20               
           END IF

           posno = 0
           negno = 0

           ! Decide the algorithm knowing which weights dominate +/-
           ! Negative weights dominate
           IF(posi == 0 .OR. ABS(possum) <= ABS(negsum) * FavorNeg) THEN

             IF(negi == 0) THEN
               PRINT *,'Negatively bonded node has no c-neighbours',ind,negsum,possum,posi
               PRINT *,'inds',inds(1:no)
               PRINT *,'cf',CF(inds(1:no))
               PRINT *,'coeffs',coeffs(1:no)
               
               projnodes =  0
               inds(1:no) = 0
               coeffs(1:no) = 0.0
               GOTO 20                
             END IF

             negcsum = 0.0
             DO i=1,no
               VALUE = coeffs(i)
               ci = inds(i)
               IF(CF(ci) == 0) CYCLE
               IF(VALUE < ProjLim * negmax) THEN
                 negno = negno + 1
                 neginds(negno) = ci
                 negcoeffs(negno) = VALUE
                 negcsum = negcsum + VALUE
               ELSE IF(VALUE > ProjLim * posmax) THEN
                 posno = posno + 1
                 posinds(posno) = ci
                 poscoeffs(posno) = VALUE                 
               END IF
             END DO

             negi = negno             
             posi = 0
             refbond = -ProjLim * negsum * negmax / (FavorNeg * negcsum)

             ! Add possible positive weights
             IF(possum > refbond) THEN
               ! Order the positive coefficients in an decreasing order
               DO j = 1, posno-1
                 DO i = 1, posno-1
                   IF(poscoeffs(i) < poscoeffs(i+1)) THEN
                     poscoeffs(posno+1) = poscoeffs(i)
                     poscoeffs(i) = poscoeffs(i+1)
                     poscoeffs(i+1) = poscoeffs(posno+1)
                     posinds(posno+1) = posinds(i)
                     posinds(i) = posinds(i+1)
                     posinds(i+1) = posinds(posno+1)
                   END IF
                 END DO
               END DO
               IF(Debug) THEN
                 IF(posno > 0) THEN
                   PRINT *,'pos connections',posno
                   PRINT *,'inds',posinds(1:posno)
                   PRINT *,'coeffs',poscoeffs(1:posno)
                 END IF
               END IF
               
               poscsum = 0.0

               ! Now go through the possible positive connections 
               DO i=1,posno
                 IF(i == 1) THEN
                   posbond = possum 
                   IF(posbond < refbond) EXIT
                 ELSE
                   posbond = possum * poscoeffs(i) / (poscsum + poscoeffs(i))                
                   IF(posbond < refbond .AND. poscoeffs(i) < 0.99 * poscoeffs(i-1)) EXIT 
                 END IF
                 posi = i
                 poscsum = poscsum + poscoeffs(posi)
               END DO
             END IF

             
           ELSE ! Positive weight dominate

             IF(posi == 0) THEN
               PRINT *,'Positively bonded node has no positive c-neighbours',ind,negsum,possum
               projnodes =  0
               inds(1:no) = 0
               coeffs(1:no) = 0.0
               GOTO 20                
             END IF

             poscsum = 0.0
             DO i=1,no
               VALUE = coeffs(i)
               ci = inds(i)
               IF(CF(ci) == 0) CYCLE
               IF( VALUE > ProjLim * posmax ) THEN
                 posno = posno + 1
                 posinds(posno) = ci
                 poscoeffs(posno) = VALUE
                 poscsum = poscsum + VALUE
               ELSE IF(VALUE < ProjLim * negmax) THEN
                 negno = negno + 1
                 neginds(negno) = ci
                 negcoeffs(negno) = VALUE
               END IF
             END DO

             posi = posno
             negi = 0             
             refbond = ProjLim * possum * posmax / ( poscsum * FavorNeg )

             IF(-negsum > refbond) THEN
               ! Order the negative coefficients in an increasing order
               DO j = 1, negno-1
                 DO i = 1,negno-1
                   IF(negcoeffs(i) > negcoeffs(i+1)) THEN
                     negcoeffs(negno+1) = negcoeffs(i)
                     negcoeffs(i) = negcoeffs(i+1)
                     negcoeffs(i+1) = negcoeffs(negno+1)
                     neginds(negno+1) = neginds(i)
                     neginds(i) = neginds(i+1)
                     neginds(i+1) = neginds(negno+1)
                   END IF
                 END DO
               END DO
               IF(Debug) THEN
                 IF(negno > 0) THEN
                   PRINT *,'neg connections after',negno
                   PRINT *,'inds',neginds(1:negno)
                   PRINT *,'coeffs',negcoeffs(1:negno)
                 END IF
               END IF

               negcsum = 0.0

               ! Now go through the possible positive connections 
               DO i=1,negno
                 IF(i == 1) THEN
                   negbond = negsum 
                   IF(-negbond < refbond) EXIT
                 ELSE
                   negbond = negsum * negcoeffs(i) / (negcsum + negcoeffs(i) )
                   IF(-negbond < refbond .AND. negcoeffs(i) > 0.99 * negcoeffs(i-1)) EXIT 
                 END IF
                 negi = i
                 negcsum = negcsum + negcoeffs(i)
               END DO

             END IF
           END IF

           projnodes = posi + negi
           IF(debug) PRINT *,'bonds',posi,negi,negcsum,neginds(1)
           
           ! Compute the weights and store them to the projection matrix
           IF(AllocationsDone) THEN
             wsum = 0.0

             IF(posi == 0) THEN
               diagsum = diagsum + possum / FavorNeg
             END IF    

             IF(negi == 0) THEN
               diagsum = (diagsum + negsum) / FavorNeg
             END IF            

             DO i=1,negi
               VALUE = -negsum * negcoeffs(i) / (negcsum * diagsum)
               ci = neginds(i)
               IF(Debug) PRINT *,'F-: Pij',ind,CF(ci),VALUE
               PCols(Prows(node)+i-1) = CF(ci)
               PValues(Prows(node)+i-1) = VALUE
               wsum = wsum + VALUE
             END DO

             DO i=1,posi
               VALUE = -possum * poscoeffs(i) / (poscsum * diagsum * FavorNeg)
               ci = posinds(i)
               IF(Debug) PRINT *,'F+: Pij',ind,CF(ci),VALUE
               PCols(Prows(node)+i+negi-1) = CF(ci)
               PValues(Prows(node)+i+negi-1) = VALUE
               wsum = wsum + VALUE
             END DO

             IF(Debug) PRINT *,'ind wsum projnodes',ind,wsum,projnodes
           END IF

           inds(1:no) = 0
           coeffs(1:no) = 0.0
         END IF

20       IF(Debug) PRINT *,'ind nodes',ind,projnodes
         PRows(node+1) = PRows(node) + projnodes
         entries = entries + projnodes

       END DO

       ! Allocate space for the projection matrix and thereafter do the loop again
       IF(.NOT. AllocationsDone) THEN

         Projector => AllocateMatrix()
         no = Fdofs / Components

         Projector % NumberOfRows = no
         ALLOCATE( PCols(entries), PValues(entries) )           

         WRITE(Message,'(A,I8)') 'Projector matrix size',no
         CALL Info('InterpolateF2C',Message)

         WRITE(Message,'(A,F8.2)') 'Projector nodes in average',1.0*entries/no
         CALL Info('InterpolateF2C',Message)

         AllocationsDone = .TRUE.
         PCols   = 0
         PValues = 0
         
         Projector % Rows   => PRows
         Projector % Cols   => PCols 
         Projector % Values => PValues
         
         GOTO 10
       END IF

       IF(Debug) THEN
         no = MAXVAL(Pcols) 
         IF( no /= Cdofs) PRINT *,'********* MAXVAL(Pcols)',no,Cdofs
         
         no = MINVAL(PCols) 
         IF( no < 1) PRINT *,'********** MINVAL(Pcols)',no      
         
         PRINT *,'Projector interval',MINVAL(PValues),MAXVAL(Pvalues)
       END IF

       DEALLOCATE(CoeffsInds)
       
     END FUNCTION InterpolateF2C
!-------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!>     A pseudo geometric version of the previous
!>     Here the strength of connections is assumed to be inversily proportional
!>     to the distance between nodes. 
!------------------------------------------------------------------------------
     FUNCTION InterpolateF2CDistance( Fmat, CF, DOFs) RESULT (Projector)
!------------------------------------------------------------------------------
       USE Interpolation
       USE CRSMatrix
       USE CoordinateSystems
!-------------------------------------------------------------------------------
       TYPE(Matrix_t), TARGET  :: Fmat
       INTEGER, POINTER :: CF(:)
       INTEGER :: DOFs
       TYPE(Matrix_t), POINTER :: Projector
!------------------------------------------------------------------------------
       INTEGER, PARAMETER :: FSIZE=1000, CSIZE=100
       INTEGER :: i, j, k, l, Fdofs, Cdofs, ind, ci, cj, Component1, Components, node
       REAL(KIND=dp), POINTER :: PValues(:), FValues(:)
       INTEGER, POINTER :: FRows(:), FCols(:), PRows(:), PCols(:), CoeffsInds(:)
       REAL(KIND=dp) :: bond, VALUE, possum, negsum, poscsum, negcsum, diagsum, &
           ProjLim, negbond, posbond, maxbond
       LOGICAL :: Debug, AllocationsDone, DirectInterpolate, Lumping
       INTEGER :: inds(FSIZE), posinds(CSIZE), neginds(CSIZE), no, diag, InfoNode, &
           posno, negno, negi, posi, projnodes, DirectLimit, InfoLevel
       REAL(KIND=dp) :: coeffs(FSIZE), poscoeffs(CSIZE), negcoeffs(CSIZE), wsum, &
           refbond, posmax, negmax, x0, y0, z0, x1, y1, z1, s2, Pow

       Debug = .FALSE.
       Components = 1
       IF(Debug) CALL Info('InterpolateF2CDistance','Starting interpolation')
       
       ProjLim = ListGetConstReal(Params,'MG Projection Limit',GotIt)
       IF(.NOT. GotIt) ProjLim = 0.5

       Pow = ListGetConstReal(Params,'MG Geometric Power',GotIt)
       IF(.NOT. GotIt) Pow = 1.0d0

       DirectInterpolate = ListGetLogical(Params,'MG Direct Interpolate',GotIt)
       DirectLimit = ListGetInteger(Params,'MG Direct Interpolate Limit',GotIt)      
       IF(.NOT. GotIt) THEN
         IF(DirectInterpolate) THEN
           DirectLimit = 0
         ELSE
           DirectLimit = HUGE(DirectLimit)
         END IF
       END IF

       Component1 = ListGetInteger(Params,'MG Determining Component',GotIt,minv=1,maxv=DOFs)
       IF(.NOT. GotIt) Component1 = 1

       InfoNode = ListGetInteger(Params,'MG Info Node',GotIt)
       InfoLevel = ListGetInteger(Params,'MG Info Level',GotIt)


       Fdofs = Fmat % NumberOfRows
       FRows   => Fmat % Rows
       FCols   => Fmat % Cols
       FValues => Fmat % Values
       
       ALLOCATE(CoeffsInds(Fdofs))
       CoeffsInds = 0
 
       ! Calculate the order of the new dofs
       Cdofs = 0
       DO ind = 1,Fdofs
         IF(CF(ind) > 0) THEN
           Cdofs = Cdofs + 1
           CF(ind) = Cdofs
         END IF
       END DO
       
       ! Initialize stuff before allocations
       AllocationsDone = .FALSE.
       ALLOCATE( PRows(Fdofs+1) )
       PRows(1) = 1


       ! Go through the fine dofs and make a projection based on the strongly coupled nodes
       ! The first time only compute the structure of matrix to be allocated
10     inds = 0
       coeffs = 0.0d0      
       posinds = 0
       poscoeffs = 0.0

       DO ind=Component1,Fdofs,Components

         node = (ind-Component1) / Components + 1

         Debug = (node == InfoNode .OR. Level == InfoLevel) 
         IF(Debug) PRINT *,'Debug is on',ind,CF(ind)

         ! For C-nodes use 1-to-1 mapping

         IF(CF(ind) > 0) THEN

           projnodes = 1
           IF(AllocationsDone) THEN
             PCols(PRows(node)) = CF(ind)
             Pvalues(PRows(node)) = 1.0d0            
           END IF

         ELSE           
           no = 0
           posmax = 0.0d0
           projnodes = 0

           IF ( Level == Solver % MultiGridTotal ) THEN
             k = ind
           ELSE
             k = AMG(Level+1) % InvCF(ind)
           END IF

           x0 = Mesh % Nodes % x(k)
           y0 = Mesh % Nodes % y(k)
           z0 = Mesh % Nodes % z(k)
                      
           DO i=FRows(ind),FRows(ind+1)-1
             ci = Fcols(i)

             IF(ci == ind) CYCLE             
             IF(MOD(ci-Component1,Components) /= 0) CYCLE
             IF(CF(ci) <= 0) CYCLE

             IF ( Level == Solver % MultiGridTotal ) THEN
               k = ci
             ELSE
               k = AMG(Level+1) % InvCF(ci)
             END IF

             x1 = Mesh % Nodes % x(k)
             y1 = Mesh % Nodes % y(k)
             z1 = Mesh % Nodes % z(k)

             s2 = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0)
             VALUE = s2 ** (-Pow/2.0_dp)

             posmax = MAX(posmax,VALUE)
             no = no + 1

             inds(no) = ci
             coeffs(no) = VALUE
           END DO


           IF(no < DirectLimit) THEN
             DO i=FRows(ind),FRows(ind+1)-1
               ci = Fcols(i)
               
               IF(ci == ind) CYCLE             
               IF(MOD(ci-Component1,Components) /= 0) CYCLE

               DO j=FRows(ci),FRows(ci+1)-1
                 cj = Fcols(j)
                 
                 IF(cj == ind) CYCLE             
                 IF(MOD(cj-Component1,Components) /= 0) CYCLE
                 IF(CF(cj) <= 0) CYCLE
                
                 IF ( Level == Solver % MultiGridTotal ) THEN
                   k = cj
                 ELSE
                   k = AMG(Level+1) % InvCF(cj)
                 END IF
                 x1 = Mesh % Nodes % x(k)
                 y1 = Mesh % Nodes % y(k)
                 z1 = Mesh % Nodes % z(k)

                 s2 = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0)
                 VALUE = s2 ** (-Pow/2.0_dp)

                 posmax = MAX(posmax,VALUE)
                 no = no + 1

                 inds(no) = cj
                 coeffs(no) = VALUE
               END DO
             END DO

             ! Eliminate redundancy in the list of course neighbours
             CALL SortF( no, inds, coeffs) 
             DO i=1, no-1
               IF( inds(i+1) == inds(i)) coeffs(i+1) = 0.0d0
             END DO             

           END IF

           IF(Debug) THEN
             PRINT *,'coeffs',coeffs(1:no)
             PRINT *,'inds',inds(1:no)
           END IF

           IF(no == 0) THEN
             projnodes =  0
             GOTO 20 
           END IF
         
           possum = 0.0d0
           DO i=1,no
             VALUE = coeffs(i)
             IF(VALUE > ProjLim * posmax) THEN
               projnodes = projnodes + 1
               posinds(projnodes) = inds(i)
               poscoeffs(projnodes) = VALUE
               possum = possum + VALUE
             END IF
           END DO
           
           IF(AllocationsDone) THEN
             DO i=1,projnodes
               VALUE =  poscoeffs(i) / possum 
               ci = posinds(i)
               IF(Debug) PRINT *,'Pij',ind,CF(ci),VALUE
               PCols(Prows(node)+i-1) = CF(ci)
               PValues(Prows(node)+i-1) = VALUE
             END DO
           END IF
         END IF

20       IF(Debug) PRINT *,'ind nodes',ind,projnodes
         PRows(node+1) = PRows(node) + projnodes

       END DO

       ! Allocate space for the projection matrix and thereafter do the loop again
       IF(.NOT. AllocationsDone) THEN
         WRITE(Message,'(A,I8)') 'Size of the projector',PRows(Fdofs+1)-1
         CALL Info('InterpolateF2CDistance',Message)

         Projector => AllocateMatrix()
         Projector % NumberOfRows = Fdofs/Components

         ALLOCATE( PCols(PRows(Fdofs/Components+1)-1), PValues(PRows(Fdofs/Components+1)-1) )           
         AllocationsDone = .TRUE.
         PCols   = 0
         PValues = 0
         
         Projector % Rows   => PRows
         Projector % Cols   => PCols 
         Projector % Values => PValues
         
         GOTO 10
       END IF

       DEALLOCATE(CoeffsInds)
       
     END FUNCTION InterpolateF2CDistance
!-------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!>     As the previous one but expects complex valued equation.
!>     The projector is built using only the absolute values. 
!------------------------------------------------------------------------------
     FUNCTION ComplexInterpolateF2C( Fmat, CF ) RESULT (Projector)
!------------------------------------------------------------------------------
       USE Interpolation
       USE CRSMatrix
       USE CoordinateSystems
!-------------------------------------------------------------------------------
       TYPE(Matrix_t), TARGET  :: Fmat
       INTEGER, POINTER :: CF(:)
       TYPE(Matrix_t), POINTER :: Projector
!------------------------------------------------------------------------------
       INTEGER, PARAMETER :: FSIZE=1000, CSIZE=100
       INTEGER :: i, j, k, l, Fdofs, Cdofs, ind, ci, cj, node
       REAL(KIND=dp), POINTER :: PValues(:), FValues(:)
       INTEGER, POINTER :: FRows(:), FCols(:), PRows(:), PCols(:), CoeffsInds(:)
       REAL(KIND=dp) :: bond, ProjLim, posbond, maxbond
       LOGICAL :: Debug, AllocationsDone, DirectInterpolate, Lumping
       INTEGER :: inds(FSIZE), posinds(CSIZE), no, diag, InfoNode, posi, &
           DirectLimit, projnodes
       REAL(KIND=dp) :: wsum, refbond, posmax

       REAL(KIND=dp) :: poscoeffs(CSIZE), VALUE, poscsum, possum, diagsum
       COMPLEX(KIND=dp) :: coeffs(FSIZE), cvalue 

       Debug = .FALSE.

       CALL Info('ComplexInterpolateF2C','Starting interpolation')

       
       ProjLim = ListGetConstReal(Params,'MG Projection Limit',GotIt)
       IF(.NOT. GotIt) ProjLim = 0.5

       Lumping = ListGetLogical(Params,'MG Projection Lumping',GotIt)

       DirectInterpolate = ListGetLogical(Params,'MG Direct Interpolate',GotIt)

       DirectLimit = ListGetInteger(Params,'MG Direct Interpolate Limit',GotIt)      
       IF(.NOT. GotIt) THEN
         IF(DirectInterpolate) THEN
           DirectLimit = 0
         ELSE
           DirectLimit = HUGE(DirectLimit)
         END IF
       END IF

       InfoNode = ListGetInteger(Params,'MG Info Node',GotIt)

       Fdofs = Fmat % NumberOfRows 
       FRows   => Fmat % Rows
       FCols   => Fmat % Cols
       FValues => Fmat % Values
       
       ALLOCATE(CoeffsInds(Fdofs))
       CoeffsInds = 0

       ! Calculate the order of the new dofs 
       Cdofs = 0
       DO ind = 1,Fdofs,2
         IF(CF(ind) > 0) THEN
           Cdofs = Cdofs + 1
           CF(ind) = Cdofs
         END IF
       END DO
       
       IF(Debug) PRINT *,'cdofs',cdofs,'fdofs',fdofs

       ! Initialize stuff before allocations
       AllocationsDone = .FALSE.
       ALLOCATE( PRows(Fdofs/2+1) )
       PRows(1) = 1

       ! Go through the fine dofs and make a projection based on the strongly coupled nodes
       ! The first time only compute the structure of matrix to be allocated
10     inds = 0
       coeffs = 0.0
       posinds = 0
       poscoeffs = 0.0

       DO node=1,Fdofs/2

         ind = 2*node-1

!         Debug = (node == InfoNode) 

         IF(debug) PRINT *,'ind',ind

         ! For C-nodes use 1-to-1 mapping (for real and complex dofs)
         IF(CF(ind) > 0) THEN
           projnodes = 1
           IF(AllocationsDone) THEN
             PCols(PRows(node)) = CF(ind)
             Pvalues(PRows(node)) = 1.0d0            
           END IF
         ELSE
           
           no = 0
           projnodes = 0
           j = 0
           
           IF(DirectInterpolate) THEN
             DO i=FRows(ind),FRows(ind+1)-1,2
               ci = Fcols(i)

               cvalue = CMPLX(FValues(i),-FValues(i+1),KIND=dp)
               IF(ABS(cvalue) < 1.0d-50) CYCLE
               no = no + 1
               inds(no) = FCols(i)
               coeffs(no) = cvalue
               IF(ci == ind) THEN
                 diag = no
               ELSE IF(CF(ci) > 0) THEN
                 j = j + 1
               END IF
             END DO
           END IF

           IF(j < DirectLimit) THEN

             IF(no > 0) THEN
               inds(1:no) = 0
               coeffs(1:no) = 0.0
               no = 0
             END IF

             ! First make the list of the C-neighbouts
             diag = 0
             DO i=FRows(ind),FRows(ind+1)-1,2
               ci = Fcols(i)

               IF(CF(ci) > 0 .OR. ci == ind) THEN
                 cvalue = CMPLX(FValues(i),-Fvalues(i+1),KIND=dp)
                 IF(ABS(cvalue) < 1.0d-50) CYCLE 
                 no = no + 1
                 inds(no) = ci
                 coeffs(no) = cvalue
                 CoeffsInds(ci) = no
                 IF(ci == ind) diag = no
               END IF
             END DO

             
             ! Then go though the F-neigbours and express them with linear combinations
             DO i=FRows(ind),FRows(ind+1)-1,2
               ci = Fcols(i)

               IF(CF(ci) > 0 .OR. ci == ind) CYCLE
               
               DO j=FRows(ci),FRows(ci+1)-1,2
                 cj = Fcols(j)
                 IF(ci == cj) CYCLE
                                  
                 cvalue = CMPLX(Fvalues(i), -Fvalues(i+1),KIND=dp) * &
                          CMPLX(Fvalues(j), -Fvalues(j+1),KIND=dp) / &
                          CMPLX(Fvalues(Fmat % diag(ci)),-Fvalues(Fmat % diag(ci)+1),KIND=dp)
                 IF(ABS(cvalue) < 1.0d-50) CYCLE
                 
                 k = CoeffsInds(cj) 
                 IF(k == 0) THEN
                   no = no + 1
                   inds(no) = cj
                   k = no
                   CoeffsInds(cj) = no
                 END IF

                 IF(k > FSIZE) THEN
                   PRINT *,'k',k,l,cj,CoeffsInds(cj)
                   CALL Fatal('InterpolateFineToCoarse','There are more neighbours than expected')
                 END IF
                 coeffs(k) = coeffs(k) - cvalue
               END DO
             END DO
             
!             IF(debug) PRINT *,'no',no,'inds',inds(1:no)
           
             CoeffsInds(inds(1:no)) = 0
           END IF
             
!          IF(Debug) THEN
!             PRINT *,'ind no diag',ind,no,diag
!             PRINT *,'coeffs',coeffs(1:no)
!             PRINT *,'inds',inds(1:no)
!           END IF

           ! Check for Dirichlet points which should not be projected
           IF(no <= 1) THEN
             IF(diag == 0) THEN
               CALL Fatal('InterpolateFineToCoarse','Diagonal seems to be zero!')
             ELSE
               projnodes =  0
               inds(1:no) = 0
               coeffs(1:no) = 0.0
               GOTO 20 
             END IF
           END IF

          
           ! Calculate the bond limit and the total sums 
           possum = 0.0
           posmax = 0.0
           posi = 0

           IF(Lumping) THEN
             coeffs(diag) = coeffs(diag) - SUM(coeffs(1:no))
           END IF

           DO i=1,no
             VALUE = ABS(coeffs(i))
             IF(i == diag) THEN
               diagsum = VALUE
             ELSE 
               possum = possum + VALUE
               ci = inds(i)
               IF(CF(ci) > 0) THEN
                 posmax = MAX(posmax, VALUE )
                 posi = posi + 1
               END IF
             END IF
           END DO

           IF(posi == 0) THEN
             PRINT *,'The node is not connected to c-neighbours' 
             PRINT *,'inds',inds(1:no)
             PRINT *,'cf',CF(inds(1:no))
             PRINT *,'coeffs',coeffs(1:no)

             projnodes =  0
             inds(1:no) = 0
             coeffs(1:no) = 0.0
             GOTO 20               
           END IF

           projnodes = 0
           poscsum = 0.0

           DO i=1,no
             VALUE = ABS(coeffs(i))
             ci = inds(i)
             IF(CF(ci) == 0) CYCLE
             IF(ABS(VALUE) > ProjLim * posmax) THEN
               projnodes = projnodes + 1
               posinds(projnodes) = ci
               poscoeffs(projnodes) = VALUE                 
               poscsum = poscsum + VALUE
             END IF
           END DO

!          IF(debug) PRINT *,'bonds',projnodes,poscsum
           
           wsum = 0.0
          ! Compute the weights and store them to a projection matrix
           IF(AllocationsDone) THEN

            wsum = 0.0
            DO i=1,projnodes
               VALUE = possum * poscoeffs(i) / (poscsum * diagsum)
               ci = posinds(i)
               IF(Debug) PRINT *,'F+: Pij',node,CF(ci),VALUE,poscoeffs(i)

               PCols(Prows(node)+i-1) = CF(ci) 
               PValues(Prows(node)+i-1) = VALUE
               
               wsum = wsum + VALUE
             END DO
             IF(Debug) PRINT *,'ind projnodes wsum',ind,projnodes,wsum,diagsum,poscsum,possum

           END IF

           inds(1:no) = 0
           coeffs(1:no) = 0.0
         END IF

20       PRows(node+1) = PRows(node) + projnodes
!         IF(Debug) PRINT *,'ind nodes',ind,projnodes,wsum

       END DO

       ! Allocate space for the projection matrix and thereafter do the loop again
       IF(.NOT. AllocationsDone) THEN
         IF(Debug) PRINT *,'allocated space',PRows(Fdofs/2+1)-1
         Projector => AllocateMatrix()
         Projector % NumberOfRows = Fdofs/2

         ALLOCATE( PCols(PRows(Fdofs/2+1)-1), PValues(PRows(Fdofs/2+1)-1) )           
         AllocationsDone = .TRUE.
         PCols   = 0
         PValues = 0
         
         Projector % Rows   => PRows
         Projector % Cols   => PCols 
         Projector % Values => PValues
         
         GOTO 10
       END IF

       DEALLOCATE(CoeffsInds)

       Cdofs = 0
       DO ind = 1,Fdofs
         IF(CF(ind) > 0) THEN
           Cdofs = Cdofs + 1
           CF(ind) = Cdofs
         END IF
       END DO
              
     END FUNCTION ComplexInterpolateF2C
!-------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   This subroutine may be used to check that mapping Xf = P Xc is accurate.
!>   The list of original nodes for each level is stored in the Cnodes vector.
!>   For the momont this only works for cases with one DOF and no permutation!
!------------------------------------------------------------------------------
    SUBROUTINE AMGTest(direction) 
!------------------------------------------------------------------------------
      INTEGER :: direction

      INTEGER :: i,j,Rounds, nods1, nods2, SaveLimit
      LOGICAL :: GotIt
      CHARACTER(LEN=MAX_NAME_LEN) :: Filename
      REAL(KIND=dp) :: RNorm
      REAL(KIND=dp), POINTER :: Ina(:), Inb(:), Outa(:), Outb(:)

      nods1 = Matrix1 % NumberOfRows
      nods2 = Matrix2 % NumberOfRows

      SaveLimit = ListGetInteger(Params,'MG Coarse Nodes Save Limit',GotIt) 
      IF(.NOT. GotIt) SaveLimit = nods2

      IF(nods2 > SaveLimit) RETURN

!      Project the fine dofs to the coarse dofs
!      ----------------------------------------

      IF(Direction == 0) THEN
        WRITE( Filename,'(a,i1,a,i1,a)' ) 'mapping', Level,'to',Level-1, '.dat'

        ALLOCATE( Ina(nods1), Inb(nods1), Outa(nods2), Outb(nods2) )

        IF(nods1 < SaveLimit) THEN
          IF ( Level == Solver % MultiGridTotal ) THEN
            WRITE( Filename,'(a,i1,a)' ) 'nodes', Solver % MultiGridTotal - Level,'.dat'
            OPEN (10,FILE=Filename)        
            DO i=1,nods1
              WRITE (10,'(3ES17.8E3)') Mesh % Nodes % X(i), Mesh % Nodes % Y(i), Mesh % Nodes % Z(i)
            END DO
          END IF
        END IF

        WRITE( Filename,'(a,i1,a)' ) 'nodes', Solver % MultiGridTotal - Level+1,'.dat'
        OPEN (10,FILE=Filename)        
        DO i=1,nods2
          WRITE (10,'(3ES17.8E3)') Mesh % Nodes % X(AMG(Level) % InvCF(i)), &
              Mesh % Nodes % Y(AMG(Level) % InvCF(i)), Mesh % Nodes % Z(AMG(Level) % InvCF(i))
        END DO
        CLOSE(10)
      END IF
  

      IF(Direction == 1) THEN
        WRITE( Filename,'(a,i1,a,i1,a)' ) 'mapping', Level,'to',Level-1, '.dat'

        ALLOCATE( Ina(nods1), Inb(nods1), Outa(nods2), Outb(nods2) )

        IF ( Level == Solver % MultiGridLevel ) THEN
          Ina = Mesh % Nodes % X
          Inb = Mesh % Nodes % Y
        ELSE
          Ina = Mesh % Nodes % X(AMG(Level+1) % InvCF )
          Inb = Mesh % Nodes % Y(AMG(Level+1) % InvCF )
        END IF

        Outa = 0.0d0
        Outb = 0.0d0
        
        CALL CRS_ProjectVector( ProjPN, Ina, Outa, 1, Trans = .FALSE. )
        CALL CRS_ProjectVector( ProjPN, Inb, Outb, 1, Trans = .FALSE. )        

        OPEN (10,FILE=Filename)        
        DO i=1,nods2
          WRITE (10,'(4ES17.8E3)') Mesh % Nodes % X(AMG(Level) % InvCF(i) ), &
              Mesh % Nodes % Y(AMG(Level) % InvCF(i) ) , Outa(i), Outb(i)
        END DO
        CLOSE(10)
      END IF
      
!      Project the coarse dofs to the fine dofs
!      ----------------------------------------

      IF(Direction == 2) THEN

        WRITE( Filename,'(a,i1,a,i1,a)' ) 'mapping', Level-1,'to',Level, '.dat'

        ALLOCATE( Ina(nods2), Inb(nods2), Outa(nods1), Outb(nods1) )
      
        Ina = 0.0d0
        Inb = 0.0d0

        Ina = Mesh % Nodes % X(AMG(Level) % InvCF )
        Inb = Mesh % Nodes % Y(AMG(Level) % InvCF )

        Outa = 0.0d0
        Outb = 0.0d0

        
        PRINT *,'Initial Interval x',MINVAL(Ina),MAXVAL(Ina)
        PRINT *,'Initial Interval y',MINVAL(Inb),MAXVAL(Inb)
        PRINT *,'Initial Mean Values',SUM(Ina)/SIZE(Ina),SUM(Inb)/SIZE(Inb)

        CALL CRS_ProjectVector( ProjQT, Ina, Outa, 1, Trans = .TRUE. )
        CALL CRS_ProjectVector( ProjQT, Inb, Outb, 1, Trans = .TRUE. )        

        PRINT *,'Final Interval x',MINVAL(Outa),MAXVAL(Outa),SUM(Outa)/SIZE(Outa)
        PRINT *,'Final Interval y',MINVAL(Outb),MAXVAL(Outb),SUM(Outb)/SIZE(Outb)
        PRINT *,'Final Mean Values',SUM(Outa)/SIZE(Outa),SUM(Outb)/SIZE(Outb)
        
        OPEN (10,FILE=Filename)        
 
        IF ( Level == Solver % MultiGridLevel ) THEN
          DO i=1,nods1
            WRITE (10,'(4ES17.8E3)') Mesh % Nodes % X(i), Mesh % Nodes % Y(i) , Outa(i), Outb(i)
          END DO
        ELSE
          DO i=1,nods1
            WRITE (10,'(4ES17.8E3)') Mesh % Nodes % X(AMG(Level+1) % InvCF(i) ), &
                Mesh % Nodes % Y(AMG(Level+1) % InvCF(i) ) , Outa(i), Outb(i)
          END DO
        END IF
        CLOSE(10)        

      END IF

      WRITE(Message,'(A,A)') 'Save mapped results into file: ',TRIM(Filename)
      CALL Info('MGTest',Message)

      DEALLOCATE(Ina, Inb, Outa, Outb)

    END SUBROUTINE AMGTest

!------------------------------------------------------------------------------



!----------------------------------------------------------------------------
!>   The following subroutines are CR (Compatible Relaxation) versions of the
!>   simple relaxation schemes. The relaxation speed of equation Ax = 0 to x=0
!>   may be used as a measure of the goodness of the chosen coarse node set.
!-----------------------------------------------------------------------------
    SUBROUTINE CR_Jacobi( A, x0, x1, f, Rounds )
!-----------------------------------------------------------------------------
       TYPE(Matrix_t), POINTER :: A
       INTEGER :: Rounds, f(:)
       REAL(KIND=dp) :: x0(:), x1(:), s
       INTEGER :: i,j,k,n
       INTEGER, POINTER :: Rows(:), Cols(:)
       REAL(KIND=dp), POINTER :: Values(:)
       
       n = A % NumberOfRows
       Rows   => A % Rows
       Cols   => A % Cols
       Values => A % Values

       DO k=1,Rounds
         IF(k > 1) x0 = x1
         DO i=1,n
           IF(f(i) > 0) CYCLE
           s = 0.0d0
           DO j=Rows(i),Rows(i+1)-1
             s = s + x0(Cols(j)) * Values(j)
           END DO

           x1(i) = x0(i) - s / A % Values(A % Diag(i))
         END DO
       END DO
     END SUBROUTINE CR_Jacobi


!------------------------------------------------------------------------------
    SUBROUTINE CR_GS( A, x0, x1, f, Rounds )
!------------------------------------------------------------------------------
       TYPE(Matrix_t), POINTER :: A
       INTEGER :: Rounds, f(:)
       REAL(KIND=dp) :: x0(:), x1(:)
       INTEGER :: i,j,k,n
       REAL(KIND=dp) :: s
       INTEGER, POINTER :: Cols(:),Rows(:)
       REAL(KIND=dp), POINTER :: Values(:)
     
       n = A % NumberOfRows
       Rows   => A % Rows
       Cols   => A % Cols
       Values => A % Values
       
       x1 = x0

       DO k=1,Rounds
         DO i=1,n
           IF(f(i) > 0) CYCLE
           s = 0.0d0
           DO j=Rows(i),Rows(i+1)-1
             s = s + x1(Cols(j)) * Values(j)
           END DO

           x1(i) = x1(i) - s / A % Values(A % Diag(i))
         END DO
         IF(k == Rounds-1) x0 = x1
       END DO

     END SUBROUTINE CR_GS


!------------------------------------------------------------------------------
     SUBROUTINE CR_SGS( A, x0, x1, f, Rounds)
!------------------------------------------------------------------------------
       TYPE(Matrix_t), POINTER :: A
       INTEGER :: Rounds, f(:)
       REAL(KIND=dp) CONTIG :: x0(:),x1(:)
       INTEGER :: i,j,k,n
       REAL(KIND=dp) :: s
       INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
       REAL(KIND=dp), POINTER CONTIG :: Values(:)
     
       n = A % NumberOfRows
       Rows   => A % Rows
       Cols   => A % Cols
       Values => A % Values
       
       x1 = x0

       DO k=1,Rounds
         DO i=1,n
           IF(f(i) > 0) CYCLE
           s = 0.0d0
           DO j=Rows(i),Rows(i+1)-1
             s = s + x1(Cols(j)) * Values(j)
           END DO
           x1(i) = x1(i) - s / A % Values(A % Diag(i))
         END DO

         DO i=n,1,-1
           IF(f(i) > 0) CYCLE
           s = 0.0d0
           DO j=Rows(i),Rows(i+1)-1
             s = s + x1(Cols(j)) * Values(j)
           END DO
           x1(i) = x1(i) - s / A % Values(A % Diag(i))
         END DO

         IF(k == Rounds-1) x0 = x1
       END DO
     END SUBROUTINE CR_SGS


!------------------------------------------------------------------------------
     SUBROUTINE CR_CSGS( n, A, rx0, rx1, f, Rounds)
!------------------------------------------------------------------------------
       TYPE(Matrix_t), POINTER :: A
       INTEGER :: n, Rounds, f(:)
       REAL(KIND=dp) CONTIG :: rx0(:),rx1(:)

       COMPLEX(KIND=dp) :: x0(n/2),x1(n/2), s
       INTEGER :: i,j,k,j2,diag,l
       INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
       REAL(KIND=dp), POINTER CONTIG :: Values(:)
     
       n = A % NumberOfRows
       Rows   => A % Rows
       Cols   => A % Cols
       Values => A % Values
              
       DO i=1,n/2
         x0(i) = CMPLX( rx0(2*i-1), rx0(2*i),KIND=dp )
       END DO

       x1 = x0      

       l = ListGetInteger(Params,'MG Info Node',GotIt)


       DO k=1,Rounds
         DO i=1,n/2
           IF(f(2*i-1) > 0) CYCLE
           IF(f(2*i) > 0) PRINT *,'f(2*i-1) /= f(2*i)',i 

           s = 0.0d0

           ! Go only through the real part of matrix
           DO j=Rows(2*i-1),Rows(2*i)-1             
             j2 = j + Rows(2*i) - Rows(2*i-1)
             
             IF(i==l) PRINT *,'i0',i,j2-j,(Cols(j)-1)/2+1,(Cols(j2)-1)/2+1

             IF(MOD(Cols(j),2) == 0) CYCLE
             s = s + x1((Cols(j)-1)/2+1) * CMPLX( Values(j), Values(j2),KIND=dp)
             IF(i == l) THEN
!               print *,'i',i,2*i-1,Cols(j),Cols(j2),(Cols(j)-1)/2+1
               PRINT *,'a',i,(Cols(j2)-1)/2+1,Values(j),Values(j2)
             END IF
           END DO
           j = A % Diag(2*i-1)
           j2 = j + Rows(2*i) - Rows(2*i-1)          
           x1(i) = x1(i) - s / CMPLX( Values(j), Values(j2),KIND=dp) 

           IF(i == l) THEN
             PRINT *,'diag',i,j,j2,s,CMPLX( Values(j), Values(j2),KIND=dp)
           END IF
         END DO

         DO i=n/2,1,-1
           IF(f(2*i-1) > 0) CYCLE
           s = 0.0d0
           DO j=Rows(2*i-1),Rows(2*i)-1
             IF(MOD(Cols(j),2) == 0) CYCLE
             j2 = j + Rows(2*i) - Rows(2*i-1)
             s = s + x1((Cols(j)-1)/2+1) * CMPLX( Values(j), Values(j2),KIND=dp)
           END DO
           j = A % Diag(2*i-1)
           j2 = j + Rows(2*i) - Rows(2*i-1)          
           x1(i) = x1(i) - s / CMPLX( Values(j), Values(j2),KIND=dp) 
         END DO

         IF(k == Rounds-1) x0 = x1
       END DO

       DO i=1,n/2
         rx0(2*i-1) =  REAL( x0(i) )
         rx0(2*i-0) =  AIMAG( x0(i) )
         rx1(2*i-1) =  REAL( x1(i) )
         rx1(2*i-0) =  AIMAG( x1(i) )
       END DO

     END SUBROUTINE CR_CSGS


!-------------------------------------------------------------------------------
  SUBROUTINE CRS_ProjectVector( PMatrix, u, v, DOFs, Trans )
!-------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: PMatrix
    REAL(KIND=dp), POINTER :: u(:),v(:)
    INTEGER :: DOFs
    LOGICAL, OPTIONAL :: Trans
!-------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,n
    REAL(KIND=dp), POINTER CONTIG :: Values(:)
    LOGICAL :: LTrans
    INTEGER, POINTER CONTIG :: Rows(:), Cols(:)
!-------------------------------------------------------------------------------
    LTrans = .FALSE.
    IF ( PRESENT( Trans ) ) LTrans = Trans

    n = PMatrix % NumberOfRows
    Rows   => PMatrix % Rows
    Cols   => PMatrix % Cols
    Values => PMatrix % Values

    IF(Trans) THEN
      IF(SIZE(u)/DOFS /= n) THEN
        PRINT *,'dofs',dofs,'u',SIZE(u),'n',n,'u/dofs',SIZE(u)/dofs        
        CALL Fatal('CRS_ProjectVector','Incompatible transpose sizes')
      END IF
    ELSE      
      IF(SIZE(v)/DOFS /= n) THEN
        PRINT *,'dofs',dofs,'v',SIZE(v),'n',n,'v/dofs',SIZE(v)/dofs
        CALL Fatal('CRS_ProjectVector','Incompatible sizes')
      END IF
    END IF

    v = 0.0d0

    IF( DOFs == 1) THEN
      IF ( LTrans ) THEN
        DO i=1,n
          DO j=Rows(i),Rows(i+1)-1
            v(Cols(j)) = v(Cols(j)) + u(i) * Values(j)
          END DO
        END DO
      ELSE
        DO i=1,n
          DO j = Rows(i), Rows(i+1)-1
            v(i) = v(i) + u(Cols(j)) * Values(j)
          END DO
        END DO
      END IF
    ELSE
      IF ( LTrans ) THEN
        DO i=1,n
          DO j=Rows(i),Rows(i+1)-1
            DO k=1,DOFs
              v(DOFs*(Cols(j)-1)+k) = &
                  v(DOFs*(Cols(j)-1)+k) + u(DOFs*(i-1)+k) * Values(j)
            END DO
          END DO
        END DO
      ELSE
        DO i=1,n
          DO j = Rows(i), Rows(i+1)-1
            DO k=1,DOFs
              v(DOFs*(i-1)+k) = &
                  v(DOFs*(i-1)+k) + u(DOFs*(Cols(j)-1)+k) * Values(j)
            END DO
          END DO
        END DO
      END IF      
    END IF
!-------------------------------------------------------------------------------
  END SUBROUTINE CRS_ProjectVector
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  SUBROUTINE CRS_ClusterProject( CF, u, v, DOFs, Trans )
!-------------------------------------------------------------------------------
    INTEGER, POINTER :: CF(:)
    REAL(KIND=dp), POINTER CONTIG :: u(:),v(:)
    INTEGER :: DOFs
    LOGICAL, OPTIONAL :: Trans
!-------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,nv,nu
    REAL(KIND=dp), POINTER CONTIG :: Values(:)
    LOGICAL :: LTrans
    INTEGER, POINTER CONTIG :: Rows(:), Cols(:)
!-------------------------------------------------------------------------------
    LTrans = .FALSE.
    IF ( PRESENT( Trans ) ) LTrans = Trans

    nu = SIZE(u)
    nv = SIZE(v)

    v = 0.0d0

    IF(Trans) THEN
      DO i=1,nv
        j = CF(i)
        IF(j > 0) v(i) = v(i) + u(j)
      END DO      
    ELSE
      DO i=1,nu
        j = CF(i)
        IF(j > 0) v(j) = v(j) + u(i)
      END DO
    END IF

!-------------------------------------------------------------------------------
  END SUBROUTINE CRS_ClusterProject
!-------------------------------------------------------------------------------





!------------------------------------------------------------------------------
!>     Project matrix A to B: B = PAR
!------------------------------------------------------------------------------
    SUBROUTINE CRS_ProjectMatrixCreate( A, P, R, B, DOFs ) 
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A,P,R,B
      INTEGER :: DOFs
!------------------------------------------------------------------------------
      INTEGER, POINTER :: L1(:), L2(:)
      REAL(KIND=dp) :: s, Epsilon
      REAL(KIND=dp), POINTER :: R1(:),R2(:)
      INTEGER :: i,ia,j,k,l,m,n,NA,NB,ci,cj,ck,cl,NoRow,i2,j2
      INTEGER :: TotalNonzeros
      LOGICAL :: AllocationsDone, GotIt 
      INTEGER, ALLOCATABLE :: Ind(:), Row(:)
!------------------------------------------------------------------------------

      AllocationsDone = .FALSE.
      NA = A % NumberOfRows
      NB = DOFs * P % NumberOfRows
      ALLOCATE( Row( NA ), Ind( NB ) )


      Epsilon = ListGetConstReal(Params,'MG Matrix Create Epsilon',GotIt)
      IF(.NOT. GotIt) Epsilon = SQRT(TINY(Epsilon))

      B => AllocateMatrix()      
      B % NumberOfRows = P % NumberOfRows * DOFs
      ALLOCATE( B % Rows( B % NumberOfRows + 1 ), &
          B % Diag( B % NumberOfRows ), &
          B % RHS( B % NumberOfRows ) )        
      B % RHS = 0.0d0
      B % Diag = 0
      B % Rows(1) = 1

      Row = 0
      Ind = 0
10    TotalNonzeros = 0

      IF(DOFs == 1) THEN
        DO i=1,P % NumberOfRows
          
          NoRow = 0
          DO j=P % Rows(i),P % Rows(i+1)-1          
            cj = P % Cols(j)
            
            DO k=A % Rows(cj), A % Rows(cj+1)-1
              ck = A % Cols(k) 
              
              DO l=R % Rows(ck), R % Rows(ck+1)-1

                cl = R % Cols(l)                
                i2 = Row(cl)
                
                IF ( i2 == 0) THEN
                  NoRow = NoRow + 1
                  Ind(NoRow) = cl
                  i2 = B % Rows(i) + NoRow - 1                   
                  Row(cl) = i2
                  
                  IF(AllocationsDone) THEN
                    IF(i == cl) B % Diag(cl) = i2
                    B % Cols(i2) = cl                     
                    B % Values(i2) = P % Values(j) * A % Values(k) * R % Values(l)
                  END IF
                ELSE IF(AllocationsDone) THEN
                  B % Values(i2) = B % Values(i2) + P % Values(j) * A % Values(k) * R % Values(l)
                END IF
              END DO
            END DO
          END DO
          
          DO j=1,NoRow
            Row(Ind(j)) = 0
          END DO
          
          B % Rows(i+1) = B % Rows(i) + NoRow
          TotalNonzeros  = TotalNonzeros + NoRow
        END DO

      ELSE  ! DOFs /= 1

        DO i=1,P % NumberOfRows
          DO m = 1,DOFs
            ia = DOFs * (i-1) + m 
            
            NoRow = 0            
            DO j=P % Rows(i),P % Rows(i+1)-1          

              cj = DOFs * (P % Cols(j)-1) + m
              
              DO k=A % Rows(cj), A % Rows(cj+1)-1

                IF( ABS(A % Values(k)) < Epsilon) CYCLE
 
                ck = A % Cols(k)
                n = MOD(ck-1,DOFs)+1
                ck = (ck-1) / DOFs + 1

                DO l=R % Rows(ck), R % Rows(ck+1)-1

                  cl = DOFs * (R % Cols(l)-1) + n
                  i2 = Row(cl)
                    
                  IF ( i2 == 0 ) THEN
                    NoRow = NoRow + 1
                    Ind(NoRow) = cl
                    i2 = B % Rows(ia) + NoRow - 1
                    Row(cl) = i2
                    
                    IF(AllocationsDone) THEN
                      IF(ia == cl) B % diag(cl) = i2 
                      B % Cols(i2) = cl                     
                      B % Values(i2) = P % Values(j) * A % Values(k) * R % Values(l)
                    END IF
                  ELSE
                    IF(AllocationsDone) THEN
                      B % Values(i2) = B % Values(i2) + P % Values(j) * A % Values(k) * R % Values(l)
                    END IF
                  END IF
                END DO
              END DO
            END DO
            
            DO j=1,NoRow
              Row(Ind(j)) = 0
            END DO
            
            B % Rows(ia+1) = B % Rows(ia) + NoRow
            TotalNonzeros  = TotalNonzeros + NoRow
          END DO
        END DO
      END IF

      IF(.NOT. AllocationsDone) THEN
        ALLOCATE( B % Cols( TotalNonzeros ), B % Values( TotalNonzeros ) )
        B % Cols = 0
        B % Values = 0.0d0
        AllocationsDone = .TRUE.
        GOTO 10 
      END IF

      DEALLOCATE( Row, Ind )

!      WRITE( Message,'(A,I8)' ) 'Projected matrix created with size',TotalNonZeros
!      CALL Info( 'CRS_ProjectMatrixCreate', Message )

      WRITE(Message,'(A,F10.3)') 'Coarse matrix reduction factor',&
          1.0 *  SIZE(Matrix1 % Cols) / TotalNonZeros
      CALL Info('CRS_ProjectMatrixCreate',Message)
       
   END SUBROUTINE CRS_ProjectMatrixCreate
!------------------------------------------------------------------------------



  SUBROUTINE SaveMatrix( A, FileName)
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    CHARACTER(LEN=*) :: FileName

    INTEGER :: i,j,k

    PRINT *,'Saving matrix ',TRIM(FileName),' of size ',A % NumberOfRows

    OPEN (10, FILE=FileName) 

    DO i=1,A % NumberOfRows
      DO j=A % Rows(i),A % Rows(i+1)-1
        WRITE(10,*) i,A % Cols(j),A % Values(j)
      END DO
    END DO

    CLOSE(10)

  END SUBROUTINE SaveMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   END SUBROUTINE AMGSolve
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Subroutine containing agglomeration or cluster multigrid solver.
!> This provides in princinple an economical approach to multilevel schemes.
!> The utilization of the routines are still not complete.
! 
!       Author: Peter R�back
!       Modified by: 
!       Date of modification: 30.10.2007
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE CMGSolve( Matrix1, Solution, &
    ForceVector, DOFs, Solver, Level, NewSystem )
!------------------------------------------------------------------------------
    USE ModelDescription
    USE Smoothers

    IMPLICIT NONE
    
    TYPE(Matrix_t), POINTER :: Matrix1
    INTEGER :: DOFs, Level
    LOGICAL, OPTIONAL :: NewSystem
    TYPE(Solver_t), TARGET :: Solver
    REAL(KIND=dp) CONTIG :: ForceVector(:), Solution(:)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER   :: Mesh
    TYPE(Matrix_t), POINTER :: Matrix2, Pmatrix
    TYPE(Solver_t), POINTER :: PSolver
   
    INTEGER :: i,j,k,l,m,n,n2,k1,k2,iter,MaxIter = 100, DirectLimit, &
        MinLevel, OrigSize=0, InvLevel, Sweeps
    LOGICAL :: Condition, Found, Parallel, EliminateDir, CoarseSave, Liter
    CHARACTER(LEN=MAX_NAME_LEN) :: str,str2,str3,FileName,LowestSolver
    INTEGER, POINTER :: CF(:), InvCF(:), Iters(:)
    
    REAL(KIND=dp), ALLOCATABLE, TARGET :: Residual(:),  Solution2(:)
    REAL(KIND=dp), POINTER CONTIG :: TmpArray(:,:), Residual2(:)
    REAL(KIND=dp) :: ResidualNorm, RHSNorm, Tolerance, ILUTOL, Alpha, Rnorm
#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: tt, tmp
#else
    REAL(KIND=dp) :: CPUTime, tt, tmp
#endif
    TYPE(ValueList_t), POINTER :: Params

    LOGICAL :: NewLinearSystem, gotit, Normalize 

    SAVE NewLinearSystem, MinLevel, OrigSize
    
!------------------------------------------------------------------------------

    WRITE(Message,'(A,I2)') 'Starting level ',Level
    CALL Info('CMGSolve',Message,Level=10)

    IF( ParEnv % PEs > 1 ) THEN
      CALL Fatal('CMGSolve','This agglomeration multigrid is not parallel')
    END IF


    Mesh => Solver % Mesh    
    Params => Solver % Values

    tt = CPUTime()
!
!      Initialize:
!      -----------
    Parallel = ParEnv % PEs > 1
    InvLevel = 1 + Solver % MultiGridTotal - Level

    n = Matrix1 % NumberOfRows
    ALLOCATE( Residual(n) )
    Residual = 0.0d0

    ! This does not seem to have any effect?
    !-------------------------------------------------------------------
    Normalize = ListGetLogical( Params,'MG Normalize RHS',GotIt)

    RHSNorm = ParallelReduction(SQRT(SUM(ForceVector**2)))
    IF( Normalize ) THEN
      Solution(1:n) = Solution(1:n) / RHSnorm
      ForceVector(1:n) = ForceVector(1:n) / RHSnorm
    END IF

    ! This is a counter that for the first full resursive round keeps the 
    ! flag NewLinearSystem true.
    IF ( Level == Solver % MultiGridLevel ) THEN
      NewLinearSystem = .TRUE.
      MinLevel = Solver % MultiGridLevel
      IF ( PRESENT( NewSystem ) ) THEN
        NewLinearSystem = NewLinearSystem .AND. NewSystem
      END IF
    END IF

    ! If at lowest level, choose a solution method and go for it
    !------------------------------------------------------------
    DirectLimit = ListGetInteger(Params,'MG Lowest Linear Solver Limit',GotIt) 
    IF(.NOT. GotIt) DirectLimit = 20

    IF ( Level <= 1 .OR. n < DirectLimit) THEN

      NewLinearSystem = .FALSE.
      IF(Level == 1 .AND. PRESENT(NewSystem)) NewSystem = .FALSE.

      CALL ListPushNamespace('mglowest:')

      LowestSolver = ListGetString(Params,'MG Lowest Linear Solver',Found)

      IF ( .NOT. Found ) THEN
        LowestSolver = 'direct'
        LIter = ListGetLogical(Params,'MG Lowest Linear Solver Iterative',Found)
        IF ( .NOT. Found .AND. Parallel ) LIter=.TRUE.
        IF ( LIter ) LowestSolver='iterative'
      END IF
      
      CALL Info('CMGSolve','Solving for the lowest level: '//TRIM(LowestSolver),Level=8)
       
      SELECT CASE( LowestSolver )
        
      CASE('iterative')
        IF ( Parallel ) THEN
          CALL ParallelIter( Matrix1, Matrix1 % ParallelInfo, DOFs, &
              Solution, ForceVector, Solver, Matrix1 % ParMatrix )
        ELSE
          CALL IterSolver( Matrix1, Solution, ForceVector, Solver )
        END IF
        
      CASE('direct')
        CALL DirectSolver( Matrix1, Solution, ForceVector, Solver )
        
      CASE('smoother')
        IF ( Parallel ) THEN
          CALL ParallelInitSolve( Matrix1, Solution, &
              ForceVector, Residual, NewLinearSystem )
        END IF
        PSolver => Solver
        tmp = MGSmooth( PSolver, Matrix1, Solver % Mesh, Solution, &
            ForceVector, Residual, Level, DOFs, LowestSmooth = .TRUE. )
        
      CASE('none') 
        CALL Info('CMGSolve','Applying no solver for coarsest level')
                
      CASE DEFAULT
        CALL Warn( 'CMGSolve', 'Unknown solver for MG lowest level: '//TRIM(LowestSolver) )
        CALL Warn( 'CMGSolve', 'Using iterative solver' )
        
        IF ( Parallel ) THEN
          CALL ParallelIter( Matrix1, Matrix1 % ParallelInfo, DOFs, &
              Solution, ForceVector, Solver, Matrix1 % ParMatrix )
        ELSE
          CALL IterSolver( Matrix1, Solution, ForceVector, Solver )
        END IF
      END SELECT

      IF( Normalize ) THEN
        Solution(1:n) = Solution(1:n) * RHSNorm
        ForceVector(1:n) = ForceVector(1:n) * RHSnorm
      END IF

      DEALLOCATE( Residual ) 
      
      CALL ListPopNamespace('mglowest:')

      CALL Info('CMGSolve','Lowest level solved',Level=9)

      RETURN

    END IF
    

!      Parallel initializations:
!      -------------------------
    IF ( Parallel ) THEN
      CALL ParallelInitSolve( Matrix1, Solution, ForceVector, Residual )
      PMatrix => ParallelMatrix( Matrix1 ) 
    END IF

!      Compute residual:
!      -----------------
    IF(SIZE(Solution) /= Matrix1 % NumberOfRows) THEN
      CALL WARN('CMGSolve','Solution and matrix sizes differ')
    END IF

    CALL MGmv( Matrix1, Solution, Residual, .TRUE. )
    Residual(1:n) = ForceVector(1:n) - Residual(1:n)

    Tolerance = ListGetConstReal( Params,'Linear System Convergence Tolerance' )
    CoarseSave = ListGetLogical(Params,'MG Cluster Save',GotIt)

!---------------------------------------------------------------------
!
!      Initialize the multilevel solve:
!      --------------------------------

    IF( .NOT. ListGetLogical(Params,'MG Recompute Projector',GotIt) ) THEN
      NewLinearSystem = .NOT. ASSOCIATED(Matrix1 % Parent)
    END IF

    IF ( NewLinearSystem ) THEN
      ! If the projection matrix is made again deallocate the old projectors
      IF(ASSOCIATED(Matrix1 % Parent)) CALL FreeMatrix(Matrix1 % Parent)
      IF(ASSOCIATED(Matrix1 % Ematrix)) CALL FreeMatrix(Matrix1 % Ematrix)

      CALL Info('CMGSolve','-------------------------------------------------')
      WRITE( Message, '(A,I3)' ) 'Creating a new matrix and projector for level',Level
      CALL Info('CMGSolve', Message)
      MinLevel = MIN(MinLevel,Level)
      
      EliminateDir = ListGetLogical(Params,'MG Eliminate Dirichlet',GotIt) 
      IF(.NOT. GotIt) EliminateDir = .TRUE.
      IF(Level /= Solver % MultiGridTotal) EliminateDir = .FALSE.
      
      CALL ChooseClusterNodes(Matrix1, Solver, DOFs, EliminateDir, CF)      

!      CALL CRS_InspectMatrix( Matrix1 ) 

      CALL CRS_ClusterMatrixCreate( Matrix1, CF, Matrix2, DOFs)       

!      CALL CRS_InspectMatrix( Matrix2 ) 

      TmpArray => ListGetConstRealArray(Params,'MG Cluster Alpha', gotIt )
      IF(gotIt) THEN
        Alpha = TmpArray(MIN(InvLevel,SIZE(TmpArray,1)),1)
        Matrix2 % Values = Matrix2 % Values / Alpha
      END IF
 
      WRITE( Message, '(A,F8.3)' ) 'MG coarse matrix creation time (s): ', CPUTime() - tt
      CALL Info( 'CMGSolve', Message, Level=5 )

      Matrix2 % Child  => Matrix1
      Matrix1 % Parent => Matrix2
      Matrix1 % Gorder => CF

      ! Make a table showing the clustering
      IF(CoarseSave) THEN
        IF(OrigSize == 0) OrigSize = SIZE(CF)
        ALLOCATE( Matrix1 % GRows(OrigSize) )
        InvCF => Matrix1 % GRows
        InvCF = 0

        IF(Level ==  Solver % MultiGridTotal) THEN
          InvCF = CF 
        ELSE
          DO i=1,OrigSize
            j = Matrix1 % Child % Grows(i)
            IF(j > 0) InvCF(i) = CF( j )
          END DO
        END IF
      END IF            

    ELSE
      ! .NOT. new linear system
      Matrix2 => Matrix1 % Parent      
      CF => Matrix1 %  Gorder 
      IF(ALLOCATED(Matrix1 % Grows)) InvCF => Matrix1 % Grows
   END IF  
 
    n  = Matrix1 % NumberOfRows
    n2 = Matrix2 % NumberOfRows

    Residual2 => Matrix2 % RHS
    ALLOCATE( Solution2(n2) )
    
!------------------------------------------------------------------------------
!      Global iteration parameters:
!      ----------------------------

    MaxIter = 1
    IF ( Level == Solver % MultiGridTotal ) THEN
      MaxIter = ListGetInteger( Params,'MG Max Iterations', Found )
      IF ( .NOT. Found ) THEN
        IF( ListGetString( Params, 'Linear System Solver', Found ) == 'multigrid') THEN
          MaxIter = ListGetInteger( Params,'Linear System Max Iterations' )
        ELSE
          MaxIter = 1
        END IF
      END IF
      Tolerance = ListGetConstReal( Params,'MG Convergence Tolerance', Found )
      IF ( .NOT. Found ) THEN
        Tolerance = ListGetConstReal( Params,'Linear System Convergence Tolerance' )
      END IF
    ELSE
      MaxIter = ListGetInteger( Params,'MG Level Max Iterations', Found )
      IF ( .NOT. Found ) MaxIter = 1         
      Tolerance = ListGetConstReal( Params,'MG Level Convergence Tolerance', Found )
      IF ( .NOT. Found ) Tolerance = HUGE(Tolerance)
    END IF
   
!   Smoothing preconditiong, if not given diagonal preconditioning is used:
!   ----------------------------------------------------------------------
    str = ListGetString( Params, 'MG Preconditioning', Found )
    IF ( .NOT. Found ) THEN
      str = ListGetString( Params,'Linear System Preconditioning', Found )
    END IF
   
    IF ( str == 'ilut' )  THEN
      IF ( NewLinearSystem ) THEN
        ILUTOL = ListGetConstReal( Params,'MG ILUT Tolerance', GotIt )
        IF ( .NOT. GotIt ) THEN
          ILUTOL = ListGetConstReal( Params,'Linear System ILUT Tolerance' )
        END IF
        
        IF ( Parallel ) THEN
          Condition = CRS_ILUT( PMatrix, ILUTOL )
        ELSE
          Condition = CRS_ILUT( Matrix1, ILUTOL )
        END IF
      END IF
      
    ELSE IF ( SEQL(str, 'ilu') ) THEN      
      IF ( NewLinearSystem ) THEN
        k = ICHAR(str(4:4)) - ICHAR('0')
        IF ( k < 0 .OR. k > 9 ) k = 0
        IF ( Parallel ) THEN
          PMatrix % Cholesky = ListGetLogical( Params, &
                 'Linear System Symmetric ILU', Found )
          Condition = CRS_IncompleteLU( PMatrix, k )
        ELSE
          Matrix1 % Cholesky = ListGetLogical( Params, &
                 'Linear System Symmetric ILU', Found )
          Condition = CRS_IncompleteLU( Matrix1, k )
        END IF
      END IF      
    END IF


!------------------------------------------------------------------------------
!      Ok, lets go:
!      ------------
    DO iter = 1,MaxIter
      ResidualNorm = CMGSweep()
      
      WRITE(Message,'(A,I0,A,I0,A,2E20.12E3)') 'MG Residual at level: ', &
          Level, ' iter: ', iter,' is:', ResidualNorm/RHSNorm, ResidualNorm
      CALL Info( 'CMGSolve', Message, Level=5 )


      IF( ResidualNorm /= ResidualNorm .OR. ResidualNorm > 1.0d50 ) THEN
         CALL Fatal('CMGSolve','We seem to have diverged')
      END IF
      
      IF( Level == Solver % MultiGridTotal ) THEN
        IF ( ResidualNorm/RHSNorm < Tolerance ) EXIT
      ELSE
        IF ( ResidualNorm < Tolerance ) EXIT
      END IF
    END DO

  
!------------------------------------------------------------------------------
!
!      Finalize:
!      ---------
    IF ( Parallel ) THEN 
      CALL ParallelUpdateResult( Matrix1, Solution, Residual )
    END IF

    IF( Normalize ) THEN
      Solution(1:n) = Solution(1:n) * RHSNorm
      ForceVector(1:n) = ForceVector(1:n) * RHSNorm
    END IF


    DEALLOCATE( Residual, Solution2 )
    
    IF ( Level == Solver % MultiGridTotal ) THEN
      WRITE( Message, '(A,F8.2)' ) 'MG iter time: ', CPUTime() - tt
      CALL Info( 'CMGSolve', Message, Level=5 )

      IF(CoarseSave) CALL SaveClusters()
    END IF

    WRITE(Message,'(A,I2)') 'Finishing level ',Level
    CALL Info('CMGSolve',Message,Level=10)

    RETURN
!------------------------------------------------------------------------------

  CONTAINS

  
!------------------------------------------------------------------------------
    RECURSIVE FUNCTION CMGSweep() RESULT(RNorm)
!------------------------------------------------------------------------------
      INTEGER :: i,j,Rounds
      LOGICAL :: GotIt
      REAL(KIND=dp) :: RNorm
!------------------------------------------------------------------------------
      INTEGER :: Sweeps
      INTEGER, POINTER :: Iters(:)
      REAL(KIND=dp), POINTER :: R1(:),R2(:)
      REAL(KIND=dp), ALLOCATABLE :: Work2(:)

!------------------------------------------------------------------------------

!      Presmoothing:
!      -------------

      
      Iters => ListGetIntegerArray( Params,'MG Sweeps',GotIt)
      IF(GotIt) THEN
        Sweeps = Iters(MIN(InvLevel,SIZE(Iters)))
      ELSE        
        Sweeps = 1
      END IF
      
      PSolver => Solver
      CALL Info('CMGSweep','Calling presmoother',Level=9)
      RNorm = MGSmooth( PSolver, Matrix1, Solver % Mesh, Solution, ForceVector, &
          Residual, Level, DOFs, PreSmooth = .TRUE., CF = CF)
      
!------------------------------------------------------------------------------
!
!      Solve (PAQ)z = Pr, x = x + Qz:
!      ==============================
!
!      Project current residual to the lower level mesh:
!      -------------------------------------------------
      R1 => Residual(1:n)
      R2 => Residual2(1:n2)

      CALL CRS_ClusterProject( CF, R1, R2, DOFs, Trans = .FALSE. )

!      Recursively solve (PAQ)z = Pr:
!      ------------------------------

!      numbers of W-cycles


! I wonder how this really should be for multiple sweeps, jpr?

      Solution2 = 0.0_dp
      IF( Sweeps > 1 ) THEN
        ALLOCATE( Work2 ( n2 ) )
        DO i=1,Sweeps
          Work2 = Solution2        
          CALL MultigridSolve( Matrix2, Solution2, Residual2, DOFs, &
              Solver, Level-1, NewLinearSystem )          
          Solution2 = Solution2 + Work2
        END DO
        DEALLOCATE( Work2 )
      ELSE
        CALL MultigridSolve( Matrix2, Solution2, Residual2, DOFs, &
            Solver, Level-1, NewLinearSystem )                  
      END IF

!      Compute x = x + Qz:
!      -------------------
      R1 => Residual (1:n)
      R2 => Solution2(1:n2)
      
      CALL CRS_ClusterProject( CF, R2, R1, DOFs, Trans = .TRUE. )

      Solution(1:n) = Solution(1:n) + Residual(1:n)

!      Post smoothing:
!      ---------------

      CALL Info('CMGSweep','Calling postmoother',Level=9)
      RNorm = MGSmooth( PSolver, Matrix1, Solver % Mesh, Solution, ForceVector, &
          Residual, Level, DOFs, CF = CF )
!------------------------------------------------------------------------------
    END FUNCTION CMGSweep
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    FUNCTION MGnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       REAL(KIND=dp) :: s
       REAL(KIND=dp) CONTIG :: x(:)
!------------------------------------------------------------------------------
       IF ( .NOT. Parallel ) THEN
          s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
       ELSE
          s = ParallelNorm( n, x )
       END IF
!------------------------------------------------------------------------------
    END FUNCTION MGnorm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE MGmv( A, x, b, Update )
!------------------------------------------------------------------------------
       REAL(KIND=dp) CONTIG :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       IF ( .NOT. Parallel ) THEN
         CALL CRS_MatrixVectorMultiply( A, x, b )
       ELSE
         IF ( PRESENT( Update ) ) THEN
           CALL ParallelMatrixVector( A,x,b,Update,ZeroNotOwned=.TRUE. )
         ELSE
           CALL ParallelMatrixVector( A,x,b,ZeroNotOwned=.TRUE. )
         END IF
       END IF
!------------------------------------------------------------------------------
    END SUBROUTINE MGmv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION MGCnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       COMPLEX(KIND=dp) :: s
       COMPLEX(KIND=dp) CONTIG :: x(:)
!------------------------------------------------------------------------------
       s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
!------------------------------------------------------------------------------
    END FUNCTION MGCnorm
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE MGCmv( A, x, b, Update )
!------------------------------------------------------------------------------
       COMPLEX(KIND=dp) CONTIG :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       CALL CRS_ComplexMatrixVectorMultiply( A, x, b )
!------------------------------------------------------------------------------
     END SUBROUTINE MGCmv
!------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
!>     Project the residuals from coarse to fine and restrict from fine to coarse
!>     in the clustering multigrid method.
!-------------------------------------------------------------------------------
  SUBROUTINE CRS_ClusterProject( CF, u, v, DOFs, Trans )
!------------------------------------------------------------------------------
    INTEGER, POINTER :: CF(:)
    REAL(KIND=dp), POINTER :: u(:),v(:)
    INTEGER :: DOFs
    LOGICAL, OPTIONAL :: Trans
!-------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,nv,nu
    REAL(KIND=dp), POINTER CONTIG :: Values(:)
    LOGICAL :: LTrans
    INTEGER, POINTER CONTIG :: Rows(:), Cols(:)
!-------------------------------------------------------------------------------
    LTrans = .FALSE.
    IF ( PRESENT( Trans ) ) LTrans = Trans


    IF(DOFs == 1) THEN
      nu = SIZE(u) 
      nv = SIZE(v) 
      IF(Trans) THEN
        ! Only one value for each v is needed
        DO i=1,nv
          j = CF(i)
          IF(j > 0) v(i) = u(j)
        END DO
      ELSE
        v = 0.0d0
        DO i=1,nu
          j = CF(i)
          IF(j > 0) v(j) = v(j) + u(i)
        END DO
      END IF
    ELSE 
      nu = SIZE(u) / DOFs
      nv = SIZE(v) / DOFs
      IF(Trans) THEN
        DO i=1,nv
          j = CF(i)
          IF(j > 0) THEN
            ! Having the dof-loop inside the condition seems to slightly faster
            DO k=1,DOFs
              v(DOFs*(i-1)+k) = u(DOFs*(j-1)+k)
            END DO
          END IF
        END DO
      ELSE
        v = 0.0d0
        DO i=1,nu
          j = CF(i)
          IF(j > 0) THEN
            DO k=1,DOFs
              v(DOFs*(j-1)+k) = v(DOFs*(j-1)+k) + u(DOFs*(i-1)+k)
            END DO
          END IF
        END DO
      END IF
    END IF

!-------------------------------------------------------------------------------
  END SUBROUTINE CRS_ClusterProject
!-------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>     Project cluster matrix A to B: B = PAR
!>     The projector is formed implicitely.
!------------------------------------------------------------------------------
   SUBROUTINE CRS_ClusterMatrixCreate( A, CF, B, Components ) 
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A,B
      INTEGER, POINTER :: CF(:)
      INTEGER :: Components
!------------------------------------------------------------------------------
      INTEGER, POINTER :: L1(:), L2(:)
      REAL(KIND=dp) :: s, Epsilon, rowsum
      REAL(KIND=dp), POINTER :: R1(:),R2(:)
      INTEGER :: i,ia,j,k,l,m,n,NA,NB,ci,cj,ck,cl,NoRow,i2,j2,indi,prevnewrow,&
          newrow,newcol
      INTEGER :: nodecj,comp,NmatA,NmatB,mati,NnodA,NnodB,nodi,TotalNonzeros
      LOGICAL :: AllocationsDone, GotIt, debug, Normalize
      INTEGER, ALLOCATABLE :: Ind(:), Row(:), Sizes(:), CumSizes(:), ClusterOrder(:), &
          Interval(:,:)
!------------------------------------------------------------------------------

      AllocationsDone = .FALSE.

      NmatA = A % NumberOfRows
      NnodA = NmatA / Components

      NnodB = MAXVAL(CF)
      NmatB = NnodB * Components

      ALLOCATE( Sizes(NnodB), CumSizes(NnodB+1), ClusterOrder(NnodA), &
          Row( NmatA ), Ind( NmatB ) )

      Epsilon = ListGetConstReal(Params,'MG Matrix Create Epsilon',GotIt)
      IF(.NOT. GotIt) Epsilon = SQRT(TINY(Epsilon))

 
      ! Order the dofs in the order of the clusters
      Sizes = 0
      DO i=1,NnodA
        j = CF(i) 
        IF( j > 0 ) THEN
          Sizes( j ) = Sizes( j ) + 1
        END IF
      END DO

      CumSizes = 0
      DO i=1,NnodB
        CumSizes(i+1) = CumSizes(i) + Sizes(i) 
      END DO

      Sizes = 0
      ClusterOrder = 0

      DO i=1,NnodA
        j = CF(i)
        IF(j == 0) CYCLE
        Sizes( j ) = Sizes( j ) + 1
        indi = Sizes( j ) + CumSizes( j )
        ClusterOrder(indi) = i
      END DO

      B => AllocateMatrix()      
      B % NumberOfRows = NmatB 
      ALLOCATE( B % Rows( NmatB + 1 ), &
          B % Diag( NmatB ), &
          B % RHS( NmatB ) )        
      B % RHS = 0.0d0
      B % Diag = 0
      B % Rows(1) = 1


      IF(.FALSE.) THEN
        PRINT *,'Initial Matrix'
        DO i=1,NnodA
          PRINT *,'i',i
          DO j=A % Rows(i),A % Rows(i+1)-1          
            PRINT *,'j',A % Cols(j),A % Values(j)
          END DO
        END DO
      END IF

10    TotalNonzeros = 0
      NoRow = 0
      Row = 0
      Ind = 0
      prevnewrow = 0

      IF( Components == 1) THEN
        NA = NnodA
        NB = NnodB

        DO indi=1,NA

          i = ClusterOrder(indi) 
          IF(i == 0) CYCLE

          newrow = CF(i)  
IF(newrow < prevnewrow ) PRINT *,'problem:',indi,i,newrow,prevnewrow        
          IF(prevnewrow /= newrow) THEN
            DO j=1,NoRow
              Row(Ind(j)) = 0
            END DO
            NoRow = 0
            prevnewrow = newrow
          END IF

          DO j=A % Rows(i),A % Rows(i+1)-1          
            cj = A % Cols(j)

            newcol = CF(cj)
            IF(newcol == 0) CYCLE
 
            IF( .FALSE. .AND. ABS(A % Values(j)) < Epsilon) CYCLE
 
            i2 = Row(newcol)

            IF ( i2 == 0) THEN
              NoRow = NoRow + 1
              TotalNonzeros = TotalNonzeros + 1
              Ind(NoRow) = newcol
              i2 = B % Rows(newrow) + NoRow - 1                   
              Row(newcol) = i2

              IF(AllocationsDone) THEN
                IF(newrow == newcol) B % Diag(newrow) = i2
                B % Cols(i2) = newcol                     
                B % Values(i2) = A % Values(j)
              END IF
            ELSE IF(AllocationsDone) THEN
              B % Values(i2) = B % Values(i2) + A % Values(j)
            END IF
          END DO
          
          B % Rows(newrow+1) = B % Rows(newrow) + NoRow
         
        END DO
      ELSE
 
        DO indi=1,NnodB          

          DO comp=1,Components
            
            DO nodi=CumSizes(indi)+1,CumSizes(indi+1)
              
              i = ClusterOrder(nodi) 
              IF(i == 0) CYCLE          
              mati = Components*(i-1) + comp
              newrow = Components*(CF(i)-1) + comp
              
              IF(prevnewrow /= newrow) THEN
                DO j=1,NoRow
                  Row(Ind(j)) = 0
                END DO
                NoRow = 0
                prevnewrow = newrow
              END IF
              
              DO j=A % Rows(mati),A % Rows(mati+1)-1          
                cj = A % Cols(j)
                nodecj = (cj-1)/Components + 1
                k = CF(nodecj)
                IF(k == 0) CYCLE
                
                IF( ABS(A % Values(j)) < Epsilon) CYCLE
                
                newcol = Components*(k-1) + MOD(cj-1,Components) + 1
                
                i2 = Row(newcol)
                
                IF ( i2 == 0) THEN
                  NoRow = NoRow + 1
                  TotalNonzeros = TotalNonzeros + 1
                  Ind(NoRow) = newcol
                  i2 = B % Rows(newrow) + NoRow - 1                   
                  Row(newcol) = i2
                  
                  IF(AllocationsDone) THEN
                    IF(newrow == newcol) B % Diag(newrow) = i2
                    B % Cols(i2) = newcol                     
                    B % Values(i2) = A % Values(j)
                  END IF
                ELSE IF(AllocationsDone) THEN
                  B % Values(i2) = B % Values(i2) + A % Values(j)
                END IF
              END DO
              B % Rows(newrow+1) = B % Rows(newrow) + NoRow         
              
            END DO
          END DO
        END DO
      END IF

      IF(.NOT. AllocationsDone) THEN
        ALLOCATE( B % Cols( TotalNonzeros ), B % Values( TotalNonzeros ) )
        B % Cols = 0
        B % Values = 0.0d0
        AllocationsDone = .TRUE.
        GOTO 10 
      END IF

      IF(.FALSE.) THEN
        PRINT *,'Created Matrix'
        DO i=1,NB
          PRINT *,'i',i
          DO j=B % Rows(i),B % Rows(i+1)-1          
            PRINT *,'j',B % Cols(j),B % Values(j)
          END DO
        END DO
      END IF

      Normalize = ListGetLogical( Params,'MG Matrix Enforce Zero Rowsum',GotIt )
      IF( Normalize ) THEN
        CALL Info('CRS_ClusterMatrixCreate','Enforcing rowsum to zero')

        DO j = 1, B % NumberOfRows
          rowsum = 0.0_dp
          DO i = B % Rows(j),B % Rows(j+1)-1
            rowsum = rowsum + B % Values( i ) 
          END DO

          k = B % Diag(j)
          IF( ABS( rowsum - B % Values(k) ) < 1.0e-3 * ABS( rowsum ) ) THEN
            PRINT *,'Dirichlet node:',j
            CYCLE
          END IF
          B % Values(k) = B % Values(k) - rowsum
        END DO
      END IF


      DEALLOCATE(  Sizes, CumSizes, ClusterOrder, Row, Ind )
      
      WRITE(Message,'(A,F10.3)') 'Coarse matrix reduction factor',&
          1.0 *  SIZE(A % Cols) / TotalNonZeros
      CALL Info('CRS_ClusterMatrixCreate',Message)

    END SUBROUTINE CRS_ClusterMatrixCreate


!------------------------------------------------------------------------------
    SUBROUTINE SaveClusters() 
!------------------------------------------------------------------------------

      LOGICAL :: Visited = .FALSE., GotIt
      TYPE(Matrix_t), POINTER :: TmpMatrix
      REAL(KIND=dp), POINTER :: Clustering(:)
      INTEGER, POINTER :: Perm(:)
      INTEGER :: i,j,k,m,istat

      SAVE Visited

      IF(Visited) RETURN
      Visited = .TRUE.

      PRINT *,'Saving clusters in levels:',Level,Solver % MultiGridTotal
      
      ! Save in ElmerPost
      IF(.TRUE.) THEN
        m = ListGetInteger(Params,'MG Cluster Save Modulo',GotIt)

        Perm => Solver % Variable % Perm
        TmpMatrix => Matrix1

        DO j=Solver % MultiGridTotal, MinLevel,-1
          k = Solver % MultiGridTotal - j + 1

	  IF( .NOT. ALLOCATED(TmpMatrix % Grows)) CYCLE

          NULLIFY( Clustering ) 
          ALLOCATE( Clustering(OrigSize), STAT=istat)
          IF ( istat /= 0 ) CALL Fatal( 'SaveClusters', 'Memory allocation error.' )           
          IF(m == 0) THEN
            Clustering = 1.0d0 * TmpMatrix % Grows
          ELSE
            Clustering = 1.0d0 * MODULO(TmpMatrix % Grows,m)
          END IF
          IF( k < 10) THEN
            WRITE(Message,'(A,I1)') 'Cluster',k
          ELSE
            WRITE(Message,'(A,I2)') 'Cluster',k
          END IF

          CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
              Solver,TRIM(Message),1,Clustering, Perm)
          TmpMatrix => TmpMatrix % Parent
        END DO        

      END IF

      ! Save in simple dat file (Matlab) 
      IF(.FALSE.) THEN
        OPEN (10,FILE='clusters.dat')        
        DO i=1,OrigSize
          WRITE (10,'(3ES17.8E3)',ADVANCE='NO') &
              Mesh % Nodes % X(i), Mesh % Nodes % Y(i), Mesh % Nodes % Z(i)
          TmpMatrix => Matrix1
          DO j=Solver % MultiGridTotal, MinLevel+1,-1
            TmpMatrix => TmpMatrix % Parent
            WRITE (10,'(I6)',ADVANCE='NO') TmpMatrix % Grows(i)
          END DO
          WRITE (10,'(A)') ' ' 
        END DO
        CLOSE(10)
      END IF

    END SUBROUTINE SaveClusters
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  END SUBROUTINE CMGSolve
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE MSolverActivate( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
     USE MeshUtils
     TYPE(Model_t)  :: Model
     TYPE(Solver_t),TARGET :: Solver
     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: dt, OrigDT, DTScal
!------------------------------------------------------------------------------
     LOGICAL :: stat, Found, GB, FirstTime=.TRUE., MeActive, NamespaceFound
     INTEGER :: i, j, n, SolverAddr, BDOFs, execi, timestep, maxdim
     REAL(KIND=dp) :: st
     TYPE(Variable_t), POINTER :: TimeVar, IterV
     TYPE(Element_t), POINTER :: CurrentElement
     CHARACTER(LEN=MAX_NAME_LEN) :: EquationName, str

     INTEGER :: comm_active, group_active, group_world, ierr
     INTEGER, ALLOCATABLE :: memb(:)
     TYPE(ValueList_t), POINTER :: Params

     SAVE TimeVar, FirstTime
!------------------------------------------------------------------------------
     CALL SetCurrentMesh( Model, Solver % Mesh )
     Model % Solver => Solver
     Params => Solver % Values

     st = ListGetConstReal( Params, 'Start Time', Found )
     IF ( Found ) THEN
       TimeVar => VariableGet( Model % Variables, 'Time' )
       IF ( TimeVar % Values(1) < st ) RETURN
     END IF

     execi = ListGetInteger( Params, 'Exec Interval', Found )
     IF ( Found ) THEN
       TimeVar => VariableGet( Model % Variables, 'Timestep' )
       execi = MOD( NINT(Timevar % Values(1))-1, execi )
       IF ( execi /= 0 ) RETURN
     END IF

!    IF ( Solver % Mesh % Changed .OR. Solver % NumberOfActiveElements <= 0 ) THEN
       Solver % NumberOFActiveElements = 0
       EquationName = ListGetString( Params, 'Equation', stat )

       IF ( Stat ) THEN
          IF ( ASSOCIATED(Solver % ActiveElements) ) DEALLOCATE( Solver % ActiveElements )
          ALLOCATE( Solver % ActiveElements( Solver % Mesh % NumberOfBulkElements + &
                       Solver % Mesh % NumberOFBoundaryElements ) )

          Maxdim = 0
          DO i=1,Solver % Mesh % NumberOfBulkElements+Solver % Mesh % NumberOFBoundaryElements
             CurrentElement => Solver % Mesh % Elements(i)
             IF ( CheckElementEquation( Model, CurrentElement, EquationName ) ) THEN
                Solver % NumberOfActiveElements = Solver % NumberOFActiveElements + 1
                Solver % ActiveElements( Solver % NumberOFActiveElements ) = i
                Maxdim = MAX( CurrentElement % TYPE % DIMENSION, Maxdim )
             END IF
          END DO
          CALL ListAddInteger( Params, 'Active Mesh Dimension', Maxdim )
       END IF
!    END IF
!------------------------------------------------------------------------------
     Solver % Mesh % OutputActive = .TRUE.
     OrigDT = dt
     DTScal = ListGetConstReal( Params,'Timestep Scale', Found )
     IF ( .NOT. Found ) DTScal = 1.0d0
     Solver % dt = DTScal * dt

     MeActive = ASSOCIATED(Solver % Matrix)
     IF ( MeActive ) &
        MeActive = MeActive .AND. (Solver % Matrix % NumberOfRows>0)
     CALL ParallelActive( MeActive )

     IF ( ASSOCIATED(Solver % Matrix) ) Solver % Matrix % Comm = ELMER_COMM_WORLD

     IF ( ParEnv % PEs>1 ) THEN
       DO i=1,ParEnv % PEs
         IF ( ParEnv % Active(i) ) THEN
           EXIT
         END IF
       END DO

       OutputPE = -1
       IF ( i-1==ParEnv % MyPE ) OutputPE=0

       n = COUNT(ParEnv % Active)
       IF ( n>0 .AND. n<ParEnv % PEs ) THEN
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
         Solver % Matrix % Comm = comm_active
       END IF
     END IF

     str = ListGetString( Params, 'Namespace', NamespaceFound )
     IF (NamespaceFound) CALL ListPushNamespace(TRIM(str))

     iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
     iterV % Values(1) = 1

     str = ListGetString( Params, 'Procedure', Found )

#ifdef SGIn32
      SolverAddr = Solver % PROCEDURE
      CALL ExecSolver( SolverAddr, Model, Solver, DTScal*dt, TransientSimulation)
#else
      CALL ExecSolver( &
             Solver % PROCEDURE, Model, Solver, DTScal*dt, TransientSimulation)
#endif

     IF(NamespaceFound) CALL ListPopNamespace()
     IF ( ASSOCIATED(Solver % Matrix) ) THEN
       IF ( Solver % Matrix % Comm /= ELMER_COMM_WORLD ) &
          CALL MPI_Comm_Free( Solver % Matrix % Comm, ierr )
     END IF

     Solver % dt = OrigDT
!------------------------------------------------------------------------------
   END SUBROUTINE MSolverActivate
!------------------------------------------------------------------------------

END MODULE Multigrid

!> \}
