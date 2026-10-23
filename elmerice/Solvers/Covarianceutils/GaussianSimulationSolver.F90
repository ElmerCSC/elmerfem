!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! ******************************************************************************
! *
! *  Authors: F. Gillet-Chaulet
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 24/06/2024
! *
! *****************************************************************************
!***********************************************************************************************
! Generate random realization from a given covariance matrix
!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************
      SUBROUTINE GaussianSimulationSolver( Model,Solver,dt,TransientSimulation )
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

      TYPE(Variable_t), POINTER :: Var,Var_b
      REAL(KIND=dp), POINTER ::  Values(:),Values_b(:)
      INTEGER, POINTER :: Perm(:),Perm_b(:)
      CHARACTER(LEN=MAX_NAME_LEN) :: Varbname
      INTEGER :: DOFs

      INTEGER :: i,k

      INTEGER :: Op

      TYPE(CovarianceState_t), POINTER :: State
      Logical :: Parallel
      LOGICAL :: Found

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="GaussianSimualtion"

      Integer,allocatable :: seed(:)
      Integer :: ssize

      ! check Parallel/Serial
      Parallel=(ParEnv % PEs > 1)

      SolverParams => GetSolverParams()
      State => GetCovarianceState(Solver)

      Var => Solver % Variable
      IF (.NOT.ASSOCIATED(Var)) &
        CALL FATAL(SolverName,'Variable not associated')
      Values => Var % Values
      Perm => Var % Perm
      DOFs = Var % DOFs

      VarbName = ListGetString(SolverParams,"Background Variable name",UnFoundFatal=.TRUE.)

      Var_b => VariableGet(Solver % Mesh % Variables,TRIM(VarbName),UnFoundFatal=.TRUE.)
      Values_b => Var_b % Values
      Perm_b => Var_b % Perm
      IF (Var_b % DOFs.GT.1) &
         CALL FATAL(SolverName,'DoFs for mean variable should be 1')

      !! some initialisation
      IF (.NOT. State % Initialized) THEN
        CALL GetActiveNodesSet(Solver,State % nn,State % ActiveNodes,State % InvPerm,State % PbDim)

        !Sanity check
        IF (ANY(Perm_b(State % ActiveNodes(1:State % nn)).LT.0)) &
          CALL FATAL(SolverName,"Pb with background variable perm")

        !! The covariance type
        State % CovType = ListGetString(SolverParams,"Covariance type",UnFoundFatal=.TRUE.)
        State % std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.)

        SELECT CASE (State % CovType)

          CASE('diagonal')
            CALL INFO(SolverName,"Using diagonal covariance",level=3)

          CASE('full matrix')
            CALL INFO(SolverName,"Using full matrix covariance",level=3)

            Op=2
            ALLOCATE(State % aap(State % nn*(State % nn+1)/2))
            CALL CovarianceInit(Solver,State % nn,State % InvPerm,State % aap,Op,State % PbDim)

          CASE('diffusion operator')
            CALL INFO(SolverName,"Using diffusion operator covariance",level=3)

            CALL CovarianceInit(Solver,State % MSolver,State % KMSolver)

        END SELECT

       allocate(State % x(State % nn),State % y(State % nn),State % rr(State % nn,DOFs))

       State % Initialized = .TRUE.
      END IF

     !  CALL random_seed()
      CALL random_seed(size=ssize)
      allocate(seed(ssize))
      seed = ListGetInteger( SolverParams , 'Random Seed',Found )
      IF (Found)  call random_seed( put=seed )
      CALL random_seed(get=seed)
      deallocate(seed)

       !Create DOFs random vectors of size n
       State % rr=0._dp
       DO k=1,DOFs
         DO i=1,State % nn
          State % rr(i,k)=NormalRandom()
         END DO
       END DO

       DO k=1,DOFs
         State % x(Perm(State % ActiveNodes(1:State % nn)))=State % rr(State % ActiveNodes(1:State % nn),k)

        SELECT CASE (State % CovType)
          CASE('diagonal')
              State % y(:) = State % std*State % x(:)

          CASE('full matrix')
              CALL SqrCovarianceVectorMultiply(Solver,State % nn,State % aap,State % x,State % y)

          CASE('diffusion operator')
             CALL SqrCovarianceVectorMultiply(Solver,State % MSolver,State % KMSolver,State % nn,State % x,State % y)

        END SELECT

        Values(DOFs*(Perm(State % ActiveNodes(1:State % nn))-1)+k)=Values_b(Perm_b(State % ActiveNodes(1:State % nn)))+&
                State % y(Perm(State % ActiveNodes(1:State % nn)))
       END DO


     END SUBROUTINE GaussianSimulationSolver
