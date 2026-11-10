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
! Comput the product of a covariance matrix with a vector
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


      TYPE(CovarianceState_t), POINTER :: State
      REAL(kind=dp) :: sigma2
      INTEGER :: Op

      LOGICAL :: Normalize
      LOGICAL :: Parallel
      LOGICAL :: Found

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="Cov.x Solver"

      ! check Parallel/Serial
      Parallel=(ParEnv % PEs > 1)

      SolverParams => GetSolverParams()
      State => GetCovarianceState(Solver)

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
      IF (.NOT. State % Initialized) THEN

        CALL GetActiveNodesSet(Solver,State % nn,State % ActiveNodes,State % InvPerm,State % PbDim)

        !Sanity check
        IF (ANY(Perm(State % ActiveNodes(1:State % nn)).LT.0)) &
          CALL FATAL(SolverName,"Pb with input variable perm")

        !! The covariance type
        State % CovType = ListGetString(SolverParams,"Covariance type",UnFoundFatal=.TRUE.)
        State % std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.)

        SELECT CASE (State % CovType)

          CASE('diagonal')
            CALL INFO(SolverName,"Using diagonal covariance",level=3)

          CASE('full matrix')
            CALL INFO(SolverName,"Using full matrix covariance",level=3)

            Op=1
            ALLOCATE(State % aap(State % nn*(State % nn+1)/2))
            CALL CovarianceInit(Solver,State % nn,State % InvPerm,State % aap,Op,State % PbDim)

          CASE('diffusion operator')
            CALL INFO(SolverName,"Using diffusion operator covariance",level=3)

            CALL CovarianceInit(Solver,State % MSolver,State % KMSolver)

        END SELECT

       allocate(State % x(State % nn),State % y(State % nn))

       IF (Normalize) THEN
          allocate(State % norm(State % nn))

          !input vector
          State % x(:) = 1._dp

          ! C . x
          SELECT CASE (State % CovType)

            CASE('diagonal')
              sigma2=State % std**2
              State % norm(:)=sigma2*State % x(:)

            CASE('full matrix')
              CALL CovarianceVectorMultiply(Solver,State % nn,State % aap,State % x,State % norm)

            CASE('diffusion operator')
              ! y = SIGMA C SIGMA . x
              CALL CovarianceVectorMultiply(Solver,State % MSolver,State % KMSolver,State % nn,State % x,State % norm)

           END SELECT
       END IF

       State % Initialized = .TRUE.
      END IF

      !input vector
      State % x(Solver%Variable%Perm(State % ActiveNodes(1:State % nn))) = Values(Perm(State % ActiveNodes(1:State % nn)))

      ! C . x
      SELECT CASE (State % CovType)

        CASE('diagonal')
          sigma2=State % std**2
          State % y(:)=sigma2*State % x(:)

        CASE('full matrix')
          CALL CovarianceVectorMultiply(Solver,State % nn,State % aap,State % x,State % y)

        CASE('diffusion operator')
          ! y = SIGMA C SIGMA . x
          CALL CovarianceVectorMultiply(Solver,State % MSolver,State % KMSolver,State % nn,State % x,State % y)

      END SELECT

      IF (Normalize) State % y(:)=State % y(:)/State % norm(:)

      Solver % Variable % Values(Solver%Variable%Perm(State % ActiveNodes(1:State % nn)))=&
                    State % y(Solver%Variable%Perm(State % ActiveNodes(1:State % nn)))

     END SUBROUTINE CovarianceVectorMultiplySolver
