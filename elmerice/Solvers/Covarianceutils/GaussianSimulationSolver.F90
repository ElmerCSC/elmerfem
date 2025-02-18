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

      TYPE(Solver_t), POINTER, SAVE :: MSolver,KMSolver
      REAL(kind=dp),allocatable,save :: aap(:) ! matrix in packed format
      REAL(kind=dp),allocatable,save :: x(:),y(:)
      REAL(kind=dp),allocatable,SAVE :: rr(:,:)

      INTEGER,SAVE  :: nn
      INTEGER, ALLOCATABLE, SAVE :: ActiveNodes(:),InvPerm(:)
      INTEGER,SAVE :: PbDim

      CHARACTER(LEN=MAX_NAME_LEN),SAVE :: CovType
      REAL(kind=dp),SAVE :: std
      INTEGER :: Op

      LOGICAL, SAVE :: Firsttime=.TRUE.
      Logical :: Parallel
      LOGICAL :: Found

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="GaussianSimualtion"

      Integer,allocatable :: seed(:)
      Integer :: ssize

      ! check Parallel/Serial
      Parallel=(ParEnv % PEs > 1)

      SolverParams => GetSolverParams()

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
      IF (Firsttime) THEN
        CALL GetActiveNodesSet(Solver,nn,ActiveNodes,InvPerm,PbDim)

        !Sanity check
        IF (ANY(Perm_b(ActiveNodes(1:nn)).LT.0)) &
          CALL FATAL(SolverName,"Pb with background variable perm")

        !! The covariance type
        CovType = ListGetString(SolverParams,"Covariance type",UnFoundFatal=.TRUE.)
        std = ListGetConstReal(SolverParams,"standard deviation",UnFoundFatal=.TRUE.)

        SELECT CASE (CovType)

          CASE('diagonal')
            CALL INFO(SolverName,"Using diagonal covariance",level=3)

          CASE('full matrix')
            CALL INFO(SolverName,"Using full matrix covariance",level=3)

            Op=2
            ALLOCATE(aap(nn*(nn+1)/2))
            CALL CovarianceInit(Solver,nn,InvPerm,aap,Op,PbDim)

          CASE('diffusion operator')
            CALL INFO(SolverName,"Using diffusion operator covariance",level=3)

            CALL CovarianceInit(Solver,MSolver,KMSolver)

        END SELECT

       allocate(x(nn),y(nn),rr(nn,DOFs))

       Firsttime=.FALSE.
      END IF

     !  CALL random_seed()
      CALL random_seed(size=ssize)
      allocate(seed(ssize))
      seed = ListGetInteger( SolverParams , 'Random Seed',Found )
      IF (Found)  call random_seed( put=seed )
      CALL random_seed(get=seed)
      deallocate(seed)

       !Create DOFs random vectors of size n
       rr=0._dp
       DO k=1,DOFs
         DO i=1,nn
          rr(i,k)=NormalRandom()
         END DO
       END DO

       DO k=1,DOFs
         x(Perm(ActiveNodes(1:nn)))=rr(ActiveNodes(1:nn),k)

        SELECT CASE (CovType)
          CASE('diagonal')
              y(:) = std*x(:)

          CASE('full matrix')
              CALL SqrCovarianceVectorMultiply(Solver,nn,aap,x,y)

          CASE('diffusion operator')
             CALL SqrCovarianceVectorMultiply(Solver,MSolver,KMSolver,nn,x,y)

        END SELECT

        Values(DOFs*(Perm(ActiveNodes(1:nn))-1)+k)=Values_b(Perm_b(ActiveNodes(1:nn)))+&
                y(Perm(ActiveNodes(1:nn)))
       END DO


     END SUBROUTINE GaussianSimulationSolver
