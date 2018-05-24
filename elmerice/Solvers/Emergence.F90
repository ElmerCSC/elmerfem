!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Authors: Thomas Zwinger,
! *  Email:   Thomas.Zwinger@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 1 April 2014
! * 
! *
! ****************************************************************************/


!-----------------------------------------------------------------------------
!>  Routine for deducing the emergence velocity
!>  as a scalar product of velcoity with a given surface normal
!-----------------------------------------------------------------------------
SUBROUTINE GetEmergenceVelocity( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE Differentials
  USE MaterialModels
  IMPLICIT NONE

  !------------------------------------------------------------------------------
  !    external variables
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t):: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation 
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------

  LOGICAL :: firstTime=.TRUE., Found, AllocationsDone = .FALSE., stat
  INTEGER :: i,j,k, t,N,NMAX, DIM, NSDOFs, istat
  INTEGER, POINTER :: FlowPerm(:), NodeIndexes(:), NPerm(:), Permutation(:)
  REAL(KIND=dp), POINTER :: FlowSolution(:), Nvector(:), EmergenceVelocity(:)
  REAL(KIND=dp), ALLOCATABLE :: Velo(:,:)
  REAL(KIND=dp) :: LocalNormalVector(3)
  CHARACTER(LEN=MAX_NAME_LEN)  :: SolverName, FlowSolName,ConvectionFlag
  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: FlowSol, NormalSolution, EmergenceSol
  TYPE(ValueList_t), POINTER :: SolverParams, Material, Equation
  TYPE(Matrix_t), POINTER :: Systemmatrix
  !---------------------------------------

  SAVE SolverName, firstTime, DIM, Velo, AllocationsDone

  IF (firstTime) THEN
    WRITE(SolverName, '(A)')  'GetEmergenceVelocity'
    WRITE(Message, '(A)') "Starting " // TRIM(SolverName)
    CALL INFO(SolverName, Message)
    DIM = CoordinateSystemDimension()
    firstTime = .FALSE.
  END IF

  IF ( (.NOT.AllocationsDone) .OR. (Solver % Mesh % Changed) ) THEN
    IF (AllocationsDone) DEALLOCATE(Velo)
    NMAX = Model % MaxElementNodes
    ALLOCATE(Velo( 3, NMAX ), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal(SolverName,'Memory allocation error 1, Aborting.')
    END IF
    AllocationsDone = .TRUE.
  END IF

  EmergenceSol => Solver % Variable
  IF (.NOT.ASSOCIATED(EmergenceSol)) &
       CALL Fatal(SolverName, 'Solver variable not found')
  Permutation  => EmergenceSol % Perm
  EmergenceVelocity => EmergenceSol % Values

  ! get normal vectors
  NormalSolution => VariableGet( Solver % Mesh % Variables, 'Normal Vector' )
  IF ( ASSOCIATED( NormalSolution ) ) THEN 
    Nvector => NormalSolution % Values
    NPerm => NormalSolution % Perm
  ELSE
    CALL Fatal(SolverName, '>Normal Vector< not found')
  END IF

  DO t=1,Solver % NumberOfActiveElements
    CurrentElement => GetActiveElement(t)
    N = GetElementNOFNodes()
    NodeIndexes => CurrentElement % NodeIndexes
    Equation => GetEquation()
    Material => GetMaterial()

    ! get flow soulution and velocity field from it
    !----------------------------------------------
    ConvectionFlag = GetString( Equation, 'Convection', Found )
    IF (.NOT. Found) &
         CALL Fatal(SolverName, 'No string for keyword > Convection < found in Equation')
    Velo = 0.0_dp
    ! constant (i.e., in section Material given) velocity
    !----------------------------------------------------
    IF ( ConvectionFlag == 'constant' ) THEN
      Velo(1,1:N) = GetReal( Material, 'Convection Velocity 1', Found )
      IF ( .NOT.Found ) &
           Velo(1,1:N) = GetReal( Equation, 'Convection Velocity 1', Found )

      Velo(2,1:N) = GetReal( Material, 'Convection Velocity 2', Found )
      IF ( .NOT.Found ) &
           Velo(2,1:N) = GetReal( Equation, 'Convection Velocity 2', Found )

      Velo(3,1:N) = GetReal( Material, 'Convection Velocity 3', Found )
      IF ( .NOT.Found ) &
           Velo(3,1:N) = GetReal( Equation, 'Convection Velocity 3', Found )
      ! computed velocity
      !------------------
    ELSE IF (ConvectionFlag == 'computed' ) THEN
      FlowSolName =  GetString( Equation,'Flow Solution Name', Found)
      IF(.NOT.Found) THEN        
        CALL WARN(SolverName,'Keyword > Flow Solution Name < not found in section >Equation<')
        CALL WARN(SolverName,'Taking default value > Flow Solution <')
      END IF
      FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName )
      IF ( ASSOCIATED( FlowSol ) ) THEN
        FlowPerm     => FlowSol % Perm
        NSDOFs     =  FlowSol % DOFs
        FlowSolution => FlowSol % Values
      ELSE
        WRITE(Message,'(A,A,A)') &
             'Convection flag set to > computed <, but no variable >',FlowSolName,'< found'
        CALL FATAL(SolverName,Message)              
      END IF
      ! get velocity profile
      IF ( ASSOCIATED( FlowSol ) ) THEN
        DO i=1,n
          j = NSDOFs*FlowPerm(NodeIndexes(i))
          IF((DIM == 2) .AND. (NSDOFs == 3)) THEN
            Velo(1,i) = FlowSolution( j-2 ) 
            Velo(2,i) = FlowSolution( j-1 ) 
            Velo(3,i) = 0.0_dp
          ELSE IF ((DIM == 3) .AND. (NSDOFs == 4)) THEN
            Velo(1,i) = FlowSolution( j-3 ) 
            Velo(2,i) = FlowSolution( j-2 ) 
            Velo(3,i) = FlowSolution( j-1 ) 
          ELSE IF ((CurrentCoordinateSystem() == CylindricSymmetric) &
               .AND. (DIM == 2) .AND. (NSDOFs == 4)) THEN
            Velo(1,i) = FlowSolution( j-3 ) 
            Velo(2,i) = FlowSolution( j-2 ) 
            Velo(3,i) = FlowSolution( j-1 ) 
          ELSE
            WRITE(Message,'(a,i0,a,i0,a)')&
                 'DIM=', DIM, ' NSDOFs=', NSDOFs, ' does not combine. Aborting'
            CALL FATAL( SolverName, Message)
          END IF
        END DO
      ELSE ! "none"
        Velo=0.0_dp
        CALL WARN(SolverName,'No variable for  velocity found - reset to zero.')
      END IF
    END IF
    DO i=1,N
      k = NPerm(NodeIndexes( i ))
      LocalNormalVector(1:DIM) = Nvector(DIM*(k-1)+1:DIM*k)
      EmergenceVelocity(Permutation(NodeIndexes( i ))) = SUM(LocalNormalVector(1:DIM) * Velo(1:DIM,i))
    END DO
  END DO
END SUBROUTINE GetEmergenceVelocity
