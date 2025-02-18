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
! *  Subroutine for saving projections to new fields.
! *
! ******************************************************************************
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 1.10.2024
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{


!------------------------------------------------------------------------------
SUBROUTINE SaveProjection_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------  
  CALL ListAddNewString( Solver % Values,'Variable',&
      '-nooutput -global SaveProjection_var') 

END SUBROUTINE SaveProjection_init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Routine for saving projections from fields as fields. 
!------------------------------------------------------------------------------
SUBROUTINE SaveProjection( Model,Solver,dt,Transient )
  USE DefUtils
  USE Types
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var, TargetVar
  INTEGER :: i,j,n,NoVar
  CHARACTER(LEN=MAX_NAME_LEN) ::  VarName, TargetName
  INTEGER, POINTER :: UnitPerm(:)
  REAL(KIND=dp) :: Nrm
  LOGICAL :: Found, Normalize, ToSlave, ToMaster
  INTEGER, POINTER :: ActiveProjs(:)

  CALL Info('SaveProjection','Creating selected projected values as fields',Level=8)

  Params => GetSolverParams()

  ! Because this solver does not have a DefaultInitialize slot the nonlinear projectors
  ! are not initialized. Do it here.
  CALL GenerateProjectors(Model,Solver,Nonlinear = .TRUE. )
 
  i = 0
  DO WHILE(.TRUE.)  
    i = i + 1
    VarName = ListGetString( Params,'Variable '//I2S(i), Found )
    IF(.NOT. Found) EXIT    
    Var => VariableGet( Model % Variables, TRIM(VarName) )
    IF(.NOT. ASSOCIATED(Var)) THEN
      CALL Warn('SaveProjection','Requested variable does not exist!')
      EXIT
    END IF
  END DO
  NoVar = i-1
  CALL Info('SaveProjection','Saving projections from '//I2S(NoVar)//' fields')
  IF(NoVar==0) RETURN

  Nrm = 0.0_dp
  DO i=1,NoVar    
    VarName = ListGetString( Params,'Variable '//I2S(i), Found )
    Var => VariableGet( Model % Variables, TRIM(VarName), ThisOnly = .TRUE.)
    CALL info('SaveProjection','Doing variable: '//TRIM(VarName),Level=8)
    
    TargetName = ListGetString( Params,'Target Variable '//I2S(i), Found )
    IF(.NOT. Found) TargetName = 'Projection '//TRIM(VarName)

    Normalize = ListGetLogical( Params,'Normalize '//I2S(i), Found )
    IF(.NOT. Found) Normalize = .TRUE.
    ToSlave = ListGetLogical( Params,'Project To Slave '//I2S(i),Found ) 
    ToMaster = ListGetLogical( Params,'Project To Master '//I2S(i),Found ) 
    
    TargetVar => VariableGet( Model % Variables, TRIM(TargetName), ThisOnly = .TRUE.)
    IF(.NOT. ASSOCIATED(TargetVar)) THEN
      IF(.NOT. ASSOCIATED(Var % Perm) ) THEN
        UnitPerm => NULL()
        ALLOCATE(UnitPerm(SIZE(Var % Values)))
        DO j=1,SIZE(Var % Values)
          UnitPerm(j) = j
        END DO
        CALL VariableAddVector( Model % Mesh % Variables, Solver % Mesh, Solver, &
            TRIM(TargetName), Var % Dofs, Perm = UnitPerm ) 
      ELSE
        CALL VariableAddVector( Model % Mesh % Variables, Solver % Mesh, Solver, &
            TRIM(TargetName), Var % Dofs, Perm = Var % Perm, Secondary = .TRUE.)
      END IF
      TargetVar => VariableGet( Model % Variables, TRIM(TargetName) )       
    END IF

    ActiveProjs => ListGetIntegerArray( Params,'Active Projectors '//I2S(i),Found )
    
    ! Do additive projection!
    CALL ProjectToVariable()
    Nrm = Nrm + SUM(TargetVar % Values**2)
  END DO
  
  Nrm = SQRT(Nrm)
  IF(SIZE(Solver % Variable % Values) == 1 ) THEN
    Solver % Variable % Values = Nrm
  END IF
  
  WRITE(Message,'(A,ES12.3)') 'Combined L2 norm of all projected fields: ', Nrm
  CALL Info('SaveProjection',Message)
  
CONTAINS


  SUBROUTINE ProjectToVariable()
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: bc, dofs, i, j, k, pi, pj, m
    INTEGER, POINTER :: Rows(:), Cols(:)
    LOGICAL :: acti, actj, AddThis
    REAL(KIND=dp) :: r1
    REAL(KIND=dp), POINTER :: Values(:)      
    REAL(KIND=dp), ALLOCATABLE :: Weight(:)
    LOGICAL, POINTER :: IsInvInvPerm(:)
    
    dofs = Var % Dofs
    TargetVar % Values = 0.0_dp
    
    IF(Normalize) THEN
      ALLOCATE(Weight(SIZE(TargetVar % Values)))
      Weight = 0.0_dp      
    END IF

    ! Go through all the projectors.
    ! There could be perhaps reason to skip some, but this will do for now.
    DO bc=1,Model % NumberOfBCs        
      IF(ASSOCIATED(ActiveProjs)) THEN
        IF(.NOT. ANY(ActiveProjs == bc)) CYCLE 
      END IF

      A => CurrentModel % BCs(bc) % PMatrix
      IF(.NOT. ASSOCIATED(A) ) THEN
        A => Solver % MortarBCs(bc) % Projector
      END IF
           
      IF(.NOT. ASSOCIATED(A)) CYCLE
      n = A % NumberOfRows
      IF(n==0) CYCLE

      CALL Info('SaveProjection','Doing projection for BC '//I2S(bc)//' of size '//I2S(n),Level=20)
      
      Rows => A % Rows
      Cols => A % Cols
      Values => A % Values

      IF(.NOT. ASSOCIATED(A % InvPerm)) THEN
        CALL Fatal('SaveProjection','InvPerm not associated!')
      END IF

      ! Create table telling which is slave/master dof. 
      m = MAXVAL(Cols)
      ALLOCATE(IsInvInvPerm(m))
      IsInvInvPerm = .FALSE.
      DO i=1,n
        pi = A % InvPerm(i)
        IF(pi > 0 .AND. pi <= m) IsInvInvPerm(pi) = .TRUE.
      END DO


      DO k=1,dofs
        DO i=1,n
          pi = A % InvPerm(i)
          IF(pi == 0 .OR. pi > m) CYCLE
          acti = IsInvInvPerm(pi)
          IF(ASSOCIATED(Var % Perm)) THEN
            IF(pi > SIZE(Var % Perm)) CYCLE
            pi = Var % Perm(pi)
            IF(pi==0) CYCLE
          END IF
          pi = dofs*(pi-1)+k
          r1 = 0.0_dp
          DO j=Rows(i),Rows(i+1)-1
            pj = Cols(j)
            IF(ASSOCIATED(Var % Perm)) THEN
              IF(pj > SIZE(Var % Perm)) CYCLE
              pj = Var % Perm(pj)
              IF(pj==0) CYCLE
            END IF
            pj = dofs*(pj-1)+k
            actj = IsInvInvPerm(Cols(j))
            
            IF( ToSlave ) THEN
              IF(.NOT. actj) THEN
                ! Project only master dofs to slave.
                r1 = -Values(j) * Var % Values(pj)
                TargetVar % Values(pi) = TargetVar % Values(pi) + r1
              END IF
            ELSE IF(.NOT. ToMaster) THEN
              ! Either project fully to slave
              r1 = -Values(j) * Var % Values(pj)
              TargetVar % Values(pi) = TargetVar % Values(pi) + r1
            END IF
            ! Normalize using just the weights from slave dofs.
            IF( Normalize .AND. actj) Weight(pi) = Weight(pi) + Values(j)                

            IF( ToMaster ) THEN
              ! Project slave dofs to master
              IF(.NOT. actj) THEN
                r1 = Values(j) * Var % Values(pi)
                TargetVar % Values(pj) = TargetVar % Values(pj) + r1               
                ! Normalize using just the weights from master dofs.
                IF(Normalize .AND. .NOT. actj) Weight(pj) = Weight(pj) + Values(j)
              END IF
            END IF                        
          END DO          
        END DO
      END DO
      DEALLOCATE(IsInvInvPerm)       
    END DO

    IF(Normalize) THEN
      WHERE(ABS(Weight) >  EPSILON(r1) )
        TargetVar % Values = TargetVar % Values / Weight
      END WHERE
    END IF
    
  END SUBROUTINE ProjectToVariable
  
END SUBROUTINE SaveProjection

!> \}
