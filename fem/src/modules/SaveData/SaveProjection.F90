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
! *  Original Date: 20 Nov 2001
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
  LOGICAL :: Found 

  CALL Info('SaveProjection','Creating selected projected values as fields')

  Params => GetSolverParams()
  
  i = 0
  DO WHILE(.TRUE.)  
    i = i + 1
    VarName = ListGetString( Params,'Variable '//I2S(i), Found )
    IF(.NOT. Found) EXIT    
  END DO
  NoVar = i-1
  CALL Info('SaveProjection','Saving projections from '//I2S(NoVar)//' fields')
  

  DO i=1,NoVar
    VarName = ListGetString( Params,'Variable '//I2S(i), Found )
    Var => VariableGet( Model % Variables, TRIM(VarName) )
    IF(.NOT. ASSOCIATED(Var)) THEN
      CALL Warn('SaveProjection','Requested variable does not exist!')
      CYCLE
    END IF
    
    TargetName = ListGetString( Params,'Target Variable '//I2S(i), Found )
    IF(.NOT. Found) VarName = 'Projection '//TRIM(VarName)
    
    TargetVar => VariableGet( Model % Variables, TRIM(TargetName) )
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
    INTEGER :: bc, dofs, i, j, k, pi, pj
    INTEGER, POINTER :: Rows(:), Cols(:)
    REAL(KIND=dp) :: r1
    REAL(KIND=dp), POINTER :: Values(:)      

    dofs = Var % Dofs
    TargetVar % Values = 0.0_dp

    ! Go through all the projectors.
    ! There could be perhaps reason to skip some, but this will do for now.
    DO bc=1,Model % NumberOfBCs        
      A => CurrentModel % BCs(bc) % PMatrix
      IF(.NOT. ASSOCIATED(A) ) THEN
        A => Solver % MortarBCs(bc) % Projector
      END IF
      
      IF(.NOT. ASSOCIATED(A)) CYCLE
      n = A % NumberOfRows
      Rows => A % Rows
      Cols => A % Cols
      Values => A % Values
      
      DO k=1,dofs
        DO i=1,n
          IF(ASSOCIATED(A % InvPerm)) THEN
            pi = A % InvPerm(i)
          ELSE
            pi = i
          END IF
          IF(ASSOCIATED(Var % Perm)) pi = Var % Perm(pi)          
          IF(pi==0) CYCLE
          pi = dofs*(pi-1)+k
          r1 = 0.0_dp
          DO j=Rows(i),Rows(i+1)-1
            pj = Cols(j)
            IF(ASSOCIATED(Var % Perm)) pj = Var % Perm(pj)
            IF(pj==0) CYCLE
            pj = dofs*(pj-1)+k
            r1 = r1 + Values(j) * Var % Values(pj)
          END DO
          
          TargetVar % Values(pi) = TargetVar % Values(pi) + r1 
        END DO
      END DO
    END DO
  END SUBROUTINE ProjectToVariable
  
END SUBROUTINE SaveProjection

!> \}
