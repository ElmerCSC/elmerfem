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
!
! ******************************************************************************
! *
! *  Authors: Iain Wheel
! *
! *  Original Date: 6/4/22
! *
! ****************************************************************************/

! create test levelset variable for parallel remeshing cube test case
SUBROUTINE TestLevelSet ( Model, Solver, dt, TransientSimulation )
  
   USE DefUtils

   IMPLICIT NONE

!-----------------------------------------------
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
!-----------------------------------------------
   TYPE(Mesh_t), POINTER :: Mesh
   TYPE(Variable_t), POINTER :: TimeVar, LSetVar
   REAL(KIND=dp), POINTER :: LSetValues(:)
   INTEGER, POINTER :: LSetPerm(:)
   INTEGER :: i,NNodes
   REAL(KIND=dp) :: TargetX,xx
   
   Mesh => Model % Mesh
   NNodes = Mesh % NumberOfNodes

   TimeVar => VariableGet( Model % Mesh % Variables, 'Timestep', ThisOnly=.TRUE.)
   IF(.NOT. ASSOCIATED(TimeVar)) THEN
     CALL Fatal('TestLevelSet','Could not find "timestep" variable!')
   END IF
   
   TargetX = TimeVar % Values(1) / 10

   LSetVar => Solver % Variable
   IF(.NOT. ASSOCIATED(LSetVar)) CALL FATAL('TestLevelSet', 'No variable associated with solver')
   LSetValues => LSetVar % Values
   LSetPerm => LSetVar % Perm
   
   DO i=1, NNodes
     xx = Mesh % Nodes % x(i)
     LSetValues(LSetPerm(i)) = xx - TargetX
   END DO

   !NULLIFY(LSetPerm, LSetValues)

END SUBROUTINE TestLevelSet
