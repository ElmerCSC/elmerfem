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
! *
! *  Authors: Olivier Gagliardini, Ga¨el Durand
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!> ZsMZsIni for variable Zs                   
!> ZsTopMZsIni for variable Zs Top         
!> ZsBottomMZsIni for variable Zs Bottom
!> DyMDyIni for any FS variable name           
FUNCTION ZsIni ( Model, nodenumber, x) RESULT(Zs)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsIni','Could not find variable >Zs<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

       Zs = Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsIni

FUNCTION ZsMzsIni ( Model, nodenumber, Zs) RESULT(mu)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: mu,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsTopMZsIni','Could not find variable >Zs<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

      mu =  Zs -  Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsMZsIni

!--------------------------------------------------------------------------------

FUNCTION ZsTopIni ( Model, nodenumber, x) RESULT(Zs)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs Top')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsIni','Could not find variable >Zs Top<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

       Zs = Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsTopIni

FUNCTION ZsTopMzsIni ( Model, nodenumber, Zs) RESULT(mu)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: mu,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs Top')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsTopMZsIni','Could not find variable >Zs Top<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

      mu =  Zs -  Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsTopMZsIni

!--------------------------------------------------------------------

FUNCTION ZsBottomIni ( Model, nodenumber, x) RESULT(Zs)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs Bottom')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsIni','Could not find variable >Zs Bottom<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

       Zs = Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsBottomIni

FUNCTION ZsBottomMzsIni ( Model, nodenumber, Zs) RESULT(mu)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: mu,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs Bottom')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsTopMZsIni','Could not find variable >Zs Bottom<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

      mu =  Zs -  Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsBottomMZsIni

!---------------------------------------------------------------------------

FUNCTION DyIni ( Model, nodenumber, x) RESULT(Zs)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x,   Zs   
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   IF (FirstTime) THEN
          FirstTime = .FALSE.
          NMAX = Model % NumberOfNodes             
          dim = CoordinateSystemDimension()
          ALLOCATE(Zs0(NMAX))
          DO i = 1, NMAX
            IF (dim==2) THEN
               Zs0(i) = Model % Nodes % y (i)
            ELSE
               Zs0(i) = Model % Nodes % z (i)
            END IF
          END DO
   END IF

       Zs = Zs0(nodenumber) 

END FUNCTION DyIni


FUNCTION DyMDyIni ( Model, nodenumber, Zs) RESULT(mu)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: mu,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   IF (FirstTime) THEN
          FirstTime = .FALSE.
          NMAX = Model % NumberOfNodes             
          dim = CoordinateSystemDimension()
          ALLOCATE(Zs0(NMAX))
          DO i = 1, NMAX
            IF (dim==2) THEN
               Zs0(i) = Model % Nodes % y (i)
            ELSE
               Zs0(i) = Model % Nodes % z (i)
            END IF
          END DO
   END IF

       mu = Zs - Zs0(nodenumber)

END FUNCTION DyMDyIni

