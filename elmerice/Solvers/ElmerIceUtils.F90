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
! *  Authors: Mondher Chekki 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 15/12/2018
!
! ****************************************************************************/
! -----------------------------------------------------------------------

MODULE ElmerIceUtils 

USE DefUtils

CONTAINS 

SUBROUTINE  ComputeWeight(Model, Solver, VarName, WeightIn)

!------------------------------------------------------------------------------
!******************************************************************************
!
!   Compute weight at Boundaries (if Weight is not associated) 
!   and Sum over all Partitions   
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
!     INPUT: Name of the variable
!
!  TYPE(Variable_t) :: WeightIn
!     INPUT/OUTPUT : Variable Associated Weight
!
!******************************************************************************
    IMPLICIT NONE

     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     CHARACTER(LEN=MAX_NAME_LEN) :: VarName
     TYPE(Variable_t), POINTER ::  WeightIn 

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName

     LOGICAL :: UnFoundFatal=.TRUE. 
     INTEGER :: nlen

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif

        FunctionName='ComputeWeight'

!------------------------------------------------------------------------------
!     Compute Boundary Weights 
!------------------------------------------------------------------------------
        nlen = LEN_TRIM(VarName)

        NULLIFY( WeightIn )

        CALL CalculateNodalWeights(Model % Solver,.TRUE.)
        WeightIn       => VariableGet( Solver % Mesh % Variables, &
                      VarName(1:nlen)//' Boundary Weights' ,UnFoundFatal=UnFoundFatal)
        IF( .NOT. ASSOCIATED( WeightIn ) ) THEN
          CALL Fatal('ComputeWeight','Weight variable not present?')
        END IF

        CALL INFO(FunctionName, 'All Done', level=3)

END SUBROUTINE  ComputeWeight

SUBROUTINE  UpdatePartitionWeight(Model, Solver, VarName, Force, Force_tmp)

!------------------------------------------------------------------------------
!******************************************************************************
!
!   Compute weight at Boundaries (if Weight is not associated) 
!   and Sum over all Partitions   
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
!     INPUT: Name of the variable
!
!  TYPE(Variable_t) :: WeightIn
!     INPUT/OUTPUT : Variable Associated Weight
!
!******************************************************************************
    IMPLICIT NONE

     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     CHARACTER(LEN=MAX_NAME_LEN) :: VarName
     TYPE(Variable_t), POINTER :: Force
     REAL(KIND=dp) :: Force_tmp(:)


     INTEGER, POINTER :: WeightPerm(:), ForcePerm(:)
     INTEGER          :: i, nlen,  istat
     REAL(KIND=dp), POINTER :: WeightValues(:), ForceValues(:)
     REAL(KIND=dp) , ALLOCATABLE :: WeightValues_tmp(:)
     TYPE(Variable_t), POINTER ::  WeightIn 

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName

     LOGICAL :: UnFoundFatal=.TRUE. 

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif

        FunctionName='UpdatePartitionWeight'

        CALL ComputeWeight(Model, Solver, VarName, WeightIn)

        WeightPerm    => WeightIn % Perm
        WeightValues  => WeightIn % Values

        ForcePerm    => Force % Perm
        ForceValues  => Force % Values

!-----  -------------------------------------------------------------------------
!       Allocate  temporary storage
!-----  -------------------------------------------------------------------------

        ALLOCATE(WeightValues_tmp(SIZE(WeightValues)),STAT=istat)
        IF ( istat /= 0 ) THEN
           CALL Fatal( FunctionName, 'Memory allocation error.' )
        END IF

        WeightValues_tmp(:)  = WeightValues(:)
        Force_tmp(:)   = ForceValues(:)

        IF ( ParEnv % PEs >1 ) THEN 
           CALL ParallelSumVector( Solver % Matrix, WeightValues)
           DO i=1, Model % Mesh % NumberOfNodes
              IF (WeightPerm(i)>0) THEN
                Force_tmp(ForcePerm(i)) = Force_tmp(ForcePerm(i)) * &
                                             (WeightValues_tmp(WeightPerm(i))/WeightValues(WeightPerm(i)))
              END IF
           END DO
        END IF

        CALL INFO(FunctionName, 'All Done', level=3)

END SUBROUTINE UpdatePartitionWeight 


SUBROUTINE  UpdatePeriodicNodes(Model, Solver, VarName, WeightIn, ThisDim)

!------------------------------------------------------------------------------
!******************************************************************************
!
!  Update values at Periodic nodes 
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
!     INPUT: Name of the variable
!
!  TYPE(Variable_t) :: WeightIn
!     INPUT/OUTPUT : Variable Associated Weight
!
!  INTEGER        :: ThisDim
!     INPUT       : Variable Component X/Y/Z
!
!******************************************************************************


     IMPLICIT NONE

     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     CHARACTER(LEN=MAX_NAME_LEN) :: VarName
     TYPE(Variable_t), POINTER ::  WeightIn 
     INTEGER    ::  ThisDim 

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Weight 
     REAL(KIND=dp), POINTER :: WeightValues(:)
     INTEGER, POINTER       :: WeightPerm(:)
     CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName

     TYPE(Matrix_t), POINTER :: Projector
     LOGICAL, ALLOCATABLE :: ActivePart(:)

     INTEGER :: iBC, ii, jj, i ,j , k, nlen, l 
     INTEGER :: PeriodicNode1, PeriodicNode2
     INTEGER :: LocalPeriodicNode1, LocalPeriodicNode2
     INTEGER :: VarDim 
     REAL(KIND=dp) :: isPeriodic, tmpValue

     LOGICAL :: UnFoundFatal=.TRUE., GotIt  

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif

     FunctionName='UpdatePeriodicNodes'


     NULLIFY( Weight )

     Weight => WeightIn

     WeightPerm    => Weight % Perm
     WeightValues  => Weight % Values

     nlen   = LEN_TRIM(VarName)
     VarDim = Weight % DOFs

     CALL INFO("Update Periodic Nodes for: "//VarName(1:nlen), 'Start', level=3)
!---------------------------------------------------------------------------
!    date Values at Periodic Nodes using a Projector
!---------------------------------------------------------------------------
     ALLOCATE(ActivePart(MAX( Model % NumberOfBodyForces,Model % NumberOfBCs)))
  
     ActivePart = .FALSE.

     DO iBC=1,Model % NumberOfBCs
       IF ( ListGetLogical( Model % BCs(iBC) % Values, &
           'Periodic BC ' // VarName(1:nlen), GotIt ) ) ActivePart(iBC) = .TRUE.  
       IF ( ListGetLogical( Model % BCs(iBC) % Values, &
           'Anti Periodic BC ' // VarName(1:nlen), GotIt ) ) ActivePart(iBC) = .TRUE.
       IF ( ListCheckPresent( Model % BCs(iBC) % Values, &
           'Periodic BC Scale ' // VarName(1:nlen) ) ) ActivePart(iBC) = .TRUE.
     END DO

     DO iBC=1,Model % NumberOfBCs
            IF (ActivePart(iBC)) THEN
               Projector => Model % BCs(iBC) % PMatrix
               IF ( ASSOCIATED( Projector ) ) THEN 
                   DO i=1,Projector % NumberOfRows
                      DO l = Projector % Rows(i), Projector % Rows(i+1)-1

                         PeriodicNode1 = Projector % InvPerm(i)
                         PeriodicNode2 = Projector % Cols(l)
                         isPeriodic = Projector % Values(l)

                         IF ( PeriodicNode1 <= 0 .OR. PeriodicNode2 <= 0 ) CYCLE

                         LocalPeriodicNode1 = WeightPerm(PeriodicNode1)
                         LocalPeriodicNode2 = WeightPerm(PeriodicNode2)

                         IF ( LocalPeriodicNode1>0 .and. LocalPeriodicNode2>0 ) THEN

                            LocalPeriodicNode1=VarDim*(LocalPeriodicNode1-1)+ThisDim
                            LocalPeriodicNode2=VarDim*(LocalPeriodicNode2-1)+ThisDim
 
                            tmpValue=WeightValues(LocalPeriodicNode1)

                            WeightValues(LocalPeriodicNode1) = tmpValue + isPeriodic*WeightValues(LocalPeriodicNode2)
                            WeightValues(LocalPeriodicNode2) = WeightValues(LocalPeriodicNode2) + isPeriodic*tmpValue

                         ENDIF !LocalPeriodicNode >0 
                      END DO !l
                   END DO !i
                END IF !Projector 
            END IF !ActivePart
  
     END DO !iBC

     CALL INFO(FunctionName, 'All Done', level=3)

END  SUBROUTINE  UpdatePeriodicNodes

SUBROUTINE  SetZeroAtPeriodicNodes(Model, Solver, VarName, WeightValues, WeightPerm, TargetPerm)

!------------------------------------------------------------------------------
!******************************************************************************
!
!  Set 0 values at Periodic nodes 
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
!     INPUT: Name of the variable
!
!******************************************************************************


     IMPLICIT NONE

     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     CHARACTER(LEN=MAX_NAME_LEN) :: VarName
     REAL(KIND=dp)  :: WeightValues(:)
     INTEGER       :: WeightPerm(:)
     INTEGER       :: TargetPerm(:)

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Weight 
     CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName

     TYPE(Matrix_t), POINTER :: Projector
     LOGICAL, ALLOCATABLE :: ActivePart(:)

     INTEGER :: iBC, ii, jj, i ,j , k, nlen, l 
     INTEGER :: PeriodicNode1
     INTEGER :: LocalPeriodicNode1
     REAL(KIND=dp) :: isPeriodic

     LOGICAL ::  GotIt  

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif

     FunctionName='SetZeroAtPeriodicNodes'

     nlen   = LEN_TRIM(VarName)
!---------------------------------------------------------------------------
!    date Values at Periodic Nodes using a Projector
!---------------------------------------------------------------------------
     ALLOCATE(ActivePart(MAX( Model % NumberOfBodyForces,Model % NumberOfBCs)))
  
     ActivePart = .FALSE.

     DO iBC=1,Model % NumberOfBCs
       IF ( ListGetLogical( Model % BCs(iBC) % Values, &
           'Periodic BC ' // VarName(1:nlen), GotIt ) ) ActivePart(iBC) = .TRUE.  
       IF ( ListGetLogical( Model % BCs(iBC) % Values, &
           'Anti Periodic BC ' // VarName(1:nlen), GotIt ) ) ActivePart(iBC) = .TRUE.
       IF ( ListCheckPresent( Model % BCs(iBC) % Values, &
           'Periodic BC Scale ' // VarName(1:nlen) ) ) ActivePart(iBC) = .TRUE.
     END DO

     DO iBC=1,Model % NumberOfBCs
            IF (ActivePart(iBC)) THEN
               Projector => Model % BCs(iBC) % PMatrix
               IF ( ASSOCIATED( Projector ) ) THEN 
                   DO i=1,Projector % NumberOfRows
                      DO l = Projector % Rows(i), Projector % Rows(i+1)-1

                         PeriodicNode1 = Projector % InvPerm(i)
                       
                         isPeriodic = Projector % Values(l)

                         IF ( PeriodicNode1 <= 0 ) CYCLE

                         LocalPeriodicNode1 = TargetPerm(PeriodicNode1)

                         IF ( WeightPerm(PeriodicNode1) >  0 )  WeightValues(LocalPeriodicNode1) = 0.0D0 

                      END DO !l
                   END DO !i
                END IF !Projector 
            END IF !ActivePart
  
     END DO !iBC

     CALL INFO(FunctionName, 'All Done', level=3)

END  SUBROUTINE SetZeroAtPeriodicNodes 

END MODULE ElmerIceUtils 
