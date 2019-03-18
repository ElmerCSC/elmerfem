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

!------------------------------------------------------------------------------
!> Saves results in a format understood by the pre-/postprocessing software GiD.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE GiDOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: Solution
  INTEGER, POINTER :: Perm(:)
  REAL(KIND=dp), POINTER :: Values(:)
  COMPLEX(KIND=dp), POINTER :: CValues(:)
  TYPE(Variable_t), POINTER :: TimeVariable
  TYPE(ValueList_t), POINTER :: SolverParams

  LOGICAL :: Found, CoordinatesWritten = .FALSE.
  LOGICAL :: EigenAnalysis = .FALSE., FirstTimeStep

  INTEGER :: i,j,k,m,n,dim, Code, body_id, ElementCounter, Nloop, Loop, ExtCount
  INTEGER :: tensorComponents, SteadyStep = 0
  INTEGER, PARAMETER :: MaxElemCode = 1000
  INTEGER :: ListElemTypes(MaxElemCode)


  CHARACTER(LEN=1024) :: OutputFile, ResFile, MshFile, Txt, Family, &
       ScalarFieldName, VectorFieldName, TensorFieldName, CompName
  CHARACTER(LEN=1024) :: Txt2, Txt3

  INTEGER :: PyramidMap613(14,4), PyramidMap605(2,4)
  INTEGER :: WedgeMap706(3,4), WedgeMap715(21,4)
  INTEGER :: GidUnit
  
!------------------------------------------------------------------------------

  PyramidMap605(1,:) = (/ 3, 5, 4, 1 /)
  PyramidMap605(2,:) = (/ 3, 5, 2, 1 /)

  PyramidMap613(1,:)  = (/ 7, 8, 3, 12 /)
  PyramidMap613(2,:)  = (/ 10, 11, 12, 5 /)
  PyramidMap613(3,:)  = (/ 10, 13, 12, 5 /)
  PyramidMap613(4,:)  = (/ 9, 10, 13, 11 /)
  PyramidMap613(5,:)  = (/ 9, 13, 11, 12 /)
  PyramidMap613(6,:)  = (/ 9, 10, 11, 1 /)
  PyramidMap613(7,:)  = (/ 9, 6, 11, 1 /)
  PyramidMap613(8,:)  = (/ 9, 6, 11, 12 /)
  PyramidMap613(9,:)  = (/ 9, 8, 12, 4 /)
  PyramidMap613(10,:) = (/ 9, 13, 12, 4 /)
  PyramidMap613(11,:) = (/ 7, 9, 8, 12 /)
  PyramidMap613(12,:) = (/ 7, 9, 6, 12 /)
  PyramidMap613(13,:) = (/ 7, 6, 11, 12 /)
  PyramidMap613(14,:) = (/ 7, 6, 11, 2 /)

  WedgeMap706(1,:) = (/ 5, 4, 3, 1 /)
  WedgeMap706(2,:) = (/ 5, 3, 2, 1 /)
  WedgeMap706(3,:) = (/ 5, 6, 4, 3 /)

  WedgeMap715(1,:) = (/ 10, 11, 5, 2 /)
  WedgeMap715(2,:) = (/ 12, 11, 6, 3 /)
  WedgeMap715(3,:) = (/ 12, 10, 4, 1 /)
  WedgeMap715(4,:) = (/ 7, 8, 11, 2 /)
  WedgeMap715(5,:) = (/ 7, 10, 11, 2 /)
  WedgeMap715(6,:) = (/ 13, 14, 11, 5 /)
  WedgeMap715(7,:) = (/ 13, 10, 11, 5 /)
  WedgeMap715(8,:) = (/ 9, 10, 8, 11 /)
  WedgeMap715(9,:) = (/ 9, 7, 10, 8 /)
  WedgeMap715(10,:) = (/ 9, 12, 10, 11 /)
  WedgeMap715(11,:) = (/ 9, 12, 10, 1 /)
  WedgeMap715(12,:) = (/ 9, 7, 10, 1 /)
  WedgeMap715(13,:) = (/ 9, 8, 11, 3 /)
  WedgeMap715(14,:) = (/ 9, 12, 11, 3 /)
  WedgeMap715(15,:) = (/ 15, 12, 10, 11 /)
  WedgeMap715(16,:) = (/ 15, 12, 10, 4 /)
  WedgeMap715(17,:) = (/ 15, 10, 14, 11 /)
  WedgeMap715(18,:) = (/ 15, 13, 10, 14 /)
  WedgeMap715(19,:) = (/ 15, 13, 10, 4 /)
  WedgeMap715(20,:) = (/ 15, 14, 11, 6 /)
  WedgeMap715(21,:) = (/ 15, 12, 11, 6 /)

  SolverParams => GetSolverParams()
  EigenAnalysis = GetLogical( SolverParams, 'Eigen Analysis', Found )

  
  ExtCount = ListGetInteger( Solver % Values,'Output Count',Found)
  IF( Found ) THEN
    SteadyStep = ExtCount 
  ELSE
    SteadyStep = SteadyStep + 1
  END IF
  FirstTimeStep = ( SteadyStep == 1 )


  OutputFile = GetString( Solver % Values, 'Output File Name', Found )
  IF(.NOT. Found) OutputFile = 'Output'

  WRITE(ResFile,'(A,A)') TRIM(OutputFile),'.flavia.res'
  WRITE(MshFile,'(A,A)') TRIM(OutputFile),'.flavia.msh'


  CALL Info('GidOutputSolver','Writing result for GiD postprocessing')
  CALL Info('GidOutputSolver','res-file = :'//TRIM(ResFile) )
  CALL Info('GidOutputSolver','msh-file = :'//TRIM(MshFile) )


  ! Write the GiD msh-file:
  !------------------------
  dim = CoordinateSystemDimension()
  IF( CoordinatesWritten ) GOTO 10

  OPEN(NEWUNIT=GidUnit, FILE=MshFile )

  ! First check how many element types are involved in the analysis:
  !-----------------------------------------------------------------
  ListElemTypes = 0
  ElementCounter = 0
!  DO i = 1, Solver % NumberOfActiveElements
!     Element => GetActiveElement(i)
  DO i = 1, Model % NumberOfBulkElements !+ Model % NumberOfBoundaryElements
     Element => Model % Mesh % Elements(i)

     Code = Element % TYPE % ElementCode
     ListElemTypes(Code) = ListElemTypes(Code)+1
  END DO
  PRINT *,'Total number of elements =',SUM(ListElemTypes)

  ! Write the different element types in different blocks:
  !-------------------------------------------------------
  DO i = 1,MaxElemCode
     IF(ListElemTypes(i) == 0) CYCLE
     PRINT *,ListElemTypes(i),'elements of type',i
     n = MOD(i,100)
     IF( INT(i/100) == 1 ) Family = 'Point'
     IF( INT(i/100) == 2 ) Family = 'Linear'
     IF( INT(i/100) == 3 ) Family = 'Triangle'
     IF( INT(i/100) == 4 ) Family = 'Quadrilateral'
     IF( INT(i/100) == 5 ) Family = 'Tetrahedra'
     IF( INT(i/100) == 6 ) THEN
        Family = 'Tetrahedra' ! PYRAMIDS WILL BE SPLITTED
        n = 4                 ! INTO LINEAR TETRAHEDRA
     END IF
     IF( INT(i/100) == 7 ) THEN
        Family = 'Tetrahedra' ! WEDGES WILL BE SPLITTED
        n = 4                 ! INTO LINEAR TETRAHEDRA
     END IF
     IF( INT(i/100) == 8 ) Family = 'Hexahedra'

     WRITE(Txt,'(A,I1,A,A,A,I0)') 'MESH "Elmer Mesh" dimension ',&
             dim,' ElemType ', TRIM(Family),' Nnode ', n

     WRITE(GidUnit,'(A)') TRIM(Txt)

     ! Write all node coordinates in the first block:
     !-----------------------------------------------
     IF( .NOT.CoordinatesWritten ) THEN
        WRITE(GidUnit,'(A)') 'Coordinates'
        DO j = 1, Model % Mesh % NumberOfNodes
           WRITE(GidUnit,'(I6,3ES16.7E3)') j, &
                Model % Mesh % Nodes % x(j), &
                Model % Mesh % Nodes % y(j), &
                Model % Mesh % Nodes % z(j)
        END DO
        WRITE(GidUnit,'(A)') 'end coordinates'
        WRITE(GidUnit,'(A)') ' '
        CoordinatesWritten = .TRUE.
     END IF

     ! Write the element connectivity tables:
     !---------------------------------------
     WRITE(GidUnit,'(A)') 'Elements'
!     DO j = 1, Solver % NumberOfActiveElements
!        Element => GetActiveElement( j )
     DO j = 1, Model % NumberOfBulkElements !+ Model % NumberOfBoundaryElements
        Element => Model % Mesh % Elements(j)

        Code = Element % TYPE % ElementCode
        IF( Code /= i ) CYCLE
        body_id = Element % BodyId

        IF( Code == 613 ) THEN
           ! 13 noded pyramids will be splitted into 14 linear tetraheda
           !------------------------------------------------------------
           DO m = 1,14
              ElementCounter = ElementCounter + 1
              WRITE(GidUnit,'(100I10)') ElementCounter, &
                   Element % NodeIndexes(PyramidMap613(m,:)), body_id
           END DO
        ELSEIF( Code == 605 ) THEN
           ! 5 noded pyramids will be splitted into 2 linear tetraheda
           !----------------------------------------------------------
           DO m = 1,2
              ElementCounter = ElementCounter + 1
              WRITE(GidUnit,'(100I10)') ElementCounter, &
                   Element % NodeIndexes(PyramidMap605(m,:)), body_id
           END DO           
        ELSEIF( Code == 706 ) THEN
           ! 6 noded wedges will be splitted into 3 linear tetraheda
           !---------------------------------------------------------
           DO m = 1,3
              ElementCounter = ElementCounter + 1
              WRITE(GidUnit,'(100I10)') ElementCounter, &
                   Element % NodeIndexes(WedgeMap706(m,:)), body_id
           END DO           
        ELSEIF( Code == 715 ) THEN
           ! 15 noded wedges will be splitted into 21 linear tetraheda
           !----------------------------------------------------------
           DO m = 1,21
              ElementCounter = ElementCounter + 1
              WRITE(GidUnit,'(100I10)') ElementCounter, &
                   Element % NodeIndexes(WedgeMap715(m,:)), body_id
           END DO           
        ELSE
           ! Standard elements are understood by GiD as such
           !------------------------------------------------
           ElementCounter = ElementCounter + 1 
           WRITE(GidUnit,'(100I10)') ElementCounter, Element % NodeIndexes, body_id
        END IF

     END DO
     WRITE(GidUnit,'(A)') 'end elements'
     WRITE(GidUnit,'(A)') ' '
  END DO

10 CONTINUE
  
  ! Write the GiD res-file:
  !------------------------
  IF( FirstTimeStep ) THEN
    OPEN(NEWUNIT=GidUnit, FILE=ResFile )
    WRITE(GidUnit,'(A)') 'GiD Post Result File 1.0'
  ELSE
    OPEN(NEWUNIT=GidUnit, FILE=ResFile, POSITION='APPEND' )
  END IF

  Nloop = 1
  IF( EigenAnalysis ) THEN
     Nloop = GetInteger( Solver % Values, 'Eigen System Values', Found )
     IF( .NOT. Found ) Nloop = GetInteger( Solver % Values, 'Number Of EigenModes', Found )
     IF( .NOT.Found ) Nloop = 1
   END IF

   DO Loop = 1, Nloop
     PRINT *,'------------'

     ! First scalar fields:
     !----------------------
     DO i = 1, 999
       WRITE(Txt,'(A,I0)') 'Scalar Field ',i
       
       ScalarFieldName = GetString( Solver % Values, TRIM(Txt), Found )
       IF(.NOT. Found) EXIT
       
       Solution => VariableGet( Model % Mesh % Variables, ScalarFieldName )
       IF( .NOT.ASSOCIATED( Solution ) ) THEN
         PRINT *,'Scalar field "',TRIM(ScalarFieldName),'" not found'
       ELSE
         PRINT *,'Scarar field',i,'= "',TRIM(ScalarFieldName),'"'
         Perm => Solution % Perm
         
         IF( .NOT.EigenAnalysis ) THEN
           Values => Solution % Values
         ELSE
           Cvalues => Solution % EigenVectors(Loop,:)
         END IF
         
         IF( TransientSimulation ) THEN
           TimeVariable => VariableGet( Model % Mesh % Variables, 'Time' )
           PRINT *,'Current time=',TimeVariable % Values(1)
           WRITE(GidUnit,'(A,A,A,ES16.7E3,A)') 'Result "',&
               TRIM(ScalarFieldName),'" "Transient analysis" ', &
               TimeVariable % Values(1) ,' Scalar OnNodes'
         ELSE
           IF( .NOT.EigenAnalysis ) THEN
             WRITE(GidUnit,'(A,A,A,I2,A)') 'Result "',&
                 TRIM(ScalarFieldName),'" "Steady analysis"',SteadyStep,' Scalar OnNodes'
           ELSE
             WRITE(GidUnit,'(A,A,A,I2,A)') 'Result "',&
                 TRIM(ScalarFieldName),'" "Eigen analysis"',Loop,' Scalar OnNodes'
           END IF
         END IF
         
         WRITE(GidUnit,'(A,A,A)') 'ComponentNames "',TRIM(ScalarFieldName),'"'
         WRITE(GidUnit,'(A)') 'Values'
         DO j = 1, Model % Mesh % NumberOfNodes
           k = Perm(j)
           IF( .NOT.EigenAnalysis ) THEN
             WRITE(GidUnit,'(I6,ES16.7E3)') j, Values(k)
           ELSE
             WRITE(GidUnit,'(I6,ES16.7E3)') j, REAL(CValues(k))
           END IF
         END DO
         WRITE(GidUnit,'(A)') 'end values'
         WRITE(GidUnit,'(A)') ' '
       END IF
     END DO

     ! Then vector fields:
     !--------------------
     DO i = 1, 999
       WRITE(Txt,'(A,I0)') 'Vector Field ',i
       
       VectorFieldName = GetString( Solver % Values, TRIM(Txt), Found )
       IF(.NOT. Found) EXIT
       PRINT *,'Vector field',i,'= "',TRIM(VectorFieldName),'"'
       
       IF( TransientSimulation ) THEN
         TimeVariable => VariableGet( Model % Mesh % Variables, 'Time' )
         PRINT *,'Current time=',TimeVariable % Values(1)
         WRITE(GidUnit,'(A,A,A,ES16.7E3,A)') 'Result "',&
             TRIM(VectorFieldName),'" "Transient analysis" ', &
             TimeVariable % Values(1) ,' Vector OnNodes'
       ELSE
         IF( .NOT.EigenAnalysis ) THEN
           WRITE(GidUnit,'(A,A,A,I2,A)') 'Result "',&
               TRIM(VectorFieldName),'" "Steady analysis"',SteadyStep,' Vector OnNodes'
         ELSE
           WRITE(GidUnit,'(A,A,A,I2,A)') 'Result "',&
               TRIM(VectorFieldName),'" "Eigen analysis"',Loop,' Vector OnNodes'
         END IF
       END IF
       
       WRITE(Txt,'(A)') 'ComponentNames '
       DO j = 1, dim
         IF(j<Dim) THEN
           WRITE(Txt,'(A,A,A,I2,A)' ) &
               TRIM(Txt), ' "', TRIM(VectorFieldName),j,'",'
         ELSE
           WRITE(Txt,'(A,A,A,I2,A)' ) &
               TRIM(Txt), ' "', TRIM(VectorFieldName),j,'"'
         END IF
       END DO
       WRITE(GidUnit,'(A)') TRIM(Txt)
       WRITE(GidUnit,'(A)') 'Values'
       
       DO j = 1, Model % Mesh % NumberOfNodes
         WRITE(Txt2,'(I10)') j
         DO k = 1,dim
           
           ! Check if vector field components have been defined explicitly:
           !----------------------------------------------------------------
           WRITE(Txt3,'(A,I1,A,I1)') 'Vector Field ',i,' component ',k
           CompName = GetString( Solver % Values, TRIM(Txt3), Found )
           IF( Found ) THEN
             WRITE(Txt,'(A)') TRIM(CompName)
           ELSE
             WRITE(Txt,'(A,I2)') TRIM(VectorFieldName), k
           END IF
           
           IF( j==1 ) PRINT *, TRIM(Txt3),' = "', TRIM(Txt),'"'
           
           Solution => VariableGet( Model % Mesh % Variables, TRIM(Txt) )
           IF( .NOT.ASSOCIATED( Solution ) ) THEN
             PRINT *,'Vector field component',k,' not found'
           ELSE
             Perm => Solution % Perm
             
             IF( .NOT.EigenAnalysis ) THEN
               Values => Solution % Values
               WRITE(Txt2,'(A,ES16.7E3)') TRIM(Txt2), Values( Perm(j) )
             ELSE
               CValues => Solution % Eigenvectors(Loop,:)
               WRITE(Txt2,'(A,ES16.7E3)') TRIM(Txt2), REAL(CValues( Perm(j) ) )
             END IF
             
           END IF
         END DO
         WRITE(GidUnit,'(A)') TRIM(Txt2)
         
       END DO
       WRITE(GidUnit,'(A)') 'end values'
     END DO
     

     ! Finally tensor fields:
     !-----------------------
     DO i = 1, 999
       WRITE(Txt,'(A,I0)') 'Tensor Field ',i
       TensorFieldName = GetString( Solver % Values, TRIM(Txt), Found )
       IF(.NOT. Found) EXIT
       
       PRINT *,'Tensor field',i,'= "',TRIM(TensorFieldName),'"'
       
       IF( TransientSimulation ) THEN
         TimeVariable => VariableGet( Model % Mesh % Variables, 'Time' )
         PRINT *,'Current time=',TimeVariable % Values(1)
         WRITE(GidUnit,'(A,A,A,ES16.7E3,A)') 'Result "',&
             TRIM(TensorFieldName),'" "Transient analysis" ', &
             TimeVariable % Values(1) ,' Matrix OnNodes'
       ELSE
         IF( .NOT.EigenAnalysis ) THEN
           WRITE(GidUnit,'(A,A,A,I2,A)') 'Result "',&
               TRIM(TensorFieldName),'" "Steady analysis"',SteadyStep,' Matrix OnNodes'
         ELSE
           WRITE(GidUnit,'(A,A,A,I2,A)') 'Result "',&
               TRIM(TensorFieldName),'" "Eigen analysis"',Loop,' Matrix OnNodes'
         END IF
         
       END IF
       
       WRITE(Txt,'(A)') 'ComponentNames '
       IF( dim == 2 ) THEN
         TensorComponents = 3
       ELSE
         TensorComponents = 6
       END IF
       
       DO j = 1, TensorComponents
         WRITE(Txt3,'(A,I1,A,I1)') 'Tensor Field ',i,' component ',j
         CompName = GetString( Solver % Values, TRIM(Txt3), Found )

         IF( Found ) THEN
           WRITE(Txt2,'(A)') TRIM(CompName)
         ELSE
           WRITE(Txt2,'(A,A,I1)') TRIM(TensorFieldName), ' ', j
         END IF

         IF(j<TensorComponents) THEN
           WRITE(Txt,'(A,A,A,A)' ) &
               TRIM(Txt), ' "', TRIM(Txt2), '",'
         ELSE
           WRITE(Txt,'(A,A,A,A)' ) &
               TRIM(Txt), ' "', TRIM(txt2), '"'
         END IF
       END DO
       WRITE(GidUnit,'(A)') TRIM(Txt)
       WRITE(GidUnit,'(A)') 'Values'
       
       DO j = 1, Model % Mesh % NumberOfNodes
         WRITE(Txt2,'(I10)') j
         DO k = 1,TensorComponents
           
           ! Check if tensor field components have been defined explicitly:
           !----------------------------------------------------------------
           WRITE(Txt3,'(A,I1,A,I1)') 'Tensor Field ',i,' component ',k
           CompName = GetString( Solver % Values, TRIM(Txt3), Found )
           IF( Found ) THEN
             WRITE(Txt,'(A)') TRIM(CompName)
           ELSE
             WRITE(Txt,'(A,A,I1)') TRIM(TensorFieldName), ' ', k
           END IF
           
           IF( j==1 ) PRINT *, TRIM(Txt3),' = "', TRIM(Txt),'"'
           
           Solution => VariableGet( Model % Mesh % Variables, TRIM(Txt) )
           IF( .NOT.ASSOCIATED( Solution ) ) THEN
             PRINT *,'Tensor field component',k,' not found'
           ELSE
             Perm => Solution % Perm
             !Values => Solution % Values
             !WRITE(Txt2,'(A,ES16.7E3)') TRIM(Txt2), Values( Perm(j) )
             
             IF( .NOT.EigenAnalysis ) THEN
               Values => Solution % Values
               WRITE(Txt2,'(A,ES16.7E3)') TRIM(Txt2), Values( Perm(j) )
             ELSE
               CValues => Solution % Eigenvectors(Loop,:)
               WRITE(Txt2,'(A,ES16.7E3)') TRIM(Txt2), REAL(CValues( Perm(j) ) )
             END IF

           END IF
         END DO
         WRITE(GidUnit,'(A)') TRIM(Txt2)
         
       END DO
       WRITE(GidUnit,'(A)') 'end values'
     END DO
     
   END DO ! Nloop


   CLOSE(GidUnit)
  
   CALL Info('GitOutputSolver','Output complete',Level=7)

!------------------------------------------------------------------------------
 END SUBROUTINE GiDOutputSolver
!------------------------------------------------------------------------------
  
