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
! *  Module information:
! *  Authors: Peter Råback, Mika Malinen, Juha Ruokolainen
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 16.3.2016
! *
! *****************************************************************************/

!-----------------------------------------------------------------------------
SUBROUTINE CraigBamptonSolver_init( Model,Solver,dt,TransientSimulation )
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
  LOGICAL :: Found
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName


  IF( ListGetLogical( Solver % Values,'Calculate Matrix Norm',Found ) ) THEN
    SolverName = ListGetString( Solver % Values, 'Equation', Found )
    IF( .NOT. ListCheckPresent( Solver % Values,'Variable') ) THEN
      CALL ListAddString( Solver % Values,'Variable',&
          '-nooutput -global '//TRIM(SolverName)//'_var')
    END IF
  END IF

END SUBROUTINE CraigBamptonSolver_init
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!> Craig-Bamption model reduction solver.
!------------------------------------------------------------------------------
SUBROUTINE CraigBamptonSolver( Model,Solver,dt,TransientSimulation )
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
  TYPE(ValueList_t), POINTER :: Params
  CHARACTER(LEN=MAX_NAME_LEN) :: Str, FileName
  INTEGER :: i,j,k,l,SolverId,Dofs,NoComponentModes,&
      NoEigenModes,NoConstraintModes, MatrixNo
  LOGICAL :: Found, SaveThis
  TYPE(Solver_t), POINTER :: ESolver
  TYPE(Variable_t), POINTER :: EVar
  REAL(KIND=dp), POINTER :: x(:)
  REAL (KIND=dp), POINTER CONTIG :: SaveValues(:)
  REAL(KIND=dp), ALLOCATABLE :: Ahat(:,:)
  REAL(KIND=dp), ALLOCATABLE :: Ax(:)
  TYPE(Matrix_t), POINTER :: A
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: val, Norm
!------------------------------------------------------------------------------

  CALL Info('CraigBamptonSolver','Performing Craig-Bampton model reduction')
  
  Params => GetSolverParams()

  str = ListGetString( Params,'Elasticity Solver Name',Found ) 
  IF( .NOT. Found ) THEN
    CALL Fatal('CraigBamptonSolver','> Elasticity Solver Name < not given!')
  END IF

  SolverId = 0
  DO i=1,Model % NumberOfSolvers
    IF( ListGetString( Model % Solvers(i) % Values,'Equation',Found ) == str ) THEN
      SolverId = i
      EXIT
    END IF
  END DO
  IF( SolverId == 0 ) THEN
    CALL Fatal('CraigBamptonSolver','No solver with Equation name: '//TRIM(str))
  END IF

  ESolver => Model % Solvers(SolverId) 
  EVar => ESolver % Variable

  A => ESolver % Matrix
  IF(.NOT. ASSOCIATED( A ) ) THEN
   CALL Fatal('CraigBamptonSolver','Matrix not associated!')
  END IF

  SaveValues => A % Values
  IF(.NOT. ASSOCIATED( A % BulkValues ) ) THEN
   CALL Fatal('CraigBamptonSolver','Bulk values of matrix not associated!')
  END IF

  A % Values => A % BulkValues
  x => EVar % Values
  Dofs = Evar % Dofs
  Mesh => ESolver % Mesh

  NoEigenModes = ESolver % NOFEigenValues
  IF( NoEigenModes == 0 ) THEN
    CALL Fatal('CraigBamptonSolver','Number of eigenmodes is zero!')
  END IF
  CALL Info('CraigBamptonSolver','Number of eigen modes: '&
      //TRIM(I2S(NoEigenModes)),Level=7)

  NoConstraintModes = EVar % NumberOfConstraintModes
  IF( SolverId == 0 ) THEN
    CALL Fatal('CraigBamptonSolver','Number of constraint modes is zero!')
  END IF
  CALL Info('CraigBamptonSolver','Number of constraint modes: '&
      //TRIM(I2S(NoConstraintModes)),Level=7)
  NoComponentModes = NoEigenModes + NoConstraintModes
  

  CALL SaveReductionDofs()

  SaveThis = ListGetLogical( Solver % Values,'Save Reduction Basis',Found)
  IF( SaveThis ) THEN
    FileName = 'ReductionModes.dat'
    CALL SaveReductionBasis()
  END IF

  SaveThis = ListGetLogical( Solver % Values,'Save Reduction Matrix',Found)
  IF(.NOT. Found ) SaveThis = .TRUE.

  IF( SaveThis ) THEN
    ALLOCATE( Ahat(NoComponentModes,NoComponentModes))
    Ahat = 0.0_dp    
    ALLOCATE( Ax(A % NumberOfRows) )

    ! 1) stiffness matrix reduction
    ! 2) mass matrix reduction
    DO MatrixNo = 1,2

      IF( MatrixNo == 1 ) THEN
        CALL Info('CraigBamptonSolver','Generating the reduced stiffness matrix',Level=5)
        A % Values => A % BulkValues
        FileName = 'ReducedStiff.dat'
      ELSE
        IF( .NOT. ASSOCIATED( A % BulkMassValues ) ) THEN
          CALL Warn('CraigBamptonSolver','Bulk mass values not present!')
          CYCLE
        END IF
        CALL Info('CraigBamptonSolver','Generating the reduced mass matrix',Level=5)
        A % Values => A % BulkMassValues
        FileName = 'ReducedMass.dat'
      END IF

      DO k=1,NoComponentModes
        IF( k <= NoEigenModes ) THEN
          x = EVar % EigenVectors(k,:)
        ELSE
          x = Evar % ConstraintModes(k-NoEigenModes,:)
        END IF 

        CALL MatrixVectorMultiply( A,x,Ax )
      
        DO l=1,NoComponentModes

          IF( l <= NoEigenModes ) THEN
            x = EVar % EigenVectors(l,:)
          ELSE
            x = Evar % ConstraintModes(l-NoEigenModes,:)
          END IF

          Ahat(l,k) = SUM( x * Ax ) 
        END DO
      END DO
      
      CALL SaveReductionMatrix()

      IF( ListGetLogical( Solver % Values,'Calculate Matrix Norm',Found ) ) THEN
        Norm = FrobeniusNorm( Ahat, NoComponentModes ) 
        IF( MatrixNo == 1 ) THEN
          Solver % Variable % Values = Norm 
          Solver % Variable % Norm = Norm
        END IF
        IF( MatrixNo == 1 ) THEN
          WRITE( Message,'(A,ES15.6)') 'Matrix norm for reduced stiffness matrix: ',Norm 
        ELSE
          WRITE( Message,'(A,ES15.6)') 'Matrix norm for reduced mass matrix: ',Norm 
        END IF
        CALL Info('CraigBamptonSolver',Message,Level=5)
      END IF

    END DO
  END IF


  A % Values => SaveValues
  DEALLOCATE( Ahat ) 
  
  CALL Info('CraigBamptonSolver','All done for now',Level=5)


CONTAINS 


  FUNCTION FrobeniusNorm( A, n ) RESULT ( Norm ) 
    REAL(KIND=dp) :: A(:,:)
    INTEGER :: n
    REAL(KIND=dp) :: Norm
    
    INTEGER :: i,j

    Norm = 0.0_dp
    DO i=1,n
      DO j=1,n
        Norm = Norm + A(i,j)**2
      END DO
    END DO
    Norm = SQRT( Norm ) 
    
  END FUNCTION FrobeniusNorm


  SUBROUTINE SaveReductionDofs()

    FileName = 'ReductionHeader.dat'
    CALL Info('CraigBamptonSolver','Saving information on the reduction basis to: '//TRIM(FileName),Level=7)
    OPEN (10, FILE=FileName )
    WRITE(10,'(I0)') NoEigenModes
    WRITE(10,'(I0)') NoConstraintModes
    CLOSE(10)

    FileName = 'ReductionIndeces.dat'
    CALL Info('CraigBamptonSolver','Saving Constraint Modes Indeces to: '//TRIM(FileName),Level=7)
    OPEN (10, FILE=FileName )
    DO i=1,A % NumberOfRows
      j = Evar % ConstraintModesIndeces(i)
      IF( j == 0 ) CYCLE
      WRITE(10,'(I0)') i 
    END DO
    CLOSE(10)

  END SUBROUTINE SaveReductionDofs


  SUBROUTINE SaveReductionBasis()

    REAL(KIND=dp) :: val

    CALL Info('CraigBamptonSolver','Saving the actual component modes to: '//TRIM(FileName),Level=7)
    OPEN (10, FILE=FileName)

    DO i=1,Mesh % NumberOfNodes
      j = Evar % Perm(i)
      IF( j == 0 ) CYCLE
      DO k=1,NoComponentModes
        WRITE(10,'(I0)',ADVANCE='NO') i
        DO l = 1, Dofs
          IF( k <= NoEigenModes ) THEN
            val = EVar % EigenVectors(k,Dofs*(j-1)+l)
          ELSE 
            val = EVar % ConstraintModes(k-NoEigenModes,Dofs*(j-1)+l)
          END IF
          WRITE(10,'(ES16.7)',ADVANCE='NO') val
        END DO
      END DO
      WRITE(10,'(A)') ' '
    END DO
    CLOSE( 10 ) 
  END SUBROUTINE SaveReductionBasis


  SUBROUTINE SaveReductionMatrix()
    CALL Info('CraigBamptonSolver','Saving the reduction matrix to: '//TRIM(FileName),Level=7)

    OPEN (10, FILE=FileName )
    DO k=1,NoComponentModes
      WRITE(10,*) Ahat(k,:)
    END DO
    CLOSE(10) 
  END SUBROUTINE SaveReductionMatrix

!------------------------------------------------------------------------------
END SUBROUTINE CraigBamptonSolver
!------------------------------------------------------------------------------
