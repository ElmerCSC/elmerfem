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
! *  Subroutine for saving 1D dependence into files.
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

SUBROUTINE SaveDependence_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  INTEGER :: NormInd
  LOGICAL :: GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  ! If we want to show a pseudonorm add a variable for which the norm
  ! is associated with.
  NormInd = ListGetInteger( Solver % Values,'Show Norm Index',GotIt)
  IF( NormInd > 0 ) THEN
    SolverName = ListGetString( Solver % Values, 'Equation',GotIt)
    IF( .NOT. ListCheckPresent( Solver % Values,'Variable') ) THEN
      CALL ListAddString( Solver % Values,'Variable',&
          '-nooutput -global '//TRIM(SolverName)//'_var')
    END IF
  END IF

END SUBROUTINE SaveDependence_init


!------------------------------------------------------------------------------
!> This subroutine saves 1D dependence of a given property.
!------------------------------------------------------------------------------
SUBROUTINE SaveDependence( Model,Solver,dt,TransientSimulation )
  
  USE Types
  USE Lists
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
  CHARACTER(LEN=MAX_NAME_LEN) :: FileName, ParName, OutputDirectory
  REAL(KIND=dp) :: x1, x0, x, w, f, Norm
  INTEGER :: i,j,n,NoPar,NormInd,IOUnit
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, GotIt
  
  IF( ParEnv % PEs > 1 ) THEN
    IF( ParEnv % MyPE > 0 ) RETURN
  END IF

  CALL Info('SaveDependence','--------------------------------')
  CALL Info('SaveDependence','Saving dependencies in a table')

  Params => GetSolverParams()

  FileName = ListGetString( Params,'Filename',Found)
  IF(.NOT. Found ) FileName = 'dep.dat'


  !------------------------------------------------------------------------------
  ! For consistency checks one compute a pseudonorm.
  !------------------------------------------------------------------------------
  NormInd = ListGetInteger( Params,'Show Norm Index',GotIt)
  Norm = 0.0_dp

  IF ( .NOT. FileNameQualified(FileName) ) THEN
    OutputDirectory = GetString( Params,'Output Directory',GotIt)
    IF( GotIt .AND. LEN_TRIM(OutputDirectory) > 0 ) THEN
      FileName = TRIM(OutputDirectory)// '/' //TRIM(Filename)
      CALL MakeDirectory( TRIM(OutputDirectory) // CHAR(0) )
    ELSE IF( LEN_TRIM(OutputPath ) > 0 ) THEN
      Filename = TRIM(OutputPath)// '/' //TRIM(Filename)
    END IF
  END IF

  IF( GetLogical(Params,'Filename Numbering',GotIt)) THEN
    Filename = NextFreeFilename( Filename )
  END IF
  

  n = ListGetInteger( Params,'Number of points',minv=2)
  x0 = ListGetCReal( Params,'Lower limit')
  x1 = ListGetCReal( Params,'Upper Limit')

  NoPar = 0
  DO j=1,100
    WRITE (ParName,'(A,I0)') 'Expression ',j
    IF( ListCheckPresent( Params, ParName ) ) THEN
      NoPar = j
    ELSE
      EXIT
    END IF
  END DO

  IF( NoPar == 0 ) THEN
    CALL Warn('SaveDependence','No parameter given!')
    RETURN
  END IF

  OPEN(NEWUNIT=IOUnit, FILE=FileName )

  DO i=1,n
    w = (1.0_dp*(i-1))/(n-1)
    x = x0 + w*(x1-x0)

    WRITE (IOUnit,'(I6,ES15.6)',ADVANCE='NO') i,x
    
    DO j=1,NoPar
      WRITE (ParName,'(A,I0)') 'Expression ',j
      f = ListGetFun( Params,ParName,x )
      WRITE (IOUnit,'(ES15.6)',ADVANCE='NO') f     

      IF( NormInd == j ) Norm = Norm + f*f
    END DO

    WRITE (IOUnit,'(A)') ' '     
  END DO
  
  CLOSE( IOUnit ) 

  IF( NormInd > 0 ) THEN
    Norm = SQRT( Norm / n )
    Solver % Variable % Values = Norm
    Solver % Variable % Norm = Norm
  END IF
  
END SUBROUTINE SaveDependence
!------------------------------------------------------------------------------

!> \}
