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
  REAL(KIND=dp) :: x1, x0, x, w, f
  INTEGER :: i,j,n,NoPar
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, GotIt

  IF( ParEnv % PEs > 1 ) THEN
    IF( ParEnv % MyPE > 0 ) RETURN
  END IF

  CALL Info('SaveDependence','Saving dependencies in a table')

  Params => GetSolverParams()

  FileName = ListGetString( Params,'Filename',Found)
  IF(.NOT. Found ) FileName = 'dep.dat'


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

  OPEN( 10, FILE=FileName )

  DO i=1,n
    w = (1.0_dp*(i-1))/(n-1)
    x = x0 + w*(x1-x0)

    WRITE (10,'(I6,ES15.6)',ADVANCE='NO') i,x
    
    DO j=1,NoPar
      WRITE (ParName,'(A,I0)') 'Expression ',j
      f = ListGetFun( Params,ParName,x )
      WRITE (10,'(ES15.6)',ADVANCE='NO') f     
    END DO

    WRITE (10,'(A)') ' '     
  END DO
  
  CLOSE( 10 ) 

END SUBROUTINE SaveDependence
!------------------------------------------------------------------------------

!> \}
