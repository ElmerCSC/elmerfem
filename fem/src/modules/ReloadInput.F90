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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 07 Oct 2002
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!>  Reaload input file from disk to get some dynamic control over solver
!>  parameters. This should be used carefully since not all information may 
!> be read without introducing problems. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ReloadInput( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE Types
  USE Lists
  USE ModelDescription

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

   CHARACTER(LEN=MAX_NAME_LEN) :: ModelName, MeshDir, MeshName
   INTEGER :: IOUnit
   
!------------------------------------------------------------------------------

   OPEN(NEWUNIT=IOUnit,file='ELMERSOLVER_REREADINFO', STATUS='OLD', ERR=10 )
   READ(IOUnit,'(a)') ModelName
   CLOSE(IOUnit)

   OPEN( InFileUnit, File=ModelName )
   CALL LoadInputFile( Model, InFileUnit, ModelName, MeshDir, MeshName, .FALSE.,.FALSE. )
   CLOSE( InFileUnit )
 
   RETURN

10 CONTINUE

   OPEN(NEWUNIT=IOUnit,file='ELMERSOLVER_REREADINFO', STATUS='OLD', ERR=20 )
   READ(IOUnit,'(a)') ModelName
   CLOSE(IOUnit)

20 CONTINUE

   RETURN

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ReloadInput
!------------------------------------------------------------------------------
