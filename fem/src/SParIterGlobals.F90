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
!/******************************************************************************
! *
! *  Authors: Jouni Malinen, Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2000
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!>  This module contains global variables (or pointers to them)
!>  needed by the parallel version of the ELMER iterative solver.
!------------------------------------------------------------------------------
!> \ingroup ElmerLib
!> \{


MODULE SParIterGlobals

  USE Types

  IMPLICIT NONE

real(kind=dp):: xxx, yyy

  TYPE HUTICtlT
     INTEGER :: Method
     INTEGER :: Precond
     DOUBLE PRECISION :: Tolerance
     INTEGER :: MaxIter
     INTEGER :: DebugLevel
  END TYPE HUTICtlT


  TYPE ErrInfoT
     INTEGER :: HUTIStatus
  END TYPE ErrInfoT

  ! Following is in correct place

  TYPE (ParEnv_t), SAVE, TARGET :: ParEnv_Common
  TYPE (ParEnv_t), POINTER, SAVE :: ParEnv => ParEnv_Common
  TYPE (SParIterSolverGlobalD_t), POINTER :: PIGpntr
  TYPE (SParIterSolverGlobalD_t), POINTER :: GlobalData

CONTAINS

!------------------------------------------------------------------------------
!> Aim the global > ParEnv < at the parallel environment of the given matrix.
!> The active partitions and the neighbours describe a matrix, not a solver: a
!> single solver may own several matrices with a differing parallel structure,
!> for example the levels of a p/algebraic multigrid preconditioner, the blocks
!> of a block system, or a collection matrix. The per matrix environment is
!> therefore owned by > Matrix % ParMatrix % ParEnv <, and > Solver % ParEnv <
!> holds a mirror of the matrix that is currently being worked on. The mirror is
!> what > ParEnv < is aimed at, so that the target stays valid also when the
!> matrix of a solver is replaced or freed.
!------------------------------------------------------------------------------
  SUBROUTINE SetMatrixParEnv( Matrix, ParMatrix )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: Matrix
    TYPE(SParIterSolverGlobalD_t), OPTIONAL, POINTER :: ParMatrix
!------------------------------------------------------------------------------
    TYPE(SParIterSolverGlobalD_t), POINTER :: p
!------------------------------------------------------------------------------
    p => Matrix % ParMatrix
    IF( PRESENT( ParMatrix ) ) p => ParMatrix

    IF( ASSOCIATED( p ) ) THEN
      ParEnv => p % ParEnv
      IF( ASSOCIATED( Matrix % Solver ) ) THEN
        Matrix % Solver % ParEnv = p % ParEnv
        ParEnv => Matrix % Solver % ParEnv
      END IF
    ELSE IF( ASSOCIATED( Matrix % Solver ) ) THEN
      ParEnv => Matrix % Solver % ParEnv
    END IF

    ParEnv % ActiveComm = Matrix % Comm
!------------------------------------------------------------------------------
  END SUBROUTINE SetMatrixParEnv
!------------------------------------------------------------------------------

END MODULE SParIterGlobals

!> \}
