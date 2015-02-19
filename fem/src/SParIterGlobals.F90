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

  TYPE (ParEnv_t), SAVE, TARGET :: ParEnv
  TYPE (SParIterSolverGlobalD_t), POINTER :: PIGpntr
  TYPE (SParIterSolverGlobalD_t), POINTER :: GlobalData
END MODULE SParIterGlobals

!> \}