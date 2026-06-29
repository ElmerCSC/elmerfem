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
! * ELMER/FEM Viewfactor computation
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \ingroup Programs
!> \{

!> \defgroup ViewFactors Program RadiatorFactors
!> \{

!------------------------------------------------------------------------------
!> A separate program that computes the radiatiave source coefficients to an
!> external file. 
!>
!> This file is later used within the ElmerSolver. If the radiative files
!> do not exist, a system call for this program is performed. 
!------------------------------------------------------------------------------
      
  PROGRAM RadiatorFactors
    INTEGER :: i,j,k
    CHARACTER(LEN=256) :: s, t

     CALL GET_COMMAND_ARGUMENT(0,t)
     i = COMMAND_ARGUMENT_COUNT()
     s = ''
     IF (i>0) CALL GET_COMMAND_ARGUMENT(i,s)

     i = LEN_TRIM(t)
     j = LEN_TRIM(s)
     DO k=i,1,-1
       IF ( t(k:k) == '/' .OR. t(k:k) == '\' .OR. t(k:k) == ':' ) EXIT 
     END DO
     k = k + 1
     t(k:k+10) = 'ViewFactors'
     t(k+11:) = ''

     i = LEN_TRIM(t)
     t = t(1:i)//' -radiators '//s(1:j)
     print*,t

     CALL system(t)
  END PROGRAM RadiatorFactors

  
!> \}
