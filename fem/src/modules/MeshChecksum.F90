!/*****************************************************************************
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
! *  A utility to compute and compare checksums for mesh input.
! *
! ******************************************************************************
! *
! *  Authors: Peter Råback
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 16.9.2020
! *
! *****************************************************************************/

SUBROUTINE MeshChecksum_init( Model,Solver,dt,Transient)

  USE DefUtils
  USE MeshUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  ! Introduce variable so that we have a place for pseudonorm
  CALL ListAddNewString( Solver % Values,'Variable',&
      '-nooutput -global meshcheck_var') 

END SUBROUTINE MeshChecksum_Init


SUBROUTINE MeshChecksum( Model,Solver,dt,Transient)

  USE DefUtils
  USE MeshUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!---------------------------------------------------------------
  LOGICAL :: Found
  INTEGER :: i, j, t, ind, tag, nsize
  TYPE(Element_t), POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  REAL(KIND=dp) :: CheckSum(8), c, PseudoNorm
  REAL(KIND=dp), POINTER :: WrkArray(:,:),RefSum(:)
  CHARACTER(*), PARAMETER :: Caller = 'MeshChecksum'

  
  CALL Info(Caller,'Checking for mesh consistency')
  
  Mesh => Solver % Mesh
  Params => Solver % Values

  CheckSum = 0.0
  CheckSum(1) = Mesh % NumberOfBulkElements + Mesh % NumberOfBulkElements + &
      Mesh % NumberOfNodes


  DO t=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements    
    Element => Mesh % Elements(t)
    IF( t <= Mesh % NumberOfBulkElements ) THEN
      ind = 2
      tag = Element % BodyId
    ELSE
      ind = 5
      tag = 0
      IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
        tag = Element % BoundaryInfo % Constraint
      END IF
    END IF
    
    CheckSum(ind) = CheckSum(ind) + Element % Type % ElementCode
    CheckSum(ind+1) = CheckSum(ind+1) + tag
    CheckSum(ind+2) = CheckSum(ind+2) + SUM( Element % NodeIndexes )
  END DO
  
  CheckSum(8) = SUM( Mesh % Nodes % x ) + &
      SUM( Mesh % Nodes % y) + SUM( Mesh % Nodes % z)
  nsize = 8
  
  PRINT *,'Checksums for file output:'
  PRINT *,'ThisResults:',NINT(CheckSum(1:7)),CheckSum(8)
  
  WrkArray => ListGetConstRealArray( Params,'Reference Values',Found )
  IF( Found ) THEN
    RefSum => WrkArray(:,1)
    PRINT *,'RefResults:',NINT(RefSum(1:7)),RefSum(8)    

    PseudoNorm = 0.0
    j = 0
    DO i=1,nsize
      IF( ABS(RefSum(i)) > EPSILON(c) ) THEN
        c = CheckSum(i) / RefSum(i)
        c = MAX( c, 1.0_dp /c ) 
      ELSE
        c = 1.0_dp + CheckSum(i) 
      END IF
      PseudoNorm = PseudoNorm + c / nsize
    END DO          
    PRINT *,'PseudoNorm:',PseudoNorm
    
    ! By construction the reference norm should now be one!
    Solver % Variable % Values = PseudoNorm
    CALL ListAddNewConstReal( Params,'Reference Norm',1.0_dp)
  END IF

END SUBROUTINE MeshChecksum
