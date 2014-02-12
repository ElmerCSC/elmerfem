!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
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
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!-----------------------------------------------------------------------------
!>  Utility routines for radiation computation
!-----------------------------------------------------------------------------

MODULE Radiation

   USE ElementUtils
   USE CoordinateSystems
   USE DefUtils

   IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
   FUNCTION ComputeRadiationLoad( Model, Mesh, Element, Temperature, &
                 Reorder, Emissivity, AngleFraction) RESULT(T)
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Model_t) :: Model
     TYPE(Element_t)  :: Element
     INTEGER :: Reorder(:)
     REAL(KIND=dp) :: Temperature(:),Emissivity
     REAL(KIND=dp), OPTIONAL :: AngleFraction
     REAL(KIND=dp) :: T

     REAL(KIND=dp) :: Asum
     TYPE(Element_t),POINTER  :: CurrentElement
     INTEGER :: i,j,n
     REAL(KIND=dp), POINTER :: Vals(:)
     INTEGER, POINTER :: Cols(:)
     REAL(KIND=dp) :: A1,A2,Emissivity1
     LOGICAL :: Found
!------------------------------------------------------------------------------

     A1 = Emissivity * ElementArea(Mesh,Element,Element % TYPE % NumberOfNodes)

     Cols => Element % BoundaryInfo % GebhardtFactors % Elements
     Vals => Element % BoundaryInfo % GebhardtFactors % Factors

     T = 0.0D0
     Asum = 0.0d0
     DO i=1,Element % BoundaryInfo % GebhardtFactors % NumberOfFactors
       CurrentElement => Mesh % Elements(Cols(i))
       n = CurrentElement % TYPE % NumberOfNodes

       Emissivity1 = SUM(ListGetReal( Model % BCs(CurrentElement % &
            BoundaryInfo % Constraint) % Values, 'Emissivity', &
            n, CurrentElement % NodeIndexes, Found) ) / n
       IF(.NOT. Found) THEN
          Emissivity1 = SUM(GetParentMatProp('Emissivity',CurrentElement)) / n 
       END IF

       A2 = Emissivity1 * ElementArea( Mesh,CurrentElement, &
              CurrentElement % TYPE % NumberOfNodes )

       T = T + A2 * ABS(Vals(i)) * &
         SUM(Temperature(Reorder(CurrentElement % NodeIndexes))/N)**4

       Asum = Asum + A2 * ABS(Vals(i))
     END DO

     T = (T/A1)**(1.0D0/4.0D0)
     IF(PRESENT(AngleFraction)) THEN
       AngleFraction = Asum / A1
     END IF

   END FUNCTION ComputeRadiationLoad
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
   FUNCTION ComputeRadiationCoeff( Model,Mesh,Element,k ) RESULT(T)
!------------------------------------------------------------------------------

     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Model_t)  :: Model
     TYPE(Element_t) :: Element
     INTEGER :: k
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: T

     TYPE(Element_t),POINTER  :: CurrentElement
     INTEGER :: i,j,n
     LOGICAL :: Found

     REAL(KIND=dp) :: Area,Emissivity
!------------------------------------------------------------------------------

     CurrentElement => Model % Elements( &
             Element % BoundaryInfo % GebhardtFactors % Elements(k) )
     n = CurrentElement % TYPE % NumberOfNodes

     Emissivity = SUM(ListGetReal(Model % BCs(CurrentElement % &
        BoundaryInfo % Constraint) % Values, 'Emissivity', &
        n, CurrentElement % NodeIndexes, Found)) / n
     IF(.NOT. Found) THEN
        Emissivity = SUM(GetParentMatProp('Emissivity',CurrentElement)) / n 
     END IF

     Area = Emissivity * ElementArea( Mesh,CurrentElement, n)

     T =  ABS(Element % BoundaryInfo % GebhardtFactors % Factors(k)) * Area

   END FUNCTION ComputeRadiationCoeff
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END MODULE Radiation
!------------------------------------------------------------------------------

!> \}