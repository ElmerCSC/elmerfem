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
                 Reorder, Emissivity, AngleFraction, Areas, Emiss ) RESULT(T)
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Model_t) :: Model
     TYPE(Element_t)  :: Element
     INTEGER :: Reorder(:)
     REAL(KIND=dp), OPTIONAL :: AngleFraction, Areas(:), Emiss(:)
     REAL(KIND=dp) :: T
     REAL(KIND=dp) :: Temperature(:), Emissivity

     REAL(KIND=dp) :: Asum
     TYPE(Element_t),POINTER  :: RadElement
     INTEGER :: i,j,n, bindex,nf 
     REAL(KIND=dp), POINTER :: Vals(:)
     INTEGER, POINTER :: Cols(:)
     REAL(KIND=dp) :: A1,A2,Emissivity1
     LOGICAL :: Found
!------------------------------------------------------------------------------

     IF( .NOT. ASSOCIATED( Element % BoundaryInfo % RadiationFactors ) ) THEN
       CALL Fatal('ComputeRadiationLoad','Gebhart factors not calculated for boundary!')
     END IF
     
     nf = Element % BoundaryInfo % RadiationFactors % NumberOfFactors          
     
     IF(PRESENT(Areas) .AND. PRESENT(Emiss)) THEN

       bindex = Element % ElementIndex - Mesh % NumberOfBulkElements
       A1  = Emiss(bIndex)

       Cols => Element % BoundaryInfo % RadiationFactors % Elements
       Vals => Element % BoundaryInfo % RadiationFactors % Factors

       T = 0._dp
       Asum = 0._dp
       DO i=1,nf
         RadElement => Mesh % Elements(Cols(i))
         n = RadElement % TYPE % NumberOfNodes
         T=T+Vals(i)*SUM(Temperature(Reorder(RadElement % NodeIndexes))/n)**4 / A1
         Asum = Asum + Vals(i)
       END DO
     ELSE
       A1 = Emissivity

       Cols => Element % BoundaryInfo % RadiationFactors % Elements
       Vals => Element % BoundaryInfo % RadiationFactors % Factors

       T = 0.0_dp
       Asum = 0.0_dp
       DO i=1,nf
         RadElement => Mesh % Elements(Cols(i))
         n = RadElement % TYPE % NumberOfNodes
         T = T + Vals(i) * &
           SUM(Temperature(Reorder(RadElement % NodeIndexes))/n)**4 / A1
         Asum = Asum + Vals(i)
       END DO
     END IF

     T = (T)**(1._dp/4._dp)

     IF(PRESENT(AngleFraction)) AngleFraction = Asum
!------------------------------------------------------------------------------
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
             Element % BoundaryInfo % RadiationFactors % Elements(k) )
     n = CurrentElement % TYPE % NumberOfNodes

     Emissivity = SUM(ListGetReal(Model % BCs(CurrentElement % &
        BoundaryInfo % Constraint) % Values, 'Emissivity', &
        n, CurrentElement % NodeIndexes, Found)) / n
     IF(.NOT. Found) THEN
        Emissivity = SUM(GetParentMatProp('Emissivity',CurrentElement)) / n 
     END IF

     Area = Emissivity * ElementArea( Mesh,CurrentElement, n)

     T =  ABS(Element % BoundaryInfo % RadiationFactors % Factors(k)) * Area

   END FUNCTION ComputeRadiationCoeff
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END MODULE Radiation
!------------------------------------------------------------------------------

!> \}
