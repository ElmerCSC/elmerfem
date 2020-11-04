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
! *  Authors: Peter Raback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2.1.2011
! *
! *****************************************************************************/

!----------------------------------------------------------------------------------
!> Module for 1D electrostatics boundary conditions. 
!> Computes electrostatic force, energy, charge and spring density on the boundary
!> assuming a 1-dimensional electrostatic model. The model may be complicated 
!> by holes or thin dielectric layer. 
!----------------------------------------------------------------------------------

MODULE StatElecBoundaryUtils

  USE DefUtils
  IMPLICIT NONE
  
CONTAINS
  
  FUNCTION StatElecBoundaryGeneric( Model, NodeNumber, dummy, Mode ) RESULT( value )
    TYPE(Model_t) :: Model
    INTEGER :: NodeNumber, Mode
    REAL(KIND=dp) :: dummy, value
    
    TYPE(Element_t), POINTER :: CurrentElement => NULL()
    INTEGER :: i, n
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: AllocationsDone, HoleCorrection, GotIt, LayerExists
    REAL(KIND=dp) :: PermittivityOfVacuum, EffectiveAperture, Alpha, Beta, Gamma
    REAL(KIND=dp), ALLOCATABLE :: PotentialDifference(:), Permittivity(:),&
        LayerThickness(:), LayerPermittivity(:), Thickness(:), &
        HoleSize(:), HoleFraction(:),ElemAperture(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: HoleType
    TYPE(ValueList_t), POINTER :: List
    

    SAVE :: AllocationsDone, n, NodeIndexes, List, CurrentElement, &
        PotentialDifference, Permittivity, LayerExists, Thickness, &
        PermittivityOfVacuum, LayerThickness, LayerPermittivity, &
        HoleSize, HoleFraction, HoleCorrection, ElemAperture

    IF(.NOT. AllocationsDone ) THEN      
      n = Model % Mesh % MaxElementNodes
      ALLOCATE( PotentialDifference(n), Permittivity(n), &
          LayerThickness(n), LayerPermittivity(n), Thickness(n), &
          HoleSize(n), HoleFraction(n), ElemAperture(n))
      AllocationsDone = .TRUE.
    END IF
    
    
    IF( .NOT. ASSOCIATED( CurrentElement, Model % CurrentElement )) THEN
      CurrentElement => Model % CurrentElement
      
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes

      IF( CurrentElement % ElementIndex > Model % Mesh % NumberOfBulkElements ) THEN
        List => GetBC(CurrentElement)
      ELSE
        List => GetMaterial(CurrentElement)
      END IF

      PotentialDifference(1:n) = ListGetReal( List, &
          'Potential Difference',n,NodeIndexes,gotIt )
      
      PermittivityOfVacuum = ListGetConstReal( Model % Constants, &
        'Permittivity of Vacuum')
      
      Permittivity(1:n) = ListGetReal( List, &
        'Permittivity',n,NodeIndexes,GotIt)
      IF(.NOT. GotIt) Permittivity(1:n) = ListGetReal( List, &
        'Relative Permittivity',n,NodeIndexes)
      
      Permittivity(1:n) = Permittivity(1:n) * PermittivityOfVacuum
      
      LayerThickness(1:n) = ListGetReal( List, &
          'Layer Thickness',n,NodeIndexes,LayerExists)
      
      IF(LayerExists) THEN
        LayerPermittivity(1:n) = ListGetReal( List, &
            'Layer Permittivity',n,NodeIndexes)
        IF(.NOT.GotIt) LayerPermittivity = 1.0d0 
        LayerPermittivity(1:n) = LayerPermittivity(1:n) * PermittivityOfVacuum
      END IF
      
      HoleType = ListGetString(List,'Hole Type',HoleCorrection)
      IF( HoleCorrection ) THEN
        HoleSize(1:n) = ListGetReal(List,'Hole Size', n, NodeIndexes)
        Thickness(1:n) = ListGetReal(List,'Hole Depth', n, NodeIndexes)
        HoleFraction(1:n) = ListGetReal( List,'Hole Fraction',n,NodeIndexes)
      END IF
      
      ElemAperture(1:n) = ListGetReal(List,'Gap Height',n,NodeIndexes)
    END IF
    
    
    DO i=1, n
      IF( NodeIndexes(i) == NodeNumber ) EXIT
    END DO
    
    
    IF(LayerExists) THEN
      IF(ElemAperture(i) > LayerThickness(i)) THEN
        EffectiveAperture = (ElemAperture(i) - LayerThickness(i)) + &
            LayerThickness(i)*Permittivity(i)/LayerPermittivity(i)
      ELSE 
        EffectiveAperture = ElemAperture(i) * Permittivity(i) / LayerPermittivity(i)
      END IF
    ELSE
      EffectiveAperture = ElemAperture(i) 
    END IF
    
    IF(HoleCorrection) THEN
      CALL ElectrostaticHoleCorrection(HoleType, HoleSize(i), Thickness(i), &
          HoleFraction(i), EffectiveAperture, Alpha, Beta, Gamma)
    END IF
    
    
    SELECT CASE ( Mode )
      
    CASE( 1 ) ! energy
      value = 0.5d0 * (PotentialDifference(i) ** 2.0d0) * &
          Permittivity(i) / EffectiveAperture            
      IF( HoleCorrection ) value = Alpha * value

      
    CASE ( 2 ) ! force
      value = -0.5d0 * (PotentialDifference(i) ** 2.0d0) * &
          Permittivity(i) /  EffectiveAperture**2.0            
      IF( HoleCorrection ) value = Beta * value      

    CASE ( 3 ) ! charge
      value = PotentialDifference(i) * &
          Permittivity(i) /  EffectiveAperture**2.0            
      IF( HoleCorrection ) value = Beta * value      
      
    CASE ( 4 ) ! spring
      value = (PotentialDifference(i) ** 2.0d0) * &
          Permittivity(i) /  EffectiveAperture**3.0
      IF( HoleCorrection ) value = Gamma * value      
      
    CASE DEFAULT
      CALL Fatal('StatElecBoundaryGeneric','Unknown mode for electrostatic boundary')
      
      
    END SELECT
  

  CONTAINS

    !------------------------------------------------------------------------------
    !> Computes the hole correction for the 1D electrostatic equation 
    !> from experimentally fitted numerical results.
    !------------------------------------------------------------------------------
    SUBROUTINE ElectrostaticHoleCorrection(holemodel, r, b, p, d, alpha, beta, gamma)
      ! r=radius, b=hole length, d=aperture, p=hole fraction
      CHARACTER(LEN=*) :: holemodel
      REAL(KIND=dp) :: r,d,b,p,alpha,beta,gamma,a,da,dda,c1,c2,dom
      
      SELECT CASE(holemodel)
        
      CASE ('slot')
        c1 = 2.3198
        c2 = 0.2284 
        
      CASE ('round')
        c1 = 4.2523d0
        c2 = 0.4133d0
        
      CASE ('square')
        c1 = 3.8434 
        c2 = 0.3148
        
      CASE DEFAULT 
        alpha = 1.0
        beta = 1.0
        gamma = 1.0
        
        CALL WARN('ComputeHoleCorrection','Unknown hole type')       
        
        RETURN
      END SELECT
      
      dom = 1.0d0 + c1*(d/r) + c2* (d/r)**2.0
      a = 1.0 - p * 1.0d0/dom
      da = p * (c1+2.0*c2*(d/r)) / dom**2.0
      dda = p * 2.0 * (c2-2.0* c1**2.0-3.0*c1*c2*(d/r)-3.0* c2**2.0 * (d/r)**2.0) / dom**3.0     
      
      alpha = a
      beta = a - da*(d/r)
      gamma = a - da*(d/r) + 0.5d0*dda*((d/r)**2.0)
      
    END SUBROUTINE ElectrostaticHoleCorrection
    !------------------------------------------------------------------------------
        
    
  END FUNCTION StatElecBoundaryGeneric
  

END MODULE StatElecBoundaryUtils




!------------------------------------------------------------------------------
!> Computes the electrostatic energy density on boundary using 1D model.
!> \ingroup UDF
!------------------------------------------------------------------------------
FUNCTION StatElecBoundaryEnergy( Model, NodeNumber, Gap ) RESULT( value )
  USE Types
  USE StatElecBoundaryUtils
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber, Mode
  REAL(KIND=dp) :: Gap, Value

  Mode = 1
  value = StatElecBoundaryGeneric( Model, NodeNumber, Gap, Mode )

END FUNCTION StatElecBoundaryEnergy


!------------------------------------------------------------------------------
!> Computes the electrostatic force density on boundary using 1D model.
!> \ingroup UDF
!------------------------------------------------------------------------------
FUNCTION StatElecBoundaryForce( Model, NodeNumber, Gap ) RESULT( value )
  USE Types
  USE StatElecBoundaryUtils
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber, Mode
  REAL(KIND=dp) :: Gap, Value

  Mode = 2
  value = StatElecBoundaryGeneric( Model, NodeNumber, Gap, Mode )
END FUNCTION StatElecBoundaryForce


!------------------------------------------------------------------------------
!> Computes the electrostatic charge density on boundary using 1D model.
!> \ingroup UDF
!------------------------------------------------------------------------------
FUNCTION StatElecBoundaryCharge( Model, NodeNumber, Gap ) RESULT( value )
  USE Types
  USE StatElecBoundaryUtils
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber, Mode
  REAL(KIND=dp) :: Gap, Value

  Mode = 3
  value = StatElecBoundaryGeneric( Model, NodeNumber, Gap, Mode )
END FUNCTION StatElecBoundaryCharge


!------------------------------------------------------------------------------
!> Computes the electrostatic spring density on boundary using 1D model.
!> \ingroup UDF
!------------------------------------------------------------------------------
FUNCTION StatElecBoundarySpring( Model, NodeNumber, Gap ) RESULT( value )
  USE Types
  USE StatElecBoundaryUtils
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber, Mode
  REAL(KIND=dp) :: Gap, Value

  Mode = 4
  value = StatElecBoundaryGeneric( Model, NodeNumber, Gap, Mode )
END FUNCTION StatElecBoundarySpring

