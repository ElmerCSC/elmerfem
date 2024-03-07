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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Oct 1996
! *
! ****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!-----------------------------------------------------------------------------
!>  Module defining coordinate systems.
!-----------------------------------------------------------------------------

MODULE CoordinateSystems

  USE Types

  IMPLICIT NONE

  INTEGER, PARAMETER :: Cartesian = 1
  INTEGER, PARAMETER :: Cylindric = 2, CylindricSymmetric = 3, AxisSymmetric = 4
  INTEGER, PARAMETER :: Polar = 5
  INTEGER :: Coordinates = Cartesian

!----------------------------------------------------------------------------
CONTAINS
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  FUNCTION CylindricalMetric( r,z,t ) RESULT(metric)

    REAL(KIND=dp) :: r,z,t
    REAL(KIND=dp), DIMENSION(3,3) :: Metric
 
    Metric = 0.0d0
    Metric(1,1) = 1.0d0
    Metric(2,2) = 1.0d0
    Metric(3,3) = 1.0d0
    IF ( r /= 0.0d0 ) Metric(3,3) = 1.0d0 / (r**2)
  END FUNCTION CylindricalMetric
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  FUNCTION CylindricalSqrtMetric( r,z,t ) RESULT(s)

    REAL(KIND=dp) :: r,z,t,s

    s = r
  END FUNCTION CylindricalSqrtMetric
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  FUNCTION CylindricalSymbols( r,z,t ) RESULT(symbols)
    REAL(KIND=dp) :: r,z,t

    REAL(KIND=dp), DIMENSION(3,3,3) :: symbols

    Symbols = 0.0d0
    Symbols(3,3,1) = -r

    IF ( r /= 0.0d0 ) THEN
       Symbols(1,3,3) = 1.0d0 / r
       Symbols(3,1,3) = 1.0d0 / r
    END IF
  END FUNCTION CylindricalSymbols
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  FUNCTION CylindricalDerivSymbols( r,z,t ) RESULT(dsymbols)
    REAL(KIND=dp) :: r,z,t

    REAL(KIND=dp), DIMENSION(3,3,3,3) :: dsymbols

    dSymbols = 0.0d0
    dSymbols(3,3,1,1) = -1.0d0

    IF ( r /= 0.0 ) THEN
       dSymbols(1,3,3,1) = -1.0d0 / (r**2)
       dSymbols(3,1,3,1) = -1.0d0 / (r**2)
    END IF
  END FUNCTION CylindricalDerivSymbols
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  FUNCTION PolarMetric( r,p,t ) RESULT(Metric)
    REAL(KIND=dp) :: r,p,t
    INTEGER :: i
    REAL(KIND=dp), DIMENSION(3,3) :: Metric

    Metric = 0.0d0
    DO i=1,3
       Metric(i,i) = 1.0d0
    END DO

    IF ( r /= 0.0d0 ) THEN
       Metric(2,2) = 1.0d0 / (r**2 * COS(t)**2)
       IF ( CoordinateSystemDimension() == 3 ) THEN
          Metric(3,3) = 1.0d0 / r**2
       END IF
    END IF
  END FUNCTION PolarMetric
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  FUNCTION PolarSqrtMetric( r,p,t ) RESULT(s)
    REAL(KIND=dp) :: r,p,t,s

    IF ( CoordinateSystemDimension() == 2 ) THEN
       s = SQRT( r**2 * COS(t)**2 )
    ELSE
       s = SQRT( r**4 * COS(t)**2 )
    END IF
  END FUNCTION PolarSqrtMetric
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  FUNCTION PolarSymbols( r,p,t ) RESULT(symbols)
    REAL(KIND=dp) :: r,p,t

    REAL(KIND=dp), DIMENSION(3,3,3) :: symbols

    Symbols = 0.0d0        
    Symbols(2,2,1) = -r * COS(t)**2
    IF ( r /= 0.0d0 ) THEN
       Symbols(1,2,2) = 1.0d0 / r
       Symbols(2,1,2) = 1.0d0 / r
    END IF

    IF ( CoordinateSystemDimension() == 3 ) THEN
       Symbols(3,3,1) = -r
       Symbols(2,2,3) = SIN(t)*COS(t)

       Symbols(2,3,2) = -TAN(t)
       Symbols(3,2,2) = -TAN(t)

       IF ( r /= 0.0d0 ) THEN
          Symbols(3,1,3) = 1.0d0 / r
          Symbols(1,3,3) = 1.0d0 / r
       END IF
    END IF
  END FUNCTION PolarSymbols
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  FUNCTION PolarDerivSymbols( r,p,t ) RESULT(dsymbols)
    REAL(KIND=dp) :: r,p,t

    REAL(KIND=dp), DIMENSION(3,3,3,3) :: dSymbols

    dSymbols = 0.0d0
    dSymbols(2,2,1,1) = -COS(t)**2
    IF ( r /= 0.0d0 ) THEN
       dSymbols(1,2,2,1) = -1.0d0 / r**2
       dSymbols(2,1,2,1) = -1.0d0 / r**2
    END IF

    IF ( CoordinateSystemDimension() == 3 ) THEN
       dSymbols(2,2,1,3) = -2*r*SIN(t)*COS(t)
       dSymbols(3,3,1,1) = -1
       dSymbols(2,2,3,3) =  COS(t)**2 - SIN(t)**2

       dSymbols(2,3,2,3) = -1.0d0 / COS(t)**2
       dSymbols(3,2,2,3) = -1.0d0 / COS(t)**2
       
       IF ( r /= 0.0d0 ) THEN
          dSymbols(1,3,3,1) = -1.0d0 / r**2
          dSymbols(3,1,3,1) = -1.0d0 / r**2
       END IF
    END IF
  END FUNCTION PolarDerivSymbols
!----------------------------------------------------------------------------



!----------------------------------------------------------------------------
  FUNCTION CoordinateSqrtMetric( X,Y,Z ) RESULT( SqrtMetric )
    REAL(KIND=dp) :: X,Y,Z,SqrtMetric

    IF ( Coordinates == Cartesian ) THEN
       SqrtMetric = 1.0d0 
    ELSE IF ( Coordinates >= Cylindric .AND. &
       Coordinates <= AxisSymmetric ) THEN
       SqrtMetric = CylindricalSqrtMetric( X,Y,Z )
    ELSE IF ( Coordinates == Polar ) THEN
       SqrtMetric = PolarSqrtMetric( X,Y,Z )
    END IF

  END FUNCTION CoordinateSqrtMetric
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  FUNCTION CurrentCoordinateSystem() RESULT(Coords)
    INTEGER :: Coords

    Coords = Coordinates
  END FUNCTION CurrentCoordinateSystem
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  SUBROUTINE CoordinateSystemInfo( Metric,SqrtMetric, &
              Symbols,dSymbols,X,Y,Z )

    REAL(KIND=dp) :: Metric(3,3),SqrtMetric, &
        Symbols(3,3,3),dSymbols(3,3,3,3)

    INTEGER :: i
    REAL(KIND=dp) :: X,Y,Z

    IF ( Coordinates == Cartesian ) THEN

       Metric  = 0.0d0
       DO i=1,3
          Metric(i,i) = 1.0d0
       END DO

       SqrtMetric = 1.0d0 
       Symbols    = 0.0d0
       dSymbols   = 0.0d0

    ELSE IF ( Coordinates >= Cylindric .AND.  &
         Coordinates <= AxisSymmetric ) THEN

       SqrtMetric = CylindricalSqrtMetric( X,Y,Z )
       Metric     = CylindricalMetric( X,Y,Z )
       Symbols    = CylindricalSymbols( X,Y,Z )
       dSymbols   = CylindricalDerivSymbols( X,Y,Z )

    ELSE IF ( Coordinates == Polar ) THEN

       SqrtMetric = PolarSqrtMetric( X,Y,Z )
       Metric     = PolarMetric( X,Y,Z )
       Symbols    = PolarSymbols( X,Y,Z )
       dSymbols   = PolarDerivSymbols( X,Y,Z )

    END IF

  END SUBROUTINE CoordinateSystemInfo
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
  FUNCTION CoordinateSystemDimension() RESULT( dim )
    INTEGER :: dim ! ,csys

!   csys = CurrentCoordinateSystem()
    dim  = CurrentModel % Dimension
  END FUNCTION CoordinateSystemDimension
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
END MODULE CoordinateSystems
!----------------------------------------------------------------------------

!> \} ElmerLib
