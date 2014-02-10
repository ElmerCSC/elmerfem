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
! *  This file contains subroutines for the wall boundary conditions of 
! *  the k-epsilon turbulence model on walls. 
! *
! ******************************************************************************
! *
! *  Authors: Jari Hämäläinen
! *  Address: VTT Energy
! *           P.O.Box 1603
! *           40101 Jyväskylä, Finland 
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 19 Jun 1996
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Solve the friction velocity of the previous iteration based 
!> on the wall law. 
!
!         Input:
!             DENSIT - Density
!             VISCOS - Viscosity 
!             DIST   - Distance from the wall
!             UT     - Tangential velocity of the previous iteration
!
!         Output:
!             UFRIC  - Friction velocity
!             DFX    - Derivative of the wall law
!
!------------------------------------------------------------------------------
      SUBROUTINE SOLVE_UFRIC(DENSIT,VISCOS,DIST,ROUGH,UT,UFRIC,DFX)

      IMPLICIT NONE
      DOUBLE PRECISION DENSIT,VISCOS,DIST,ROUGH,UT,UFRIC,DFX,TAUW,  &
      YPLUS, FX, WALL_LAW, D_WALL_LAW

      INTEGER :: ITER 
      INTEGER :: MAXITER=100
      DOUBLE PRECISION ::  TOL=1.0D-14
 
! Default value:
      TAUW = UT / DIST
      UFRIC = DSQRT( TAUW / DENSIT )

      DO ITER=1,MAXITER
         FX  = WALL_LAW( UFRIC,UT,DENSIT,VISCOS,DIST,ROUGH )
         DFX = D_WALL_LAW( UFRIC,UT,DENSIT,VISCOS,DIST,ROUGH )

! Newton step:
         IF (DFX.EQ.0.0d0) STOP 'dfx=0'
         UFRIC = UFRIC - FX/DFX
         YPLUS = DENSIT * UFRIC * DIST / VISCOS
         IF ( DABS(FX) <= TOL ) EXIT
      END DO

      IF ( DABS(FX) > 1.0d-9 ) WRITE(*,*)'Problems in SOLVE_UFRIC, FX=',FX

      RETURN
      END
      


!----------------------------------------------------------------------------
!> Give difference between the tangential velocity given by 
!> Reichardt´s wall law and the tangential velocity of the previous 
!> iteration.
!
!         Input:
!             UFRIC  - Friction velocity
!             UT     - Tangential velocity
!             DENSIT - Density
!             VISCOS - Viscosity
!             DIST   - Distance
!
!         Output:
!
!----------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION WALL_LAW(UFRIC,UT,DENSIT, &
                 VISCOS,DIST,ROUGH)

      IMPLICIT NONE
      DOUBLE PRECISION DENSIT,VISCOS,DIST,ROUGH,UT,UFRIC,DFX, &
      YPLUS

      DOUBLE PRECISION :: DKAPPA = 0.41D0, RAJA


      YPLUS = DENSIT*UFRIC*DIST / VISCOS

! Log-law:
!      RAJA=11.2658567D0 
!      IF(YPLUS.GE.RAJA) THEN
!         WALL_LAW=(UFRIC/DKAPPA)*DLOG(ROUGH*YPLUS)-UT
!      ELSE
!         WALL_LAW=UFRIC*YPLUS-UT
!      ENDIF

! Reichardts law:
      WALL_LAW=(UFRIC/DKAPPA)*DLOG(1.0D0+0.4D0*YPLUS)  &
            + UFRIC*7.8D0  &
              *(           &
                 1.0D0-DEXP(-YPLUS/11.0D0) &
                 -(YPLUS/11.0D0)*DEXP(-0.33D0*YPLUS) &
               ) - UT

      RETURN
      END


!----------------------------------------------------------------------------
!> Calculate derivative of the wall law.
!
!         Input:
!             UFRIC  - Friction velocity
!             UT     - Tangential velocity
!             DENSIT - Density
!             VISCOS - Viscosity
!             DIST   - Distance
!
!         Output:
!
!----------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION D_WALL_LAW( UFRIC,UT, DENSIT,  &
                    VISCOS,DIST,ROUGH )

      IMPLICIT NONE

      DOUBLE PRECISION DENSIT,VISCOS,DIST,ROUGH,UT,UFRIC,DFX,  &
      YPLUS

      DOUBLE PRECISION :: DKAPPA = 0.41D0, RAJA
      
      YPLUS=DENSIT*UFRIC*DIST/VISCOS

! Log-law:
!      RAJA=11.2658567D0 
!      IF(YPLUS.GE.RAJA) THEN
!         D_WALL_LAW=(1.0D0/DKAPPA)* &
!             ( DLOG(ROUGH*DENSIT*UFRIC*DIST/VISCOS) + 1.0D0 ) 
!      ELSE
!         D_WALL_LAW=DENSIT*DIST*2.0D0*UFRIC/VISCOS
!      ENDIF
 
! Reichardts law:
      D_WALL_LAW=DLOG(1.0D0 + 0.4D0*YPLUS)/DKAPPA  &
          + (0.4D0/DKAPPA)*YPLUS/(1.0D0 + 0.4D0*YPLUS) &
          + 7.8D0*( 1.0D0 - DEXP(-YPLUS/11.0D0) -   &
          (YPLUS/11.0D0)*DEXP(-0.33*YPLUS) )  &
          + 7.8D0*(YPLUS/11.0D0)  &
          *(  &
          DEXP(-YPLUS/11.0D0)-DEXP(-0.33*YPLUS)  &
          +0.33D0*YPLUS*DEXP(-0.33*YPLUS) &
          )

      RETURN
      END

!-----------------------------------------------------------------------------
!> Calculate the boundary values of turbulent kinetic energy
!> and its dissipation based on the wall law.
!
!         Input:
!             UT     - Tangential velocity of the previous iteration
!             DIST   - Distance from the wall
!             VISCOS - Viscosity
!             DENSIT - Density
!
!         Output:
!             TK     - Turbulent kinetic energy 
!             TEPS   - Turbulent kinetic energy dissipation
!
!----------------------------------------------------------------------------
      SUBROUTINE KEWALL (TK, TEPS, TOMG, UT, DIST, ROUGH, VISCOS, DENSIT )

      IMPLICIT NONE

      DOUBLE PRECISION TK, TEPS, TOMG, UT, DIST, VISCOS, DENSIT, ROUGH
      DOUBLE PRECISION UFRIC, DFX, UTLOCAL
      DOUBLE PRECISION :: CMYY   = 0.09D0
      DOUBLE PRECISION :: KARMAN = 0.41D0
      DOUBLE PRECISION :: SMALL  = 1.0D-10
      DOUBLE PRECISION :: OmegaWallPlus,Yplus, KsPlus, OmegaPlus, TomgL,TomgT,t,Alpha

      UTLOCAL = DMAX1( UT,SMALL )
      CALL SOLVE_UFRIC(DENSIT,VISCOS,DIST,ROUGH,UTLOCAL,UFRIC,DFX)

      yplus = densit*ufric*dist/viscos
      alpha = MIN(1.0d0, yplus/10)

      TK   = ( UFRIC**2 ) /  DSQRT( CMYY ) * alpha
      TEPS = ( UFRIC**3 ) / ( KARMAN * DIST ) * &
            MIN(1d0,alpha+0.2d0*Karman*(1-alpha**2)/SQRT(CMYY))

      OmegaPlus = 6/(0.072d0*yplus**2)
      TOMGL = Densit*Ufric**2*OmegaPlus/Viscos
      TOMGT = UFRIC / ( SQRT(CMYY) * KARMAN * DIST )

      IF (yplus<4 ) THEN
        TOMG = TOMGL
      ELSE IF ( yplus<32 ) THEN
        TOMG = SQRT(tomgL**2+tomgT**2)
      ELSE
        TOMG = TOMGT
      END IF

      RETURN
      END

!-----------------------------------------------------------------------------
