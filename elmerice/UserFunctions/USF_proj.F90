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
! *  Authors: fabien Gillet-Chaulet
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 3 May 2022
! *
! ******************************************************************************/
!--------------------------------------------------------------------------------
!>  Generic ElmerIce user function for geographic projections
!>  Use ProjUtils in under elmerice/Utils
!--------------------------------------------------------------------------------
      FUNCTION xy2Lon(Model,nodenumber,VarIn) RESULT(VarOut)
      USE ProjUtils
      implicit none
      !-----------------
      TYPE(Model_t) :: Model
      INTEGER :: nodenumber
      REAL(kind=dp) :: VarIn(2) !x,y
      REAL(kind=dp) :: VarOut   !Longitude
      REAL(kind=dp) :: x,y,Lon,Lat
        
        x=VarIn(1)
        y=VarIn(2)
        CALL xy2LonLat(x,y,Lon,Lat) 
        VarOut=Lon

      End FUNCTION xy2Lon

      FUNCTION xy2Lat(Model,nodenumber,VarIn) RESULT(VarOut)
      USE ProjUtils
      implicit none
      !-----------------
      TYPE(Model_t) :: Model
      INTEGER :: nodenumber
      REAL(kind=dp) :: VarIn(2) !x,y
      REAL(kind=dp) :: VarOut   !Latitude
      REAL(kind=dp) :: x,y,Lon,Lat

        x=VarIn(1)
        y=VarIn(2)
        CALL xy2LonLat(x,y,Lon,Lat) 
        VarOut=Lat

      End FUNCTION xy2Lat

      FUNCTION LonLat2x(Model,nodenumber,VarIn) RESULT(VarOut)
      USE ProjUtils
      implicit none
      !-----------------
      TYPE(Model_t) :: Model
      INTEGER :: nodenumber
      REAL(kind=dp) :: VarIn(2) !Lon,Lat
      REAL(kind=dp) :: VarOut   ! x
      REAL(kind=dp) :: x,y,Lon,Lat

        Lon=VarIn(1)
        Lat=VarIn(2)
        CALL LonLat2xy(Lon,Lat,x,y) 
        VarOut=x

      End FUNCTION LonLat2x

      FUNCTION LonLat2y(Model,nodenumber,VarIn) RESULT(VarOut)
      USE ProjUtils
      implicit none
      !-----------------
      TYPE(Model_t) :: Model
      INTEGER :: nodenumber
      REAL(kind=dp) :: VarIn(2) !Lon,Lat
      REAL(kind=dp) :: VarOut   !y
      REAL(kind=dp) :: x,y,Lon,Lat

        Lon=VarIn(1)
        Lat=VarIn(2)
        CALL LonLat2xy(Lon,Lat,x,y) 
        VarOut=y

      End FUNCTION LonLat2y
