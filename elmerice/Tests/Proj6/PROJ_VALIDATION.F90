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
!/******************************************************************************
! *   Check the cordinate transformation implemented in the ProjUtils module
! *   Read an input file with x(m),y(m),longitude(degrees),latitude(degrees)
! *   Check the transformations (x,y) -> (Lon,Lat) -> (x,y)
! *   Result in Fatal errors if the differences are larger than a given treshold
! ******************************************************************************
! *
! *  Authors: Fabien Gillet-Chaulet
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *
! *  Original Date: 2/06/2026
! *
! *****************************************************************************/
    SUBROUTINE PROJ_VALIDATION_init( Model,Solver,dt,TransientSimulation )
      !------------------------------------------------------------------------------
      USE DefUtils

      IMPLICIT NONE
      !------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
      !------------------------------------------------------------------------------
      ! Local variables
      !------------------------------------------------------------------------------
      CHARACTER(LEN=MAX_NAME_LEN) :: Name
      LOGICAL :: Gotit

      Name = ListGetString( Solver % Values, 'Equation',GotIt)
      IF(.NOT. GotIt) Name = "Check"

      IF( .NOT. ListCheckPresent( Solver % Values,'Variable') ) THEN
        CALL ListAddString( Solver % Values,'Variable',&
            '-nooutput -global '//TRIM(Name)//'_var')
      END IF

      END SUBROUTINE


      SUBROUTINE PROJ_VALIDATION( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE ProjUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
      
      REAL(Kind=dp),DIMENSION(2) :: xy,xyr
      REAL(Kind=dp),DIMENSION(2) :: LonLat,LonLatr
      INTEGER,PARAMETER :: IO=12
      INTEGER :: nv
      INTEGER :: ok
      REAL(kind=dp) :: d
      REAL(kind=dp) :: errmax,errav
      REAL(kind=dp) :: dlmax,dlav
      REAL(Kind=dp),parameter :: Check=5.0e-4
      CHARACTER(LEN=MAX_NAME_LEN) :: InputFile
      LOGICAL :: Verbose,Found
      CHARACTER(LEN=*),PARAMETER :: Caller="PROJ_VALIDATION"

      InputFile = ListGetString(Solver % Values,'Input Filename',UnFoundFatal=.true. )

      open(IO,file=trim(InputFile),status = 'old',iostat = ok)

      CALL ProjINIT

      errmax=1.0d-12
      errav=0._dp
      dlav=0._dp
      dlmax=1.0d-12
      nv=0
      do while(ok == 0)
       !# reference x,y (m) and longitude,latitue (degrees)
       read(io,*,iostat = ok) xy(1:2),LonLat(1:2)
       if (ok == 0) THEN
          nv = nv + 1
          !xy -> lonlat
          CALL xy2LonLat(xy(1),xy(2),LonLatr(1),LonLatr(2))
          !lonlat -> xy
          CALL LonLat2xy(LonLatr(1),LonLatr(2),xyr(1),xyr(2))

          ! error xy -> LonLat -> xy
          d=sqrt((xy(1)-xyr(1))**2+(xy(2)-xyr(2))**2)
          errav=errav+d
          errmax=max(errmax,d)

          ! error xy -> LonLat   vs Ref. LonLat
          dlav=dlav+abs(LonLat(1)-LonLatr(1))+abs(LonLat(2)-LonLatr(2))
          dlmax=max(dlmax,abs(LonLat(1)-LonLatr(1)))
          dlmax=max(dlmax,abs(LonLat(2)-LonLatr(2)))

          !> Print errors
          CALL INFO(Caller,"#####################",Level=4)
          CALL INFO(Caller," - Point n "//I2S(nv),Level=4)
          WRITE(Message,*) "   xy in  : ",xy(1:2)
          CALL INFO(Caller,Message,Level=4)
          WRITE(Message,*) "   xy out : ",xyr(1:2)
          CALL INFO(Caller,Message,Level=4)
          WRITE(Message,*) "   LonLat in : ",LonLat(1:2)
          CALL INFO(Caller,Message,Level=4)
          WRITE(Message,*) "   LonLat out : ",LonLatr(1:2)
          CALL INFO(Caller,Message,Level=4)

       endif
      end do
     close(IO)

      !> Get Avaraged results
      errav=errav/nv
      dlav=dlav/(2*nv)

      CALL INFO(Caller,"#####################",Level=3)
      CALL INFO(Caller,"Averaged errors:",Level=3)
      WRITE(Message,*) "errav=",errav
      CALL INFO(Caller,Message,Level=3)
      WRITE(Message,*) "errmax=",errmax
      CALL INFO(Caller,Message,Level=3)
      WRITE(Message,*) "dlav=",dlav
      CALL INFO(Caller,Message,Level=3)
      WRITE(Message,*) "dlmax=",dlmax
      CALL INFO(Caller,Message,Level=3)
      CALL INFO(Caller,"#####################",Level=3)

      IF ((errav > check).OR.(errmax > check).OR.&
         (dlav > check).OR.(dlmax > check)) &
         CALL FATAL(Caller,"Test did not pass")

      ! Test passed if we are here; otherwise would have stop with a fatal
      Solver % Variable % Norm = 1._dp
      Solver % Variable % Values = 1._dp

      END SUBROUTINE PROJ_VALIDATION
