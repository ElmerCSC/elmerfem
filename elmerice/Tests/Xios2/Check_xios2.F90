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
! *  Subroutine to check that output from XIOS was successful
! *  
! ******************************************************************************
! *
! *  Authors: Fabien Gillet-Chaulet
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *
! *  Original Date: 4/05/2022
! *
! *****************************************************************************/
      SUBROUTINE Check_init( Model,Solver,dt,TransientSimulation )
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
      
      SUBROUTINE Check( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE Netcdf

      IMPLICIT NONE
      !------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      TYPE(Model_t) :: Model

      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
      !------------------------------------------------------------------------------
      ! Local variables
      !------------------------------------------------------------------------------
      TYPE(Element_t),POINTER :: Element
      INTEGER :: NetcdfStatus,varid,ncid

      INTEGER :: ntime
      INTEGER :: t
      REAL(KIND=dp), Allocatable :: Volume(:)
      REAL(KIND=dp) :: area
      REAL(KIND=dp) :: tv,av,vv
      LOGICAL :: Parallel

      REAL(KIND=dp),PARAMETER :: Tol=1.0d-6
      LOGICAL :: success


      Parallel=(ParEnv%PEs>1)

      ! open output_ugrid.nc
      NetCDFstatus = NF90_OPEN("output_Scalar.nc",NF90_NOWRITE,ncid)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
         CALL Fatal("Check",'Unable to open NETCDF File')
      END IF

      NetCDFstatus = nf90_inq_dimid(ncid, "time" , varid)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
         CALL Fatal("Check",'Unable to get time')
      END IF
      NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=ntime)

      ALLOCATE(Volume(ntime))

      NetCDFstatus = nf90_inq_varid(ncid,"volume",varid)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal("Check",'Unable to get volume id')
      END IF

      NetCDFstatus = nf90_get_var(ncid, varid,Volume)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal("Check",'Unable to get volume')
      END IF

      NetCDFstatus = nf90_inq_varid(ncid,"area",varid)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal("Check",'Unable to get area id')
      END IF

      NetCDFstatus = nf90_get_var(ncid, varid,area)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal("Check",'Unable to get area')
      END IF

      NetCDFstatus = NF90_Close(ncid)


      OPEN(12,File="f.dat")
      DO t=1,ntime
         read(12,*) tv,av,vv
         success=( ((abs(av-area)/area).LT.Tol).OR.&
                 ((abs(vv-volume(t))/volume(t)).LT.Tol) )
         IF (.NOT.success) &
           CALL FATAL("Check","Pb with area")
      END DO
      CLOSE(12)
      ! Test passed if we are here; otherwise would have stop with a fatal
      Solver % Variable % Norm = 1._dp
      Solver % Variable % Values = 1._dp

      !! End
      DEALLOCATE(Volume)

      END SUBROUTINE Check
