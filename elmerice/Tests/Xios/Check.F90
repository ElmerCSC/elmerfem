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
! *  here we check the connectivity.
! *   Netcdf file name and variable names are hard-coded and thus might
! *   be specific to this test
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

      INTEGER :: nface,nvertex
      INTEGER :: t
      INTEGER :: n
      INTEGER :: EIndex
      INTEGER, Allocatable :: connectivity(:,:),Vertices(:)
      LOGICAL :: Parallel

      REAL(KIND=dp) :: Norm=1._dp,Change=0._dp


      Parallel=(ParEnv%PEs>1)

      ! open output_ugrid.nc
      NetCDFstatus = NF90_OPEN("output_ugrid.nc",NF90_NOWRITE,ncid)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
         CALL Fatal("Check",'Unable to open NETCDF File')
      END IF

      NetCDFstatus = nf90_inq_dimid(ncid, "nmesh2D_face" , varid)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
         CALL Fatal("Check",'Unable to get <nmesh2D_face>')
      END IF
      NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=nface)

      NetCDFstatus = nf90_inq_dimid(ncid, "nmesh2D_vertex" , varid)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
         CALL Fatal("Check",'Unable to get <nmesh2D_vertex>')
      END IF
      NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=nvertex)

      ALLOCATE(connectivity(nvertex,nface),Vertices(Model%MaxElementNodes))

      NetCDFstatus = nf90_inq_varid(ncid,"mesh2D_face_nodes",varid)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal("Check",'Unable to get <mesh2D_face_nodes> id')
      END IF

      NetCDFstatus = nf90_get_var(ncid, varid,connectivity)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal("Check",'Unable to get <mesh2D_face_nodes>')
      END IF

      NetCDFstatus = NF90_Close(ncid)

      ! connectivity is zero-based
      connectivity=connectivity+1

      DO t=1,Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes(Element)

        IF (Parallel) THEN
          EIndex=Element%GElementIndex
          Vertices(1:n)=Model%Mesh % ParallelInfo % GlobalDOFs(Element%NodeIndexes(1:n))
        ELSE
          EIndex=Element%ElementIndex
          Vertices(1:n)=Element%NodeIndexes(1:n)
        ENDIF
        IF (.NOT.(ALL(Vertices(1:n).EQ.connectivity(:,EIndex)))) THEN
           PRINT *,"ElementIndex",Element%ElementIndex
           PRINT *,"NodeIndexes",Element%NodeIndexes(:)
           PRINT *,"connectivity",connectivity(:,Element%ElementIndex)
           CALL Fatal("Check",'Error in connectivity array')
        END IF
      END DO

      ! Test passed if we are here; otherwise would have stop with a fatal
      Solver % Variable % Norm = 1._dp
      Solver % Variable % Values = 1._dp

      !! End
      DEALLOCATE(connectivity,Vertices)

      END SUBROUTINE Check
