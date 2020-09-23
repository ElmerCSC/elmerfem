!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Samuel Cook 
! *  Email:   sc690@cam.ac.uk
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - Scientific Computing Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 30 July 2018
! *
! *****************************************************************************/
!> A minor alteration to the Restart() subroutine in ElmerSolver.F90 to allow a
!> secondary hydrology mesh to be restarted without getting caught up in all the
!> model setup routines
   SUBROUTINE HydroRestart( Model,Solver,Timestep,TransientSimulation )
!******************************************************************************
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: Timestep
!     INPUT: Timestep size for time dependent simulations
!
!******************************************************************************
!------------------------------------------------------------------------------
     USE DefUtils

     IMPLICIT NONE
 
     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp) :: Timestep
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Solver_t), POINTER :: HydroSolver, ChannelSolver, ThicknessSolver
     LOGICAL :: Gotit
     INTEGER :: k, i, HPSolver, ChanSolver, ThickSolver, SolverID
     CHARACTER(LEN=MAX_STRING_LEN) :: OutputName, RestartFile
!------------------------------------------------------------------------------

     RestartFile = ListGetString( CurrentModel % Simulation, &
         'Restart File', GotIt )

     DO i=1,Model % NumberOfSolvers
       IF(Model % Solvers(i) % Variable % Name == 'hydraulic potential') THEN
         HPSolver = i
         HydroSolver => Model % Solvers(HPSolver)
         EXIT
       END IF
     END DO
     DO i=1,Model % NumberOfSolvers
       IF(Model % Solvers(i) % Variable % Name == 'channel area') THEN
         ChanSolver = i
         ChannelSolver => Model % Solvers(ChanSolver)
         EXIT
       END IF
     END DO
     DO i=1,Model % NumberOfSolvers
       IF(Model % Solvers(i) % Variable % Name == 'sheet thickness') THEN
         ThickSolver = i
         ThicknessSolver => Model % Solvers(ThickSolver)
         EXIT
       END IF
     END DO

     SolverID = Solver % SolverId

     IF ( GotIt ) THEN
       !Because all the actual variable output will be on this solver mesh, not
       !on the identical NO meshes defined in the other solvers for
       !initialisation and result output purposes
       OutputName = TRIM(HydroSolver % Mesh % Name) // '/' // TRIM(RestartFile)
       k = ListGetInteger( CurrentModel % Simulation,'Restart Position',GotIt, &
                         minv=0 )

       IF ( ParEnv % PEs > 1 ) &
         OutputName = TRIM(OutputName) // '.' // TRIM(i2s(ParEnv % MyPe))

       CALL ListPushNameSpace('hp:')
       Mesh => HydroSolver % Mesh
       CALL LoadRestartFile( OutputName,k,Mesh, SolverId = SolverID )
       CALL ListPopNameSpace()

       CALL ListPushNameSpace('channel:')
       Mesh => ChannelSolver % Mesh
       CALL LoadRestartFile( OutputName,k,Mesh, SolverId = SolverID )
       CALL ListPopNameSpace()

       CALL ListPushNameSpace('sheet:')
       Mesh => ThicknessSolver % Mesh
       CALL LoadRestartFile( OutputName,k,Mesh, SolverId = SolverID )
       CALL ListPopNameSpace()
  
     END IF
     NULLIFY(HydroSolver, ChannelSolver, ThicknessSolver, Mesh)
!------------------------------------------------------------------------------
   END SUBROUTINE HydroRestart
!------------------------------------------------------------------------------
