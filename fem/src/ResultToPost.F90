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

!> \ingroup Programs 
!> \{

!> \defgroup ResultToPost Program ResultToPost
!> \{

!------------------------------------------------------------------------------
!>  Stand-alone program for Elmer results file to Elmer post processing 
!> file conversion.
!------------------------------------------------------------------------------
   PROGRAM ResultToPost
!------------------------------------------------------------------------------

     USE MainUtils

!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     INTEGER :: i,j,k,n,l,t,k1,k2,iter,Ndeg,Time,istat,nproc

     REAL(KIND=dp) :: s,dt
     REAL(KIND=dp), TARGET :: SimulationTime(1)

     TYPE(Element_t),POINTER :: CurrentElement
     REAL(KIND=dp), POINTER ::  TimeVariable(:)

     LOGICAL :: GotIt,TransientSimulation,LastSaved

     INTEGER :: TimeIntervals,interval,timestep,SavedSteps,TimeStepMaxIter

     REAL(KIND=dp), POINTER :: TimestepSizes(:,:)
     INTEGER, POINTER :: Timesteps(:),OutputIntervals(:)

     LOGICAL :: KEModelSolved = .FALSE., SteadyStateReached = .FALSE.

     TYPE(ElementType_t),POINTER :: elmt

     TYPE(ParEnv_t), POINTER :: ParallelEnv

     CHARACTER(:), ALLOCATABLE :: eq,OutputFile,PostFile,OutputName,PostName

     CHARACTER(LEN=MAX_NAME_LEN) :: ModelName

     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Variable_t), POINTER :: Var
     TYPE(Solver_t), POINTER :: Solver

!------------------------------------------------------------------------------
!    Read input file name and whether parallel execution is requested
!------------------------------------------------------------------------------

     OPEN( 1,file='ELMERSOLVER_STARTINFO')
       READ(1,'(a)') ModelName
     CLOSE(1)

!------------------------------------------------------------------------------
!    If parallel execution requested, initialize parallel environment
!------------------------------------------------------------------------------
     ParEnv % PEs  = 1 
     ParEnv % MyPE = 0

!------------------------------------------------------------------------------
!    Read element definition file, and initialize element types
!------------------------------------------------------------------------------
     CALL InitializeElementDescriptions
!------------------------------------------------------------------------------
!    Read Model and mesh from Elmer mesh data base
!------------------------------------------------------------------------------
     IF ( ParEnv % MyPE == 0 ) THEN
       PRINT*, ' '
       PRINT*, ' '
       PRINT*, '-----------------------'
       PRINT*,'Reading Model ...       '
     END IF

     CurrentModel => LoadModel( ModelName,.FALSE.,ParEnv % PEs,ParEnv % MyPE )

     IF ( ParEnv % MyPE == 0 ) THEN
       PRINT*,'... Done               '
       PRINT*, '-----------------------'
     END IF

!------------------------------------------------------------------------------
!    1 point integration rule is not enough for complete integration of
!    time derivative terms (or any other zeroth order space derivative terms for
!    that  matter, .i.e. the KE model) for linear elements. It is enough
!    (for the time derivative term), if mass matrix  lumping is used though.
!    Anyway, to be on the safe side, if the simulation  is time dependent,
!    change the number of integration points  here.
!
!    NOTE: THIS DOESN'T FIX THE PROBLEM FOR THE KE Model
!------------------------------------------------------------------------------
     eq = ListGetString( CurrentModel % Simulation, 'Simulation Type' )
     TransientSimulation = .FALSE.

     IF ( eq == 'transient' ) THEN
       TransientSimulation= .TRUE.

       elmt => GetElementType( 303 )
       IF ( elmt % GaussPoints == 1 ) elmt % GaussPoints = 3

       elmt => GetElementType( 504 )
       IF ( elmt % GaussPoints == 1 ) elmt % GaussPoints = 4
     END IF
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!    Figure out requested coordinate system
!------------------------------------------------------------------------------
     eq = ListGetString( CurrentModel % Simulation, 'Coordinate System', GotIt )

     Coordinates = Cartesian
     IF ( eq == 'cartesian 2d' ) THEN

       CurrentModel % DIMENSION = 2
       Coordinates = Cartesian

     ELSE IF ( eq == 'cartesian 3d' ) THEN

       CurrentModel % DIMENSION = 3
       Coordinates = Cartesian

     ELSE IF ( eq == 'axi symmetric' ) THEN

       CurrentModel % DIMENSION = 2
       Coordinates = AxisSymmetric

     ELSE IF( eq == 'cylindric symmetric' ) THEN

       CurrentModel % DIMENSION = 2
       Coordinates = CylindricSymmetric

     ELSE IF( eq == 'cylindrical' ) THEN

       CurrentModel % DIMENSION = 3
       Coordinates = Cylindric

     ELSE IF( eq == 'polar 2d' ) THEN

       CurrentModel % DIMENSION = 2
       Coordinates = Polar

     ELSE IF( eq == 'polar 3d' ) THEN

       CurrentModel % DIMENSION = 3
       Coordinates = Polar

     ELSE

       PRINT*,'Solver: ERROR: Unknown global coordinate system: ',TRIM(eq),' Aborting'
       STOP

     END IF

!------------------------------------------------------------------------------
!   Figure out what (flow,heat,stress,...) should be computed, and get
!   memory for the dofs
!------------------------------------------------------------------------------
     DO i=1,CurrentModel % NumberOfSolvers
       eq = ListGetString( CurrentModel % Solvers(i) % Values,'Equation' )

       Solver => CurrentModel % Solvers(i)
       CALL AddEquationBasics( Solver, eq, TransientSimulation )
       CALL AddEquationSolution( Solver, TransientSimulation )
     END DO

!------------------------------------------------------------------------------
!    Add coordinates to list of variables so that coordinate dependent
!    parameter computing routines can ask for them...
!------------------------------------------------------------------------------
     TimeVariable => SimulationTime
     Mesh => CurrentModel % Meshes 
     DO WHILE( ASSOCIATED( Mesh ) )
       CALL VariableAdd(Mesh % Variables,Mesh,NULL(),'Coordinate 1',1,Mesh % Nodes % x )
       CALL VariableAdd(Mesh % Variables,Mesh,NULL(),'Coordinate 2',1,Mesh % Nodes % y )
       CALL VariableAdd(Mesh % Variables,Mesh,NULL(),'Coordinate 3',1,Mesh % Nodes % z )
       CALL VariableAdd( Mesh % Variables,Mesh, NULL(),'Time',1,TimeVariable )
       Mesh => Mesh % Next
     END DO


!------------------------------------------------------------------------------
!    Convert results file to post processing file, if requested
!------------------------------------------------------------------------------
     PostFile = ListGetString( CurrentModel % Simulation,'Post File',GotIt )
     OutputFile = ListGetString(CurrentModel % Simulation,'Output File',GotIt)
     IF ( GotIt ) THEN
       IF ( ParEnv % PEs > 1 ) THEN
         PostFile = PostFile // '.' // i2s(ParEnv % MyPE)
       END IF
       Mesh => CurrentModel % Meshes
       DO WHILE( ASSOCIATED( Mesh ) )
         IF ( LEN_TRIM( Mesh % Name ) > 0 )  THEN
           OutputName = TRIM(Mesh % Name) // '/' // TRIM(OutputFile)
           PostName   = TRIM(Mesh % Name) // '/' // TRIM(PostFile)
         ELSE
           PostName   = PostFile
           OutputName = OutputFile
         END IF
         CALL SetCurrentMesh( CurrentModel, Mesh )
         CALL WritePostFile( PostName,OutputName,CurrentModel,10000 )
         Mesh => Mesh % Next
       END DO
     END IF
!------------------------------------------------------------------------------
!  THIS IS THE END (...,at last, the end, my friend,...)
!------------------------------------------------------------------------------
     PRINT*,'*** Result To Post: ALL DONE ***'

!------------------------------------------------------------------------------
  END PROGRAM ResultToPost
!------------------------------------------------------------------------------

!> \}
!> \}
