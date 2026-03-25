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
! *  Author: F. Gillet-Chaulet (IGE)
! *  Email:  fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:    http://elmerice.elmerfem.org
! *
! *  Original Date: 13-07-2017, 
! *****************************************************************************
!!!  
!!!  Mesh adaptation solver based on the MMG2D library
!!!
!!! Use the master branch: (tested with update on 13/07/2017) 
!!! https://github.com/MmgTools/mmg
!!!  USE CURRENT MESH AND
!!!  IMPLEMENT ISO (Scalar variable=mesh file)
!!!  AND ANISO Remeshing (variable DIM=3; Metric matrix (1,1),(2,2),(1,2))
!!!  SAVE Output mesh in elmer format
!!!   -  SERIAL ONLY
!!!   -  2D ONLY
!!!  
!!!   - replace current mesh by NewMesh to have remeshing during the simulation:
!!!     
!!!    TODO: Variable RELEASE_MESH (default false) allow to release prev
!mesh if true (seems that not everything is deallocated); 
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MMG2DSolver_init( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      IMPLICIT NONE
      !------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
      !--------------------------------------------------------------------------
      CHARACTER(LEN=MAX_NAME_LEN) :: Name
      TYPE(ValueList_t), POINTER :: SolverParams
      LOGICAL :: GotIt
  
      SolverParams => Solver % Values 

      Name = ListGetString( SolverParams, 'Equation',GotIt)
      IF( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
        CALL ListAddString( SolverParams,'Variable',&
           '-nooutput '//TRIM(Name)//'_var')
      ENDIF

      CALL ListAddNewLogical(SolverParams,'Optimize Bandwidth',.FALSE.)

      END SUBROUTINE MMG2DSolver_init
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE MMG2DSolver( Model,Solver,dt,TransientSimulation)
!------------------------------------------------------------------------------
        USE Types
        USE Lists
        USE MeshUtils
        USE MeshRemeshing
!------------------------------------------------------------------------------
      IMPLICIT NONE

#ifdef HAVE_MMG
#include "mmg/common/libmmgtypesf.h"
#endif
      
      TYPE(Model_t) :: Model
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Solver_t), POINTER ::PSolver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

#ifdef HAVE_MMG
      TYPE(Mesh_t),POINTER :: Mesh
      TYPE(Mesh_t),POINTER :: NewMesh,PrevMesh
      TYPE(Variable_t),POINTER :: MeshSize
      TYPE(ValueList_t), POINTER :: SolverParams
      TYPE(Element_t),POINTER :: Element

      REAL(KIND=dp) :: Hmin, Hmax,hsiz
      REAL(KIND=dp) :: Change,Tolerance
 
      INTEGER, SAVE :: MeshNumber=0

      INTEGER :: ier
      INTEGER :: ii,kk
      INTEGER :: MinIter
      INTEGER, SAVE :: nVisit=0

      CHARACTER(len=300) :: filename
      CHARACTER(LEN=1024) :: Path
      CHARACTER(LEN=MAX_NAME_LEN) :: OutPutFileName
      CHARACTER(LEN=MAX_NAME_LEN) :: MeshSizeName
      CHARACTER(LEN=MAX_NAME_LEN) :: DefaultFileName='MMGMesh' 

      LOGICAL :: DEBUG=.FALSE.,Scalar,Found
      LOGICAL :: SAVE_MMG_INI=.FALSE.,SAVE_MMG_FINAL=.FALSE.,SAVE_ELMER_MESH=.TRUE.
      LOGICAL :: RELEASE_MESH=.FALSE.
      LOGICAL :: UnFoundFatal=.TRUE.
      LOGICAL :: IncrementMeshNumber
      LOGICAL :: UniformSize
      LOGICAL :: TestConvergence
      LOGICAL :: Parallel

      CALL Info('MMGSolver','Starting 2D mesh refinement',Level=5)
              Mesh => Solver % Mesh
      
      Parallel = (ParEnv % PEs > 1)

      nVisit = nVisit+1

      IF ((Parallel).AND.(nVisit>1)) THEN
        CALL Fatal('MMGSolver',&
            'As interpolation is not ready in Parallel; no point &
            &to visit this solver more than once?')
      END IF
        
      SolverParams => Solver % Values 
      IF(.NOT. ASSOCIATED(SolverParams)) THEN
        CALL Fatal('MMGSolver','Solver paramaters not associated!')
      END IF
      
      hsiz = ListGetConstReal( SolverParams, 'hsiz', UniformSize)
      
      ! NewMesh Name
      OutPutFileName = ListGetString(SolverParams,'Output file name',Found)
      IF (.NOT.Found) OutPutFileName = DefaultFileName

      IncrementMeshNumber = ListGetLogical(SolverParams,'Increment Mesh Number',Found)
      IF (.NOT.Found) IncrementMeshNumber = .TRUE.
      IF( IncrementMeshNumber ) MeshNumber = MeshNumber + 1
      
      ! Get serial mesh
      IF (Parallel) THEN
        Mesh => GATHER_PMesh(Solver % Mesh)
      ELSE
        Mesh => Solver % Mesh
      END IF
      
      ! GET REQUIRED VARIABLE
      IF (.NOT.UniformSize) THEN
        MeshSizeName = ListGetString( SolverParams, &
                        'Metric Variable Name',UnFoundFatal=UnFoundFatal )
        MeshSize => VariableGet(Solver%Mesh%Variables,TRIM(MeshSizeName),UnFoundFatal=UnFoundFatal)
        IF (MeshSize%DOFS /= 1 .AND. MeshSize%DOFS /= 3) THEN
          CALL Fatal('MMGSolver','Variable <ElementSize> should have 1 or 3 DOFs')
        END IF
        Scalar = (MeshSize % DOFS == 1)
        IF (Parallel) MeshSize => GATHER_PVar(Solver % Mesh,Mesh,MeshSize,MeshSizeName)
      END IF

      IF (Parallel .AND. ParEnv % MyPE /= 0) RETURN
             
      !INITIALISATION OF MMG MESH AND SOL STRUCTRURES
      CALL Info('MMGSolver','Initialization of MMG',Level=20)
      mmgMesh = 0
      mmgSol  = 0
      CALL MMG2D_Init_mesh(MMG5_ARG_start, &
             MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
             MMG5_ARG_end)
      
      ! SET PARAMETERS 
      CALL Info('MMGSolver','Set MMG2D Parameters',Level=20)
      CALL SET_MMG2D_PARAMETERS(SolverParams)

      !> COPY MESH TO MMG FORMAT      
      CALL Info('MMGSolver','Copy Mesh to MMG format',Level=20)
      CALL SET_MMG2D_MESH(Mesh)
      
      IF (.NOT.UniformSize) THEN
        CALL Info('MMGSolver','Set the SOL',Level=20)      
        CALL SET_MMG2D_SOL(Mesh,MeshSize,Scalar)

        ! (not mandatory): check if the number of given entities match with mesh size
        CALL Info('MMGSolver','Check the mesh data',Level=20)              
        CALL MMG2D_Chk_meshData(mmgMesh,mmgSol,ier)
        IF ( ier == 0 ) CALL Fatal('MMGSolver',&
                   'CALL TO MMG2D_Chk_meshData FAILED')
      END IF

      ! (not mandatory): save the mesh
      IF (SAVE_MMG_INI) THEN
        filename="MMGini.mesh"
        CALL MMG2D_SaveMesh(mmgMesh,TRIM(ADJUSTL(filename)), &
            LEN(TRIM(ADJUSTL(filename))),ier)
        filename="MMGini.sol"
        CALL MMG2D_SaveSol(mmgMesh,mmgSol,TRIM(ADJUSTL(filename)), &
            LEN(TRIM(ADJUSTL(filename))),ier)
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MESH OPTIMIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL MMG2D_mmg2dlib(mmgMesh,mmgSol,ier)
      IF ( ier == MMG5_STRONGFAILURE ) THEN
        PRINT*,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH"
        STOP MMG5_STRONGFAILURE
      ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
        PRINT*,"BAD ENDING OF MMG3DLIB"
      ENDIF
      IF (DEBUG) PRINT *,'--**-- MMG2D_mmg2dlib DONE'
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! (not mandatory): save the new mesh
      IF (SAVE_MMG_FINAL) THEN
        filename="MMGsortie.mesh"
        CALL MMG2D_SaveMesh(mmgMesh,TRIM(ADJUSTL(filename)), &
            LEN(TRIM(ADJUSTL(filename))),ier)
      END IF
      
!! GET THE NEW MESH
      NewMesh => GET_MMG2D_MESH(MeshNumber,OutputFilename)
      
!! MIMIC COMPUTE CHANGE STYLE
      Change=2.*(NewMesh%NumberOfNodes-Mesh%NumberOfNodes)/float(NewMesh%NumberOfNodes+Mesh%NumberOfNodes)
      Change=ABS(Change)
      WRITE( Message, '(a,i0,g15.8,a)') &
          'SS (ITER=1) (NRM,RELC): (',NewMesh%NumberOfNodes, Change,&
          ' ) :: MMG2DSolver'
      CALL Info( 'ComputeChange', Message, Level=3 )
!! IF CONVERGENCE CRITERIA IS GIVEN AND REACHED, SET EXIT CONDITION TO TRUE
!! (can not used the internal criteria as solver is executed After timestep)
      Tolerance = &
          ListGetCReal( SolverParams,'Steady State Convergence Tolerance',Found)
      IF (Found) THEN
        MinIter=&
            ListGetInteger(SolverParams,'Steady State Min Iterations',Found)
        IF (Found) TestConvergence = ( nVisit >= MinIter )
        IF ((TestConvergence).AND.(Change < Tolerance)) THEN
          CALL ListAddConstReal(Model % Simulation,'Exit Condition',1.0_dp)
        END IF
      END IF

!  SAVE MESH TO DISK
      IF (SAVE_ELMER_MESH) THEN
        Path = TRIM(NewMesh % Name)
        CALL MakeDirectory( TRIM(path) // CHAR(0) )
        CALL WriteMeshToDisk( NewMesh, Path )
      END IF

      ! Interpolate only in serial
      IF (.NOT.Parallel) THEN
        
! Get previous Mesh to do interpolation
!---------------------------------------
        PrevMesh => Model % Meshes
      
! Add Variable to NewMesh Structure + Model%mesh => NewMesh
!----------------------------------------------------------
        CALL AddValueToMesh(NewMesh,PrevMesh) 

! Add the new mesh to the global list of meshes
!----------------------------------------------
        Model % Meshes  => NewMesh
      
! Update Solver Meshes
!---------------------
        DO ii=1,Model % NumberOfSolvers
          PSolver => Model%Solvers(ii)
          IF (.NOT.ASSOCIATED(PSolver)) CYCLE
          IF (.NOT.ASSOCIATED(PSolver%Matrix)) THEN
            PSolver % Mesh => NewMesh
            IF (ASSOCIATED(PSolver%Variable)) THEN
              PSolver % Variable => VariableGet( NewMesh % Variables,&
                  Solver % Variable % Name, ThisOnly = .TRUE. )
            ENDIF
          ELSE
            CALL UpdateSolverMesh(PSolver,NewMesh) ! call distrib version in MeshUtils.F90
            WRITE(Message,'(A,A,A)') 'Updated for variable: ',TRIM(PSolver%Variable%Name)
            CALL INFO('MMGSolver',TRIM(Message),level=5)
          END IF
        END DO

! Release previous mesh (should be optional)
!-------------------
        RELEASE_MESH=ListGetLogical(SolverParams,'Release previous mesh',Found)
        IF (RELEASE_MESH) THEN
          CALL ReleaseMesh(PrevMesh)
          Deallocate(PrevMesh)
          Model % Meshes % Next => NULL()
        END IF

!! End Serial interpolation
!-----------------
      END IF

      NewMesh % Changed = .TRUE.
! Print info
!-----------
      write(Message,'(A)') '--**--New mesh ready'
      CALL INFO('MMGSolver',trim(Message),level=5)
      write(Message,'(A,I0)') '--**--    NumberOfNodes:',NewMesh%NumberOfNodes
      CALL INFO('MMGSolver',trim(Message),level=5)
      write(Message,'(A,I0)') '--**--    NumberOfBulkElements: ',NewMesh%NumberOfBulkElements
      CALL INFO('MMGSolver',trim(Message),level=5)
      write(Message,'(A,I0)') '--**--    NumberOfEdges: ',NewMesh%NumberOfEdges
      CALL INFO('MMGSolver',trim(Message),level=5)
      write(Message,'(A,I0)') '--**--    NumberOfBoundaryElements: ',NewMesh%NumberOfBoundaryElements
      CALL INFO('MMGSolver',trim(Message),level=5)
      
!!!! Free the MMG3D5 structures
      CALL MMG2D_Free_all(MMG5_ARG_start, &
          MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
          MMG5_ARG_end)      
     

!-------------------------------------
!! Subroutine
!-------------------------------------
    CONTAINS

      !Gather Parallel mesh to Part 0
      FUNCTION GATHER_PMesh(Mesh) RESULT(GatheredMesh)
        USE MeshPartition
        IMPLICIT NONE
        TYPE(Mesh_t), POINTER :: Mesh,GatheredMesh
        INTEGER :: n,allocstat

        n = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

        IF( ASSOCIATED( Mesh % RePartition ) ) THEN
          IF( SIZE( Mesh % RePartition ) < n ) THEN
            DEALLOCATE(Mesh % RePartition)
            Mesh % RePartition => Null()
          END IF
        END IF

        IF(.NOT. ASSOCIATED( Mesh % RePartition ) ) THEN
          ALLOCATE( Mesh % RePartition( n ), STAT = allocstat)
          IF( allocstat /= 0 ) THEN
            CALL Fatal('GATHER_PMesh','Allocation error')
          END IF
        END IF

        Mesh % RePartition(:) = 1
        GatheredMesh => RedistributeMesh(Model, Mesh, .TRUE., .FALSE.)

        IF (DEBUG) THEN
          IF (ParEnv%MyPE==0) THEN 
            PRINT *,"-- GatheredMesh --"
            PRINT *,ParEnv%MyPE,GatheredMesh%NumberOfNodes,&
                GatheredMesh%NumberOfBulkElements,&
                GatheredMesh%NumberOfBoundaryElements
            GatheredMesh % Name = "Gathered"
            Path = TRIM(GatheredMesh % Name)
            CALL MakeDirectory( TRIM(path) // CHAR(0) )
            CALL WriteMeshToDisk( GatheredMesh, Path )
          END IF
        ENDIF
      END FUNCTION GATHER_PMesh

      !Gather Parellel variable for the mesh size
      FUNCTION GATHER_PVar(Mesh,GatheredMesh,Var,MeshSizeName) RESULT(GVar)
        IMPLICIT NONE
        TYPE(Variable_t),POINTER :: Var,GVar
        TYPE(Mesh_t),POINTER :: Mesh,GatheredMesh
        CHARACTER(LEN=MAX_NAME_LEN) :: MeshSizeName
        REAL(KIND=dp),POINTER :: Values(:)
        INTEGER,POINTER :: Perm(:)
        INTEGER,ALLOCATABLE :: node_count(:),disps(:),PGDOFs_send(:)
        REAL(KIND=dp),ALLOCATABLE :: PVar_send(:),local_var(:)

        INTEGER :: DOFs
        INTEGER :: ntot
        INTEGER :: ierr
        INTEGER :: i,k

        DOFs = Var % DOFs

        ALLOCATE(node_count(ParEnv%PEs),disps(ParEnv%PEs),local_var(Mesh%NumberOfNodes))
        node_count = 0

        ! Total number of nodes (including shared nodes)
        CALL MPI_Gather(Mesh%NumberOfNodes, 1, MPI_INTEGER, node_count, 1, &
            MPI_INTEGER, 0, ELMER_COMM_WORLD,ierr)
        ntot=SUM(node_count)
        
        ! Create disps table for GatherV
        ! and do some allocation
        IF(ParEnv % MyPE == 0) THEN

          disps(1) = 0
          DO i=2,ParEnv%PEs
            disps(i) = disps(i-1) + node_count(i-1)
          END DO

          ALLOCATE(PGDOFs_send(ntot),PVar_send(ntot))
          ALLOCATE(Values(DOFs*GatheredMesh % NumberOfNodes),Perm(GatheredMesh % NumberOfNodes))
        ELSE
          ALLOCATE(PGDOFs_send(1),PVar_send(1))
        END IF

        ! Gather Global DOFs
        CALL MPI_GatherV(Mesh % ParallelInfo % GlobalDOFs, &
            Mesh%NumberOfNodes, MPI_INTEGER, PGDOFs_send, node_count, &
            disps, MPI_INTEGER, 0, ELMER_COMM_WORLD ,ierr)

        ! Gather Variable values
        DO k=1,DOFs
          DO i=1,Mesh%NumberOfNodes
            local_var(i)=Var % Values(DOFs*(Var%Perm(i)-1)+k)
          END DO
          CALL MPI_GatherV(local_var, &
              Mesh%NumberOfNodes, MPI_DOUBLE_PRECISION, PVar_send, node_count, &
              disps, MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD ,ierr)

          IF(ParEnv % MyPE == 0) THEN
            Do i=1,ntot
              Values(DOFs*(PGDOFs_send(i)-1)+k)=PVar_send(i)
            End do
          END IF

        END DO

        ! Create Variable on gathered mesh
        IF(ParEnv % MyPE == 0) THEN
          Do i=1,GatheredMesh % NumberOfNodes
            Perm(i) = i
          End do
          CALL VariableAdd( GatheredMesh % Variables, GatheredMesh, &
              Name=TRIM(MeshSizeName),DOFs=DOFs,Values=Values,Perm=Perm)
          GVar => VariableGet( GatheredMesh % Variables,TRIM(MeshSizeName),ThisOnly=.TRUE.,UnFoundFatal=.TRUE. )
        END IF

        DEALLOCATE(node_count,disps,local_var,PGDOFs_send,PVar_send)
        
      END FUNCTION GATHER_PVar
      
!------------------------------------------------------------------------------
      SUBROUTINE AddValueToMesh(NewMesh,RefMesh)
!------------------------------------------------------------------------------
        implicit none
        TYPE(Solver_t) :: Solver
        TYPE(Mesh_t), POINTER :: NewMesh, RefMesh
        TYPE(Variable_t), POINTER :: Var,NewVar
        TYPE( Matrix_t ), POINTER :: NewMatrix
        INTEGER :: ii,k,jj
        INTEGER, POINTER :: Permutation(:)

        CALL Info('AddVauleToMesh','Enter',Level=1)
        NewMesh % MaxBDOFs = RefMesh % MaxBDOFs
        ! Initialize local variables for the new mesh:
        NULLIFY( NewMesh % Variables )
        CALL VariableAdd( Newmesh % Variables, Newmesh, Solver, &
            'Coordinate 1', 1, NewMesh % Nodes % x )
        CALL VariableAdd( Newmesh % Variables, NewMesh, Solver, &
            'Coordinate 2', 1, Newmesh % Nodes % y )
        CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
            'Coordinate 3', 1, Newmesh % Nodes % z )

        ! Time must always be there:
        ! --------------------------
        Var => VariableGet( RefMesh % Variables,'Time',ThisOnly=.TRUE. )
        CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
            'Time', 1, Var % Values )

        Var => VariableGet( RefMesh % Variables,'Timestep',&
            ThisOnly=.TRUE.)
        CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
            'Timestep', 1, Var % Values ) 

        Var => VariableGet( RefMesh % Variables,'Timestep size',&
            ThisOnly=.TRUE. )
        CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
            'Timestep size', 1, Var % Values )

        Var => VariableGet( RefMesh % Variables,&
            'Timestep interval',ThisOnly=.TRUE. )
        CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
            'Timestep interval', 1, Var % Values )

        Var => VariableGet( RefMesh % Variables,'Coupled iter',&
            ThisOnly=.TRUE. )
        CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
            'Coupled iter', 1, Var % Values )

        Var => VariableGet( RefMesh % Variables,'Nonlin iter',&
            ThisOnly=.TRUE. )
        CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
            'Nonlin iter', 1, Var % Values )

        ! Set Mesh for model
        ! !/r510/home/ltavard/work1/MyElmerIce_WorkProgress/elmerice/fem/src/MeshUtils.F90 
        !-------------------
        CALL SetCurrentMesh(Model,NewMesh)
        
        ! Initialize the field variables for the new mesh. These are
        ! interpolated from the old meshes variables. Vector variables
        ! are in the variable lists in two ways: as vectors and as
        ! vector components.
        ! We MUST update the vectors (i.e. DOFs>1) first!
        ! ------------------------------------------------
        Var => RefMesh % Variables
        DO WHILE( ASSOCIATED( Var ) )
          IF ( Var % DOFs > 1 ) THEN
            NewVar => VariableGet( NewMesh % Variables,Var %Name,.FALSE. )
            NewVar % PrimaryMesh => NewMesh
            ! 
            kk = SIZE( NewVar % Values )
            IF ( ASSOCIATED( NewVar % Perm ) ) THEN
              Kk = COUNT( NewVar % Perm > 0 )
            END IF
            IF ( GetVarName( NewVar ) == 'flow solution' ) THEN
              NewVar % Norm = 0.0d0
              DO ii=1,NewMesh % NumberOfNodes
                DO jj=1,NewVar % DOFs-1
                  NewVar % Norm = NewVar % Norm + &
                      NewVar % Values( NewVar % DOFs*(ii-1)+jj )**2
                END DO
              END DO
              NewVar % Norm = SQRT( NewVar % Norm / kk )
            ELSE
              NewVar % Norm = SQRT( SUM(NewVar % Values**2) / kk )
            END IF
          END IF
          Var => Var % Next
        END DO
        
        !   Second time around, update scalar variables and
        !   vector components:
        !   -----------------------------------------------
        Var => RefMesh % Variables
        DO WHILE( ASSOCIATED( Var ) )
          IF( SIZE( Var % Values ) == Var % DOFs ) THEN
            NewVar => VariableGet( NewMesh % Variables,Var %Name,.TRUE. )
            IF (.NOT.ASSOCIATED(NewVar)) &
                CALL VariableAdd( NewMesh % Variables, NewMesh, Var % Solver, &
                TRIM(Var % Name), Var % DOFs , Var % Values )
            Var => Var % Next
            CYCLE
          END IF
          SELECT CASE( Var % Name )
          CASE('coordinate 1','coordinate 2','coordinate 3','time',&
              'timestep', 'timestep size', 'timestep interval', &
              'coupled iter', 'nonlin iter' )
          CASE DEFAULT
            IF ( Var % DOFs == 1 ) THEN
              Found = .FALSE.
              IF ( Found ) THEN
                kk = Solver % Variable % NameLen
                IF ( Var % Name(1:kk) /= Solver % Variable % Name(1:kk)) THEN
                  Var => Var % Next
                  CYCLE
                END IF
              END IF

              NewVar => VariableGet( NewMesh % Variables, Var % Name,.FALSE. )
              ! added by fab:
              NewVar % PrimaryMesh => NewMesh
              !---------------
              kk = SIZE( NewVar % Values )
              IF ( ASSOCIATED( NewVar % Perm ) ) THEN
                kk = COUNT( NewVar % Perm > 0 )
              END IF
              NewVar % Norm = SQRT( SUM(NewVar % Values**2) / kk )
            END IF
          END SELECT
          Var => Var % Next
        END DO
        
        ! update solver structure to use the new mesh
        !--------------------------------------------
        !Solver % Mesh => NewMesh ! deleted by fab
        CALL MeshStabParams(NewMesh)
        
        ! Nothing computed on this mesh
        !------------------------------          
        NewMesh % SavesDone = 0 ! start new output file -> to check
        NewMesh % OutputActive = .TRUE.
        NewMesh % changed = .TRUE. 
        NewMesh % Next => RefMesh
        
!------------------------------------------------------------------------------
      END SUBROUTINE AddValueToMesh ! Subroutine AddVlueToMesh
!------------------------------------------------------------------------------
 
#else
      CALL Fatal('MMG2DSolver',&
          'Remeshing utility MMG2DSolver has not been installed')
#endif
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END SUBROUTINE MMG2DSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
