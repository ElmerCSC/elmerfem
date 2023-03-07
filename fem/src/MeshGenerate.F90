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
! *  Authors: Juha Ruokolainen, Mikko Lyly
! *  Email:   Juha.Ruokolainen@csc.fi, Mikko.Lyly@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: Autumn 2000
! *
! *****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!--------------------------------------------------------------------------------------------------------
!> Module for determinishting meshing routines without adaptivivity.
!--------------------------------------------------------------------------------------------------------
MODULE MeshGenerate

  USE GeneralUtils
  USE SolverUtils
  USE ModelDescription
  USE LoadMod
  USE MeshUtils
  USE MeshRemeshing

  IMPLICIT NONE

  
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE ReMesh( Model,Solver)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Solver_t), TARGET :: Solver
    TYPE( Model_t ) :: Model
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER   :: RefMesh,NewMesh, Mesh, sMesh
    TYPE( Nodes_t ) :: Nodes
    TYPE( Matrix_t ), POINTER :: NewMatrix
    INTEGER, POINTER :: Permutation(:)
    TYPE( Element_t ), POINTER :: RefElement
    INTEGER :: i,j,k,n,nn,MarkedElements
    TYPE( Variable_t ), POINTER :: HVar, Var, Var1, NewVar
    REAL(KIND=dp) :: MaxChangeFactor, &
        TotalTime,RemeshTime,s,FinalRef,t
    LOGICAL :: BandwidthOptimize, Found, Coarsening, MeshNumbering, DoIt
    INTEGER :: MaxDepth, MinDepth, NLen
    CHARACTER(:), ALLOCATABLE :: Path, VarName
    REAL(KIND=dp), POINTER  :: Time(:), PrevValues(:), &
         Hvalue(:), ptr(:), tt(:)
    REAL(KIND=dp), POINTER  :: eRef(:), hRef(:), Work(:)
    LOGICAL :: NoInterp, Parallel, AdaptInit
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Solver_t), POINTER :: pSolver
    CHARACTER(*), PARAMETER :: Caller = 'ReMesh'

!---------------------------------------------------------------------------------
!
!   Initialize:
!   -----------
    CALL Info( Caller, ' ', Level=5 )
    CALL Info( Caller, &
        '----------- M E S H   R E F I N E M E N T --------------', Level=5 )
    TotalTime = CPUTime()
    RemeshTime = 0.0d0

    RefMesh => Solver % Mesh
    IF( RefMesh % DiscontMesh ) THEN
      CALL Fatal(Caller,'No remeshing possible for discontinuous mesh!')
    END IF

    Params => Solver % Values 
    
    AdaptInit = ( RefMesh % AdaptiveDepth == 0 ) 
    IF( AdaptInit ) THEN
      CALL Info(Caller,'Initializing stuff on coarsest level!')
    END IF
        
    ! Interpolation is costly in parallel. Do it by default only in serial. 
    Parallel = ( ParEnv % PEs > 1 )
    NoInterp = ListGetLogical( Params,'Adaptive Interpolate',Found )
    IF(.NOT. Found) NoInterp = Parallel 
    
!   Compute the local error indicators:
!   -----------------------------------
    t = CPUTime()

!   Add nodal average of the h-value to the mesh variable list:
!   -----------------------------------------------------------

    NN = RefMesh % NumberOfNodes


!   Generate the new mesh:
!   ----------------------
    t = RealTime()
#ifdef HAVE_MMG
    CALL Info(Caller,'Using MMG libary for mesh refinement',Level=5)
    NewMesh => MMG_ReMesh( RefMesh )
#else
    CALL Fatal( Caller,'Remeshing requested with MMG but not compiled with!')
#endif          
    RemeshTime = RealTime() - t
    WRITE( Message, * ) 'Remeshing time (real-secs):                      ',RemeshTime
    CALL Info( Caller, Message, Level=6 )

    IF ( .NOT.ASSOCIATED( NewMesh ) ) THEN
      CALL Info( Caller,'Current mesh seems fine. Nothing to do.', Level=6 )
      RefMesh % OUtputActive = .TRUE.
      RefMesh % Parent % OutputActive = .FALSE.
      GOTO 10
    ELSE
      CALL SetMeshMaxDofs(NewMesh)
    END IF

    CALL Info( Caller,'The new mesh consists of: ', Level=5 )
    CALL Info( Caller,'Nodal points: '&
        //TRIM(I2S(NewMesh % NumberOfNodes)),Level=5)
    CALL Info( Caller,'Bulk elements: '&
        //TRIM(I2S(NewMesh % NumberOfBulkElements)),Level=5)
    CALL Info( Caller,'Boundary elements: '&
        //TRIM(I2S(NewMesh % NumberOfBoundaryElements)),Level=5)

!-------------------------------------------------------------------

!   All the mesh geometry related tables are ready now,
!   next we update model and solver related tables:
!   ----------------------------------------------------

    t = CPUTime()

!   Add the new mesh to the global list of meshes:
!   ----------------------------------------------
    NewMesh % Next   => Model % Meshes 
    Model % Meshes   => NewMesh
    RefMesh % Child  => NewMesh
    NewMesh % Parent => RefMesh
    NewMesh % Child => NULL()

    NewMesh % MaxBDOFs = RefMesh % MaxBDOFs

    NewMesh % Name = ListGetString( Params,'Remesh Mesh Name', Found )
    IF ( .NOT. Found ) NewMesh % Name = 'DeformedMesh'
    
    MeshNumbering = ListGetLogical( Params,'Remesh Mesh Numbering', Found )
    IF(.NOT. Found ) MeshNumbering = .TRUE.
    
    NewMesh % AdaptiveDepth = RefMesh % AdaptiveDepth + 1
    IF( MeshNumbering ) THEN
      NewMesh % Name = TRIM( NewMesh % Name ) // I2S(NewMesh % AdaptiveDepth)
    END IF
    
    IF ( ListGetLogical( Params, 'Remesh Save Mesh', Found ) ) THEN
      Nlen = LEN_TRIM(OutputPath)
      IF ( Nlen > 0 ) THEN
        Path = OutputPath(1:Nlen) // '/' // TRIM(NewMesh % Name)
      ELSE
        Path = TRIM(NewMesh % Name)
      END IF
      CALL MakeDirectory( TRIM(path) // CHAR(0) )
      IF( ParEnv % PEs > 1 ) THEN
        CALL WriteMeshToDisk2( Model, NewMesh, Path, ParEnv % MyPe )
      ELSE
        CALL WriteMeshToDisk( NewMesh, Path )
      END IF
    END IF
    
!   Initialize local variables for the new mesh:
!   --------------------------------------------
    NULLIFY( NewMesh % Variables )
    
    CALL TransferCoordAndTime( RefMesh, NewMesh ) 
    
    ! Initialize the field variables for the new mesh. These are
    ! interpolated from the old meshes variables. Vector variables
    ! are in the variable lists in two ways: as vectors and as
    ! vector components. We MUST update the vectors (i.e. DOFs>1)
    ! first!!!!!
    ! -----------------------------------------------------------
    CALL SetCurrentMesh( Model, NewMesh )

#if 0
    IF( .NOT. ListGetLogical( Params,'Skip Interpolation', Found ) ) THEN
      
      CALL Info(Caller,'Interpolate vectors from old mesh to new mesh!',Level=7)    
      Var => RefMesh % Variables
      DO WHILE( ASSOCIATED( Var ) )
        ! This cycles global variable such as time etc. 
        IF( SIZE( Var % Values ) == Var % DOFs ) THEN
          Var => Var % Next
          CYCLE
        END IF

        IF ( Var % DOFs > 1 ) THEN
          NewVar => VariableGet( NewMesh % Variables,Var % Name,.FALSE. )
          k = SIZE( NewVar % Values )
          IF ( ASSOCIATED( NewVar % Perm ) ) THEN
            k = COUNT( NewVar % Perm > 0 )
          END IF
          IF ( GetVarName( NewVar ) == 'flow solution' ) THEN
            NewVar % Norm = 0.0d0
            DO i=1,NewMesh % NumberOfNodes
              DO j=1,NewVar % DOFs-1
                NewVar % Norm = NewVar % Norm + &
                    NewVar % Values( NewVar % DOFs*(i-1)+j )**2
              END DO
            END DO
            NewVar % Norm = SQRT( NewVar % Norm / k )
          ELSE
            NewVar % Norm = SQRT( SUM(NewVar % Values**2) / k )
          END IF
        END IF
        Var => Var % Next
      END DO
      CALL Info(Caller,'Interpolation to new mesh done!',Level=20)    

      !   Second time around, update scalar variables and
      !   vector components:
      !   -----------------------------------------------
      CALL Info(Caller,'Interpolate scalars from old mesh to new mesh!',Level=7)    
      Var => RefMesh % Variables
      DO WHILE( ASSOCIATED( Var ) )

        ! This cycles global variable such as time etc. 
        IF( SIZE( Var % Values ) == Var % DOFs ) THEN
          Var => Var % Next
          CYCLE
        END IF

        SELECT CASE( Var % Name )

        CASE( 'coordinate 1', 'coordinate 2', 'coordinate 3' )
          CONTINUE

        CASE DEFAULT
          IF ( Var % DOFs == 1 ) THEN
            IF( NoInterp ) THEN
              ! Add field without interpolation
              NewVar => VariableGet( NewMesh % Variables, Var % Name, .TRUE. )
              IF(.NOT. ASSOCIATED(NewVar) ) THEN
                CALL VariableAddVector( NewMesh % Variables, NewMesh, Solver,&
                    Var % Name,Var % DOFs)
                NewVar => VariableGet( NewMesh % Variables, Var % Name, .TRUE. )
              END IF
              NewVar % Norm = Var % Norm
            ELSE
              ! Interpolate scalar variables using automatic internal interpolation 
              NewVar => VariableGet( NewMesh % Variables, Var % Name, .FALSE. )
              k = SIZE(NewVar % Values)
              IF ( ASSOCIATED( NewVar % Perm ) ) THEN
                k = COUNT( NewVar % Perm > 0 )
              END IF
              NewVar % Norm = SQRT(SUM(NewVar % Values**2)/k)
            END IF
          END IF
        END SELECT
        Var => Var % Next
      END DO
    
!-------------------------------------------------------------------    
      WRITE( Message, * ) 'Mesh variable update time (cpu-secs):            ',CPUTime()-t
      CALL Info( Caller, Message, Level = 6 )
!-------------------------------------------------------------------    
    END IF
#endif
    
    
    ! Adaptive meshing is trailing. This is usually preceeding and hance we immediately change the
    ! Mesh to the active one.
    !----------------------------------
    RefMesh % OutputActive = .FALSE.
    NewMesh % SavesDone = 0  
    NewMesh % OutputActive = .TRUE.
    NewMesh % Changed = .TRUE.

!
!   Create matrix structures for the new mesh:
!   ------------------------------------------    
    t = CPUTime()


!   Try to account for the reordering of DOFs
!   due to bandwidth optimization:
!   -----------------------------------------    
    CALL Info( Caller,'Updating primary solver structures: '//TRIM(Solver % Variable % Name))
    CALL UpdateSolverMesh( Solver, NewMesh, NoInterp )          
    CALL ParallelInitMatrix( Solver, Solver % Matrix )

    DO i=1,Model % NumberOfSolvers
      pSolver => Model % Solvers(i)
      IF( .NOT. ASSOCIATED( pSolver ) ) CYCLE
      IF( ASSOCIATED( pSolver, Solver ) ) CYCLE      
      IF( .NOT. ASSOCIATED( pSolver % Mesh, RefMesh ) ) CYCLE

      DoIt = .TRUE.
      IF( .NOT. ASSOCIATED( pSolver % Variable ) ) DoIt = .FALSE.
      IF( DoIt ) THEN
        IF( .NOT. ASSOCIATED( pSolver % Variable % Values ) ) DoIt = .FALSE.
      END IF
      IF( DoIt ) THEN
        IF( SIZE( pSolver % Variable % Values ) == pSolver % Variable % Dofs ) DoIt = .FALSE.
      END IF
      IF( DoIt ) THEN
        CALL Info( Caller,'Updating other solver structures: '//TRIM(pSolver % Variable % Name))
        CALL UpdateSolverMesh( pSolver, NewMesh, NoInterp )          
        CALL ParallelInitMatrix( pSolver, pSolver % Matrix )        
      ELSE
        pSolver % Mesh => NewMesh
      END IF
    END DO

    WRITE( Message, * ) 'Matrix structures update time (cpu-secs):        ',CPUTime()-t
    CALL Info( Caller, Message, Level=6 )

    !   Update Solver structure to use the new mesh:
    !   ---------------------------------------------    
    CALL MeshStabParams( NewMesh )
            
!   Release previous meshes. Keep only the original mesh, and
!   the last two meshes:
!   ---------------------------------------------------------
    n = 0
    Mesh => RefMesh % Parent
    DO WHILE( ASSOCIATED(Mesh) )
      sMesh => Mesh % Parent
      n = n+1
      IF ( Mesh % AdaptiveDepth /= 0 ) THEN
        IF ( ASSOCIATED( Mesh % Parent ) ) THEN
          Mesh % Parent % Child => Mesh % Child                        
        END IF

        IF ( ASSOCIATED(Mesh % Child) ) THEN
          Mesh % Child % Parent => Mesh % Parent
          ! Eliminate the mesh to be released also from here!
          Mesh % Child % Next => Mesh % Next 
        END IF

        CALL Info(Caller,'Releasing mesh: '//TRIM(Mesh % Name),Level=8)
        CALL ReleaseMesh( Mesh )
      END IF

      Mesh => sMesh
    END DO

!------------------------------------------------------------------------------

10  CONTINUE

!   Comment the next calls, if you want to keep the edge tables:
!   ------------------------------------------------------------
    IF(.NOT.isPelement(RefMesh % Elements(1))) THEN
      CALL ReleaseMeshEdgeTables( RefMesh )
      CALL ReleaseMeshFaceTables( RefMesh )
    END IF
    
    CALL SetCurrentMesh( Model, NewMesh )
    Model % Solver % Mesh => NewMesh
    
20  CONTINUE

    WRITE( Message, * ) 'Mesh refine took in total (cpu-secs):           ', &
        CPUTIme() - TotalTime 
    CALL Info( Caller, Message, Level=6 )
    IF ( RemeshTime > 0 ) THEN
      WRITE( Message, * ) 'Remeshing took in total (real-secs):            ', &
          RemeshTime
      CALL Info( Caller, Message, Level=6 )
    END IF
    CALL Info( Caller,'----------- E N D   M E S H   R E F I N E M E N T --------------', Level=5 )
    
    
CONTAINS

  
#ifdef HAVE_MMG

!------------------------------------------------------------------------------
  FUNCTION MMG_ReMesh( RefMesh ) RESULT( NewMesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: NewMesh, RefMesh
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh, TmpMesh
    INTEGER :: i,j,k,n
    CHARACTER(LEN=MAX_NAME_LEN) :: MeshInputFile
    LOGICAL :: Success, Rebalance
    TYPE(Solver_t), POINTER :: pSolver
    LOGICAL :: Visited = .FALSE., UsePerm
!------------------------------------------------------------------------------
    
    Var => VariableGet( RefMesh % Variables, 'Hvalue', ThisOnly=.TRUE. )      

    IF( RefMesh % MeshDim == 2 ) THEN
      CALL Info('MMG_Remesh','Calling serial remeshing routines in 2D',Level=10)
      pSolver => Solver
      UsePerm = ASSOCIATED( Solver % Variable ) 
      IF( UsePerm ) THEN
        UsePerm = ASSOCIATED( Solver % Variable % Perm )
      END IF
      IF( UsePerm ) THEN
        UsePerm = ANY( Solver % Variable % Perm(1:RefMesh % NumberOfNodes) == 0)
      END IF
      IF( UsePerm ) THEN
        CALL Info('MMG_Remesh','Masking meshing where solver is active!')        
        NewMesh => MMG2D_ReMesh( RefMesh, Var, pSolver )
      ELSE
        CALL Info('MMG_Remesh','Performing meshing everywhere!')
        NewMesh => MMG2D_ReMesh( RefMesh, Var )
      END IF
    ELSE
      IF( ParEnv % PEs > 1 ) THEN
        CALL Info('MMG_Remesh','Calling parallel remeshing routines in 3D',Level=10)
        CALL DistributedRemeshParMMG(Model, RefMesh, TmpMesh,&
            Params = Solver % Values, HVar = Var )
        CALL RenumberGElems(TmpMesh)
        Rebalance = ListGetLogical(Model % Solver % Values, "Adaptive Rebalance", Found, DefValue = .TRUE.)
        IF(Rebalance) THEN
          CALL Zoltan_Interface( Model, TmpMesh, StartImbalanceTol=1.1_dp, TolChange=0.02_dp, MinElems=10 )          
          NewMesh => RedistributeMesh(Model, TmpMesh, .TRUE., .FALSE.)
          CALL ReleaseMesh(TmpMesh)
        ELSE
          NewMesh => TmpMesh          
        END IF
      ELSE              
        CALL Info('MMG_Remesh','Calling serial remeshing routines in 3D',Level=10)
        CALL RemeshMMG3D(Model, RefMesh, NewMesh,Params = Solver % Values, &
            HVar = Var, Success = Success )
      END IF
      CALL Info('MMG_Remesh','Finished MMG remeshing',Level=20)      
    END IF

    CALL CheckMeshInfo(NewMesh)
    
    Visited = .TRUE.
    
!------------------------------------------------------------------------------
  END FUNCTION MMG_Remesh
!------------------------------------------------------------------------------
#endif

  
 
!------------------------------------------------------------------------------
END SUBROUTINE ReMesh
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END MODULE MeshGenerate
!-----------------------------------------------------------------------------

!> \} 
