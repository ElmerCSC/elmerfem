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
    LOGICAL :: NoInterp, Parallel
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Solver_t), POINTER :: pSolver
    INTEGER :: VisitedCount = 0, RemeshInterval    
    CHARACTER(*), PARAMETER :: Caller = 'ReMesh'

    INTERFACE
      SUBROUTINE InterpolateMeshToMesh( OldMesh, NewMesh, OldVariables, &
          NewVariables, UseQuadrantTree, Projector, MaskName, UnfoundNodes )
        USE Types
        TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
        TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
        LOGICAL, OPTIONAL :: UseQuadrantTree
        LOGICAL, POINTER, OPTIONAL :: UnfoundNodes(:)
        CHARACTER(LEN=*),OPTIONAL :: MaskName
        TYPE(Projector_t), POINTER, OPTIONAL :: Projector
      END SUBROUTINE InterpolateMeshToMesh
    END INTERFACE
    
    SAVE VisitedCount
    
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
    VisitedCount = VisitedCount + 1

    RemeshInterval = ListGetInteger( Params,'Remesh Interval',Found )
    IF( Found ) THEN
      IF( MODULO( VisitedCount, RemeshInterval ) == 0 ) THEN
        CALL Info( Caller,'Visited Count is '//I2S(VisitedCount)//', remeshing is active!')
      ELSE
        CALL Info( Caller,'Visited Count is '//I2S(VisitedCount)//', doing nothing!')
        RETURN
      END IF
    END IF         
    
    ! Interpolation is best done all in one sweep at the end
    ! Hence this is omitted at the moment. 
    NoInterp = .TRUE. 
    
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
    
!   Initialize glocal variables for the new mesh:
!   --------------------------------------------
    NULLIFY( NewMesh % Variables )    
    CALL TransferCoordAndTime( RefMesh, NewMesh ) 
    
    ! Set the current mesh.
    ! -----------------------------------------------------------
    CALL SetCurrentMesh( Model, NewMesh )

    ! Change the active mesh to be the new one.
    !------------------------------------------
    RefMesh % OutputActive = .FALSE.
    NewMesh % SavesDone = 0  
    NewMesh % OutputActive = .TRUE.
    NewMesh % Changed = .TRUE.
    NewMesh % Projector => NULL()
    
!   Create matrix structures for the new mesh:
!   ------------------------------------------    
    t = CPUTime()
    CALL Info( Caller,'Updating primary solver structures: '//TRIM(Solver % Variable % Name))
    CALL UpdateSolverMesh( Solver, NewMesh, NoInterp )
    NewMesh % Changed = .TRUE.
    CALL ParallelInitMatrix( Solver, Solver % Matrix )

    ! This is no longer valid. We could call the routine to update this
    ! but it is not visible in this module now so just deallocate it.
    IF(Solver % NumberOfActiveElements > 0 ) THEN
      DEALLOCATE( Solver % ActiveElements )
      Solver % NumberOfActiveElements = 0
    END IF

    ! Update also other variables that use the same mesh but in different solvers.    
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
        IF(pSolver % NumberOfActiveElements > 0 ) THEN
          DEALLOCATE(pSolver % ActiveElements )
          pSolver % NumberOfActiveElements = 0
        END IF
      ELSE
        pSolver % Mesh => NewMesh
      END IF
    END DO

    WRITE( Message, * ) 'Matrix structures update time (cpu-secs):        ',CPUTime()-t
    CALL Info( Caller, Message, Level=6 )

    !   Update Solver structure to use the new mesh:
    !   ---------------------------------------------    
    CALL MeshStabParams( NewMesh )

    Model % Solver % Mesh => NewMesh
        
    ! Here is all the real interpolation work (serial & parallel)
    ! Do it at the end since now all variables have been properly initialized
    ! and mainly interpolation remains to be done. 
    !---------------------------------------------------------------------
    IF( .NOT. ListGetLogical( Params,'Skip Interpolation', Found ) ) THEN
      CALL Info(Caller,'Mapping all fields in old mesh to new mesh!',Level=7)
      CALL InterpolateMeshToMesh( RefMesh, NewMesh, &
          RefMesh % Variables, NewMesh % Variables )
      CALL Info(Caller,'Done mapping all fields in old mesh to new mesh!',Level=15)

      ! Compute the norms so that we have better estimate for the accuracy
      DO i=1,CurrentModel % NumberOfSolvers
        pSolver => CurrentModel % Solvers(i)
        IF(.NOT. ASSOCIATED( NewMesh, pSolver % Mesh ) ) CYCLE
        IF(.NOT. ASSOCIATED( pSolver % Variable) ) CYCLE
        IF(.NOT. ASSOCIATED( pSolver % Matrix ) ) CYCLE
        CALL Info(Caller,'Computing norm for mapped field of solver: '//I2S(pSolver % SolverId),Level=12)
        pSolver % Variable % Norm = ComputeNorm(pSolver, &
            SIZE(pSolver % Variable % Values), pSolver % Variable % Values ) 
      END DO
    END IF    

    IF( ListGetLogical( Params,'Elemental Interpolation', Found ) ) THEN
      IF( ASSOCIATED( NewMesh % InvPerm ) ) THEN
        CALL InterpolateMeshToMeshElemental(RefMesh, NewMesh, NewMesh % InvPerm)
        DEALLOCATE( NewMesh % InvPerm )
      END IF
    END IF

    
    ! Release previous meshes. Keep only the original mesh, and
    ! the last two meshes:
    !---------------------------------------------------------
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
    
20  CONTINUE

    WRITE( Message, * ) 'Mesh alteration took in total (cpu-secs):           ', &
        CPUTIme() - TotalTime 
    CALL Info( Caller, Message, Level=6 )
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
    CHARACTER(LEN=MAX_NAME_LEN) :: Hname
    LOGICAL :: Success, Rebalance
    TYPE(Solver_t), POINTER :: pSolver
    TYPE(Variable_t), POINTER :: hVar
    LOGICAL :: Visited = .FALSE., UsePerm
!------------------------------------------------------------------------------
    
    Var => VariableGet( RefMesh % Variables, 'Hvalue', ThisOnly=.TRUE. )      

    hName = ListGetString( Params, 'Metric Variable Name',Found )
    IF(.NOT. Found) hName = "hvalue"
    hVar => VariableGet(RefMesh % Variables,hname)
    IF( ASSOCIATED( hVar ) ) THEN
      IF (hvar % dofs /= 1 .AND. hvar % dofs /= 3) THEN
        CALL Fatal('MMG_ReMesh','Variable "'//TRIM(hname)//'" should have 1 or 3 DOFs!')
      END IF
      CALL Info('MMG_Remesh','Using externally provided mesh metric',Level=10)
    END IF

    UsePerm = .FALSE.
    IF( ListGetLogical( Params,'Remesh Active Regions', Found ) ) THEN      
      pSolver => Solver
      UsePerm = ASSOCIATED( Solver % Variable ) 
      IF( UsePerm ) THEN
        IF(.NOT. ASSOCIATED( Solver % Variable % Perm ) ) THEN
          CALL Fatal('MMG_Remesh','Requesting "Remesh Active Regions" but Perm not associated!')
        END IF
      END IF
      IF( UsePerm ) THEN
        UsePerm = ANY( Solver % Variable % Perm(1:RefMesh % NumberOfNodes) == 0)
      END IF
    END IF
    IF( UsePerm ) THEN
      CALL Info('MMG_Remesh','Remeshing only active regions')
    END IF

    IF( RefMesh % MeshDim == 2 ) THEN
      CALL Info('MMG_Remesh','Calling serial remeshing routines in 2D',Level=10)
      IF( UsePerm ) THEN
        CALL Info('MMG_Remesh','Masking meshing where solver is active!')        
        NewMesh => MMG2D_ReMesh( RefMesh, hVar, pSolver )
      ELSE
        CALL Info('MMG_Remesh','Performing meshing everywhere!')
        NewMesh => MMG2D_ReMesh( RefMesh, hVar )
      END IF
    ELSE
      IF( ParEnv % PEs > 1 ) THEN
        CALL Info('MMG_Remesh','Calling parallel remeshing routines in 3D',Level=10)
        CALL DistributedRemeshParMMG(Model, RefMesh, TmpMesh,&
            Params = Solver % Values, HVar = hVar )
        CALL RenumberGElems(TmpMesh)
        Rebalance = ListGetLogical(Params, 'Adaptive Rebalance', Found, DefValue = .TRUE.)
        IF(Rebalance) THEN
          CALL Zoltan_Interface( Model, TmpMesh, StartImbalanceTol=1.1_dp, TolChange=0.02_dp, MinElems=10 )          
          NewMesh => RedistributeMesh(Model, TmpMesh, .TRUE., .FALSE.)
          CALL ReleaseMesh(TmpMesh)
        ELSE
          NewMesh => TmpMesh          
        END IF
      ELSE              
        CALL Info('MMG_Remesh','Calling serial remeshing routines in 3D',Level=10)
        IF( UsePerm ) THEN
          CALL RemeshMMG3D(Model, RefMesh, NewMesh,Params = Params, &
              HVar = hVar, Solver = pSolver, Success = Success )
        ELSE
          CALL RemeshMMG3D(Model, RefMesh, NewMesh,Params = Params, &
              HVar = hVar, Success = Success )
        END IF
      END IF
      CALL Info('MMG_Remesh','Finished MMG remeshing',Level=20)      
    END IF

    !For debugging
    !CALL CheckMeshInfo(NewMesh)
    
    Visited = .TRUE.
    
!------------------------------------------------------------------------------
  END FUNCTION MMG_Remesh
!------------------------------------------------------------------------------
#endif

!------------------------------------------------------------------------------
!> Interpolate the values that are in the elements that stay the same. Then
!> the interpolation becomes easy involving only mapping of the indexes of old
!> mesh to the new mesh.
!------------------------------------------------------------------------------
  SUBROUTINE InterpolateMeshToMeshElemental( OldMesh, NewMesh, InvPerm )
!------------------------------------------------------------------------------
       TYPE(Mesh_t), TARGET  :: OldMesh   !< Old mesh structure
       TYPE(Mesh_t), TARGET  :: NewMesh   !< New mesh structure
       INTEGER, POINTER :: InvPerm(:)
       !------------------------------------------------------------------------------
       TYPE(Element_t), POINTER :: Element, OldElement
       TYPE(Variable_t), POINTER :: Var, OldVar
       INTEGER :: i,j,k,ii,jj,kk,m,np, oldnp, dofs
       INTEGER, POINTER :: pIndexes(:), OldpIndexes(:)
       TYPE(Solver_t), POINTER :: pSolver
       LOGICAL :: DoIt
!------------------------------------------------------------------------------
       
       IF ( OldMesh % NumberOfNodes == 0 ) RETURN       
       IF ( OldMesh % NumberOfBulkElements == 0 ) RETURN       
             
!------------------------------------------------------------------------------
! Loop over all bulk elements in the new mesh
!------------------------------------------------------------------------------
       DO i=1,NewMesh % NumberOfBulkElements
         
         j = InvPerm(i)
         IF(j==0) CYCLE
         Element => NewMesh % Elements(i)         
         
         ! Find the same element in the old mesh         
         OldElement => OldMesh % Elements(j)

         ! Go through all variables to be interpolated:
         Var => NewMesh % Variables
         DO WHILE( ASSOCIATED( Var ) )            
           DoIt = .TRUE.
           
           ! There are many reasons to skip interpolation of this field
           IF( SIZE( Var % Values ) == Var % DOFs ) DoIt = .FALSE.            
           IF( Var % Secondary ) DoIt = .FALSE.            
           IF( Var % Name(1:10) == 'coordinate') DoIt = .FALSE.
           IF( Var % Secondary ) DoIt = .FALSE.
           
           IF( DoIt ) THEN
             ! Get the old variable with the same name
             dofs = Var % Dofs 
             OldVar => VariableGet( OldMesh % Variables, Var % Name, .TRUE. ) 
             m = 0
             IF( ASSOCIATED( OldVar ) ) THEN
               IF( ASSOCIATED( OldVar % PrevValues ) ) THEN
                 m = SIZE(OldVar % PrevValues, 2)
                 IF(.NOT. ASSOCIATED( Var % PrevValues ) ) THEN
                   ALLOCATE( Var % PrevValues(SIZE(Var % Values),m))
                   Var % PrevValues = 0.0_dp
                 END IF
               END IF
             END IF
               
             IF(.NOT. ASSOCIATED(OldVar) ) THEN
               CONTINUE
               
             ELSE IF( Var % TYPE == Variable_on_nodes ) THEN
               DO ii=1,Element % TYPE % NumberOfNodes
                 jj = Var % Perm( Element % NodeIndexes(ii) )
                 IF(jj==0) CYCLE
                 kk = OldVar % Perm( OldElement % NodeIndexes(ii) )
                 IF(kk==0) CYCLE
                 
                 Var % Values( dofs*(jj-1)+1:dofs*jj ) = OldVar % Values(dofs*(kk-1)+1:dofs*kk )
                 IF(m>0) Var % PrevValues( dofs*(jj-1)+1:dofs*jj,1:m ) = &
                     OldVar % PrevValues( dofs*(kk-1)+1:dofs*kk,1:m )
               END DO
                 
             ELSE IF (Var % TYPE == Variable_on_nodes_on_elements ) THEN
               DO ii=1,Element % TYPE % NumberOfNodes
                 jj = Var % Perm( Element % DGIndexes(ii) )
                 IF(jj==0) CYCLE
                 kk = OldVar % Perm( OldElement % DGIndexes(ii) )
                 IF(kk==0) CYCLE
                 
                 Var % Values( dofs*(jj-1)+1:dofs*jj ) = OldVar % Values( dofs*(kk-1)+1:dofs*kk )
                 IF(m>0) Var % PrevValues( dofs*(jj-1)+1:dofs*jj,1:m ) = &
                     OldVar % PrevValues( dofs*(kk-1)+1:dofs*kk,1:m )
               END DO
                 
             ELSE IF( Var % Type == Variable_on_elements ) THEN
               jj = Var % Perm( i )
               IF(jj==0) CYCLE
               kk = OldVar % Perm( j )
               IF(kk==0) CYCLE
               
               Var % Values( dofs*(jj-1)+1:dofs*jj ) = OldVar % Values( dofs*(kk-1)+1:dofs*kk ) 
               IF(m>0) Var % PrevValues( dofs*(jj-1)+1:dofs*jj,1:m ) = &
                   OldVar % PrevValues( dofs*(kk-1)+1:dofs*kk,1:m )
             ELSE
               pSolver => Var % Solver
               IF( ASSOCIATED( pSolver ) ) THEN
                 np = mGetElementDOFs(pIndexes,Element,USolver=pSolver,UMesh=Mesh)
                 oldnp = mGetElementDOFs(OldpIndexes,OldElement,USolver=pSolver,UMesh=Mesh)
                 
                 DO ii=1,np
                   jj = Var % Perm(pIndexes(ii))
                   IF(jj==0) CYCLE
                   kk = OldVar % perm(OldPIndexes(ii))
                   IF(kk==0) CYCLE
                   
                   Var % Values( dofs*(jj-1)+1:dofs*jj ) = OldVar % Values( dofs*(kk-1)+1:dofs*kk ) 
                   IF(m>0) Var % PrevValues( dofs*(jj-1)+1:dofs*jj,1:m ) = &
                       OldVar % PrevValues( dofs*(kk-1)+1:dofs*kk,1:m )
                 END DO
               END IF
             END IF
           END IF

           Var => Var % next
         END DO
       END DO
       
     END SUBROUTINE InterpolateMeshToMeshElemental

!------------------------------------------------------------------------------
   END SUBROUTINE ReMesh
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END MODULE MeshGenerate
!-----------------------------------------------------------------------------

!> \} 
