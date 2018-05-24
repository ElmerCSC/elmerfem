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
!!!    TODO: Varaible RELEASE_MESH (default false) allow to release prev
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

      IF(.NOT. ListCheckPresent(SolverParams,'Optimize Bandwidth')) &
        CALL ListAddLogical(SolverParams,'Optimize Bandwidth',.FALSE.)

      END SUBROUTINE MMG2DSolver_init
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE MMG2DSolver( Model,Solver,dt,TransientSimulation)
!------------------------------------------------------------------------------
      USE DefUtils
!------------------------------------------------------------------------------
      IMPLICIT NONE

#ifdef HAVE_MMG
#include "mmg/mmg2d/libmmg2df.h"
#endif
      TYPE(Model_t) :: Model
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Solver_t), POINTER ::PSolver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

#ifdef HAVE_MMG
      MMG5_DATA_PTR_T  :: mmgMesh
      MMG5_DATA_PTR_T  :: mmgSol

      TYPE(Mesh_t),POINTER :: Mesh
      TYPE(Mesh_t),POINTER :: NewMesh,PrevMesh
      TYPE(Variable_t),POINTER :: MeshSize
      TYPE(ValueList_t), POINTER :: SolverParams

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

      nVisit=nVisit+1

      SolverParams => GetSolverParams()

      hsiz = GetConstReal( SolverParams, 'hsiz', UniformSize)
!!
      MeshNumber = MeshNumber + 1

! NewMesh Name
      OutPutFileName = ListGetString(SolverParams,'Output file name',Found)
      IF (.NOT.Found) OutPutFileName = DefaultFileName
      IncrementMeshNumber=&
                 ListGetLogical(SolverParams,'Increment Mesh Number',Found)
      IF (.NOT.Found) IncrementMeshNumber=.TRUE.

! GET REQUIRED VARIABLE
      IF (.NOT.UniformSize) THEN
        MeshSizeName = ListGetString( SolverParams, &
                        'Metric Variable Name',UnFoundFatal=UnFoundFatal )
        MeshSize => VariableGet(Solver%Mesh%Variables,TRIM(MeshSizeName),UnFoundFatal=UnFoundFatal)
        IF (MeshSize%DOFS.NE.1) THEN 
         IF (MeshSize%DOFS.NE.3) &
               CALL FATAL('MMGSolver',&
                       'Variable <ElementSize> should have 1 or 3 DOFs')
        END IF
        Scalar=(MeshSize%DOFS.EQ.1)
      END IF  

!INITIALISATION OF MMG MESH AND SOL STRUCTRURES
      IF (DEBUG) PRINT *,'--**-- INITIALISATION OF MMG '
      mmgMesh = 0
      mmgSol  = 0
      CALL MMG2D_Init_mesh(MMG5_ARG_start, &
             MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
             MMG5_ARG_end)

! SET PARAMETERS 
      IF (DEBUG) PRINT *,'--**-- SET MMG2D PARAMETERS '
      CALL SET_MMG2D_PARAMETERS(SolverParams)

      Mesh => Solver % Mesh

!> COPY MESH TO MMG FORMAT
      IF (DEBUG) PRINT *,'--**-- COPY MESH TO MMG FORMAT'
      CALL SET_MMG2D_MESH(Mesh)

      IF (.NOT.UniformSize) THEN
!> SET THE SOL 
        CALL SET_MMG2D_SOL(Mesh,MeshSize)

! (not mandatory): check if the number of given entities match with mesh size
        CALL MMG2D_Chk_meshData(mmgMesh,mmgSol,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver',&
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
      NewMesh => GET_MMG2D_MESH()

!! MIMIC COMPUTE CHANGE STYLE
      Change=2.*(NewMesh%NumberOfNodes-Mesh%NumberOfNodes)/float(NewMesh%NumberOfNodes+Mesh%NumberOfNodes)
      Change=abs(Change)
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
        IF ((TestConvergence).AND.(Change.LE.Tolerance)) THEN
          CALL ListAddConstReal(Model%Simulation,'Exit Condition',1.0_dp)
        END IF
      END IF

!  SAVE MESH TO DISK
      IF (SAVE_ELMER_MESH) THEN
       Path = TRIM(NewMesh % Name)
       CALL MakeDirectory( TRIM(path) // CHAR(0) )
       CALL WriteMeshToDisk( NewMesh, Path )
      END IF

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

      FUNCTION GET_MMG2D_MESH() RESULT(NewMesh)
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Mesh_t), POINTER :: NewMesh
!------------------------------------------------------------------------------
      TYPE(Element_t),POINTER ::  Element
      INTEGER, POINTER :: NodeIndexes(:)
      INTEGER :: np,nt,na,ier
      INTEGER :: ref,corner,required,ridge
      INTEGER :: parent,ied
      INTEGER :: tt, jj, kk, ll


     !> a) get the size of the mesh: vertices,  triangles, edges
      CALL MMG2D_Get_meshSize(mmgMesh,np,nt,na,ier)
      IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                           'CALL TO MMGS_Get_meshSize FAILED')
      IF (DEBUG) PRINT *,'--**-- MMG2D_Get_meshSize DONE'    

! INITIALISE THE NEW MESH STRUCTURE
      NewMesh => AllocateMesh()
      IF (IncrementMeshNumber) THEN
         write(NewMesh % Name,'(A,A,I0)') TRIM(OutPutFileName),'_N',MeshNumber
      ELSE
         NewMesh % Name=TRIM(OutPutFileName)
      END IF

      NewMesh%MaxElementNodes=3
      NewMesh%MaxElementDOFs=3
      NewMesh%MeshDim=Solver%Mesh%MeshDim

      NewMesh % NumberOfNodes = np
      NewMesh % NumberOfBulkElements = nt
      NewMesh % NumberOfBoundaryElements = na

      CALL AllocateVector( NewMesh % Nodes % x, NewMesh % NumberOfNodes)
      CALL AllocateVector( NewMesh % Nodes % y, NewMesh % NumberOfNodes)
      CALL AllocateVector( NewMesh % Nodes % z, NewMesh % NumberOfNodes)
      CALL AllocateVector( NewMesh % Elements, NewMesh % NumberOfBulkElements + &
               NewMesh % NumberOfBoundaryElements )
      
!! GET NEW VERTICES
      NewMesh % Nodes % z = 0._dp
      Do ii=1,np
         CALL MMG2D_Get_vertex(mmgMesh,&
                                       NewMesh % Nodes % x(ii),&
                                       NewMesh % Nodes % y(ii),&
                                       ref,corner,required,ier)
         IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                           'CALL TO  MMG2D_Get_vertex FAILED')
      End do
      IF (DEBUG) PRINT *,'MMG2D_Get_vertex DONE'

!! GET NEW TRIANGLES
      Do tt=1,NewMesh % NumberOfBulkElements
         Element => NewMesh % Elements(tt)
         Element % TYPE => GetElementType( 303 )
         Element % NDOFs = Element % TYPE % NumberOfNodes
         Element % ElementIndex = tt
         Element % PartIndex = ParEnv % myPE
         CALL AllocateVector(Element % NodeIndexes, 3)
         NodeIndexes => Element % NodeIndexes
         CALL MMG2D_Get_triangle(mmgMesh, &
                                         NodeIndexes(1), &
                                         NodeIndexes(2), &
                                         NodeIndexes(3), &
                                         Element % BodyId, &
                                         required,ier)
      End do
      IF (DEBUG) PRINT *,'MMG2D_Get_triangle DONE'
      
!! Get BC Elements
      kk=NewMesh % NumberOfBulkElements
      Do ii=1,na
         kk = kk + 1

         Element => NewMesh % Elements(kk)

         Element % TYPE => GetElementType( 202 )
         Element % NDOFs = Element % TYPE % NumberOfNodes
         Element % ElementIndex = kk
         Element % PartIndex = ParEnv % myPE

         CALL AllocateVector(Element % NodeIndexes, 2)
         NodeIndexes => Element % NodeIndexes
         CALL MMG2D_Get_edge(mmgMesh, &
                                     NodeIndexes(1), &
                                     NodeIndexes(2), &
                                     ref   , &
                                     ridge,required,ier)
         IF ( ier /= 1 ) CALL FATAL('MMGSolver',&
                          'CALL TO  MMG2D_Get_Edge FAILED')
        
         Allocate(Element % BoundaryInfo)
         Element % BoundaryInfo % Constraint=ref

      End Do

      kk=NewMesh % NumberOfBulkElements
      Do ii=1,na
         kk = kk + 1

         Element => NewMesh % Elements(kk)

         CALL MMG2D_GET_TRIFROMEDGE(mmgMesh, &
                                    ii,parent,ied,ier)

         IF ( ier /= 1 ) CALL FATAL('MMGSolver',&
                          'CALL TO  MMG2D_Get_TRIFROMEDGE FAILED')

         Element % BoundaryInfo % Left => NewMesh % Elements(parent)
      End do

      END FUNCTION GET_MMG2D_MESH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_MMG2D_SOL(Mesh,MeshSize)
      IMPLICIT NONE
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Variable_t), POINTER :: MeshSize
    
      REAL(KIND=dp) :: M11,M22,M12

      INTEGER :: NVert
      INTEGER :: ier
      INTEGER :: ii

      NVert= Mesh%NumberOfNodes

      IF (Scalar) THEN
         CALL MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NVert,MMG5_Scalar,ier)
      ELSE
         CALL MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NVert,MMG5_Tensor,ier)
      END IF
      IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                           'CALL TO MMG2D_Set_SolSize FAILED')
      IF (DEBUG) PRINT *,'--**-- MMG2D_Set_solSize DONE'

      Do ii=1,NVert
        IF (Scalar) THEN
          CALL MMG2D_Set_scalarSol(mmgSol,&
                     MeshSize%Values(MeshSize%Perm(ii)), &
                     ii,ier)
        ELSE
          M11=MeshSize%Values(3*(MeshSize%Perm(ii)-1)+1)
          M22=MeshSize%Values(3*(MeshSize%Perm(ii)-1)+2)
          M12=MeshSize%Values(3*(MeshSize%Perm(ii)-1)+3)
          CALL MMG2D_Set_tensorSol(mmgSol,&
                     M11,M12,M22,ii,ier)
        ENDIF
        IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                             'CALL TO MMG2D_Set_scalarSo FAILED')
      End Do
      IF (DEBUG) PRINT *,'--**-- MMG2D_Set_tensorSol DONE'
      END SUBROUTINE SET_MMG2D_SOL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_MMG2D_MESH(Mesh)
      IMPLICIT NONE
      TYPE(Mesh_t), POINTER :: Mesh

      TYPE(Element_t),POINTER :: Element
      INTEGER, POINTER :: NodeIndexes(:)

      INTEGER :: NVert,NEle,NEdge
      INTEGER :: n
      INTEGER :: ier
      INTEGER :: ii,tt
      
      NVert=Mesh%NumberOfNodes
      NEle=Mesh%NumberOfBulkElements
      NEdge=Mesh%NumberOfBoundaryElements

      CALL MMG2D_Set_meshSize(mmgMesh,NVert,NEle,NEdge,ier)
      IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                        'CALL TO MMG2D_Set_meshSize FAILED')
      IF (DEBUG) PRINT *,'--**-- MMG2D_Set_meshSize DONE'

      Do ii=1,NVert
        CALL MMG2D_Set_vertex(mmgMesh, Mesh%Nodes%x(ii), &
                                       Mesh%Nodes%y(ii), &
                                       0, ii , ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                        'CALL TO MMG2D_Set_vertex FAILED')
      End Do
      IF (DEBUG) PRINT *,'--**-- MMG2D_Set_vertex DONE'

      Do tt=1,NEle
         Element => Solver % Mesh % Elements(tt)
         NodeIndexes => Element % NodeIndexes

         IF (Element % TYPE % ElementCode .NE. 303) &
             CALL FATAL('MMGSolver','Work only with 303 elements')
         n = GetElementNOFNodes(Element)

         CALL MMG2D_Set_triangle(mmgMesh, NodeIndexes(1), &
                                          NodeIndexes(2), &
                                          NodeIndexes(3), &
                                          Element % BodyId, tt, ier)
         IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                             'CALL TO MMG2D_Set_triangle FAILED')
      End do
      IF (DEBUG) PRINT *,'--**-- MMG2D_Set_triangle DONE'

      Do tt=1,NEdge
         Element => GetBoundaryElement(tt)
         NodeIndexes => Element % NodeIndexes

         IF (Element % TYPE % ElementCode .NE. 202) &
                 CALL FATAL('MMGSolver',&
                      'Work only with 202 boundary elements')
         n = GetElementNOFNodes(Element)

         CALL MMG2D_Set_Edge(mmgMesh, NodeIndexes(1), &
                                      NodeIndexes(2), &
                                      Element % BoundaryInfo % Constraint, &
                                      tt, ier)
         IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                              'CALL TO MMG2D_Set_Edge FAILED')
      End Do
      IF (DEBUG) PRINT *,'--**-- MMG2D_Set_Edge DONE'

      END SUBROUTINE SET_MMG2D_MESH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_MMG2D_PARAMETERS(SolverParams)
      IMPLICIT NONE
      TYPE(ValueList_t), POINTER :: SolverParams

      REAL(KIND=dp) :: hsiz,Pval
      INTEGER :: ier
      LOGICAL :: AngleDetect
      INTEGER :: verbosity,MeMIncrease,Bucket,GMSHoption     
      LOGICAL :: DebugMode,NoInsert,NoSwap,NoMove,NoSurf
      LOGICAL :: Found

      ! Minimal mesh size:  hmin
      Hmin = GetConstReal( SolverParams, 'hmin', Found)
      IF (Found) THEN
        CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_hmin,&
                 Hmin,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver', &
                'CALL TO MMG2D_SET_DPARAMETER <hmin> Failed')
      END IF
      ! Maximal mesh size - hmax
      Hmax = GetConstReal( SolverParams, 'hmax', Found)
      IF (Found) THEN
        CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_hmax,&
                Hmax,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver', &
                'CALL TO MMG2D_SET_DPARAMETER <hmax> Failed')
      END IF

      hsiz = GetConstReal( SolverParams, 'hsiz', Found)
      IF (Found) THEN
         CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_hsiz,&
                 hsiz,ier)
         IF ( ier == 0 ) CALL FATAL('MMGSolver', &
                'CALL TO MMG2D_SET_DPARAMETER <hsiz> Failed')
      END IF

!!! PARAMS: generic options (debug, mem, verbosity)
        ! [val] Set the verbosity level to n
      Verbosity=GetInteger( SolverParams,'verbosity',Found)
      IF (Found) THEN
        CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_verbose, &
               Verbosity,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                               'CALL TO MMG2D_SET_IPARAMETER Failed')
      END IF
        ! [val] Set the maximal memory size to n MBytes.
      MemIncrease=GetInteger(SolverParams,'Increase Memory',Found)
      IF (FOUND) THEN
        CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_mem,&
                MemIncrease,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver', &
                               'CALL TO MMG2D_SET_IPARAMETER Failed')
      END IF
        ! [0/1] Turn on the debug mode.
      DebugMode=GetLogical(SolverParams,'Debug Mode',Found)
      IF (Found .AND. DebugMode) THEN
        CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_debug,1, &
               ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                               'CALL TO MMG2D_SET_IPARAMETER Failed')
      END IF
!!! PARAMS
      ! Control global Hausdorff distance (on all the boundary surfaces of the mesh)
      ! MMG2D_DPARAM_hausd default est 0.01 semble bien trop petit;
      !  il semble qu'il faille mettre une taille > taille des elements.
      Pval=GetConstReal( SolverParams, 'hausd', Found)
      IF (Found) THEN
         CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_hausd,&
                 Pval,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver', &
                'CALL TO MMG2D_SET_DPARAMETER <hausd> Failed')
      END IF
      ! Control gradation
      Pval=GetConstReal( SolverParams, 'hgrad', Found)
      IF (Found) THEN
         CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_hgrad,&
                 Pval,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver', &
                'CALL TO MMG2D_SET_DPARAMETER <hgrad> Failed')
      END IF

!!! OTHER PARAMETERS: NOT ALL TESTED
      Pval = GetConstReal( SolverParams, 'Angle detection',Found)
      IF (Found) THEN
        CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_angleDetection,&
                 Pval,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver', &
                'CALL TO MMG2D_SET_DPARAMETER <Angle detection> Failed')
      ENDIF
      ! !< [1/0], Avoid/allow surface modifications */ 
      AngleDetect=GetLogical(SolverParams,'No Angle detection',Found)
      IF (Found.AND.AngleDetect) THEN
        CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_angle, &
                1,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver', &
              'CALL TO MMG2D_SET_IPARAMETER <No Angle detection> Failed')
      END IF
      ! [1/0] Avoid/allow point insertion
      NoInsert=GetLogical(SolverParams,'No insert',Found)
      IF (Found .AND. NoInsert) THEN
        CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_noinsert,&
                1,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver', &
                      'CALL TO MMG2D_SET_IPARAMETER <No insert> Failed') 
      END IF
      ! [1/0] Avoid/allow edge or face flipping
      NoSwap=GetLogical(SolverParams,'No swap',Found)
      IF (Found .AND. NoSwap) THEN
        CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_noswap,&
                 1,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                        'CALL TO MMG2D_SET_IPARAMETER <No swap>Failed')
      END IF
      ! [1/0] Avoid/allow point relocation
      NoMove=GetLogical(SolverParams,'No move',Found)
      IF (Found .AND. NoMove) THEN
        CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_nomove,&
                1,ier)
        IF ( ier == 0 ) CALL FATAL('MMGSolver',&
                     'CALL TO MMG2D_SET_IPARAMETER <No move> Failed')
      END IF
      ! [1/0] Avoid/allow surface modifications
      NoSurf=GetLogical(SolverParams,'No surf',Found)
      IF (Found .AND. NoSurf) THEN
        CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_nosurf,&
                0,ier)
        IF ( ier == 0 ) CALL FATAL('M/MGSolver',&
                       'CALL TO MMG2D_SET_IPARAMETER <No surf> Failed')
      END IF
!!!
      END SUBROUTINE SET_MMG2D_PARAMETERS

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
          LOGICAL :: BandwidthOptimize, GlobalBubbles 
 
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
        End ! Subroutine AddVlueToMesh
!------------------------------------------------------------------------------
 
#else
     CALL FATAL('MMG2DSolver',&
        'Remeshing utility MMG2DSolver has not been installed')
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END SUBROUTINE MMG2DSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
