!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!---------------------------------------------------------------
!> Subroutine for solving either gebhardt factors or radiosities.
!> The radiation may be resolved on a coarser level of the mesh
!> hierarchy, hence restore the current mesh when done.
!---------------------------------------------------------------
   SUBROUTINE RadiationFactors( TSolver, TopoCall, Newton )

     USE DefUtils
     IMPLICIT NONE

     LOGICAL :: TopoCall
     LOGICAL :: Newton
     TYPE(Solver_t) :: TSolver

     CALL RadiationFactorsMesh( TSolver, TopoCall, Newton )
     CALL SetCurrentMesh( CurrentModel, TSolver % Mesh )

   END SUBROUTINE RadiationFactors


   SUBROUTINE RadiationFactorsMesh( TSolver, TopoCall, Newton )

     USE DefUtils
     IMPLICIT NONE

     LOGICAL :: TopoCall
     LOGICAL :: Newton
     TYPE(Solver_t) :: TSolver

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(Mesh_t), POINTER :: Mesh, FineMesh

     ! Radiosity may be computed on a coarser level of the "Mesh Levels"
     ! hierarchy (Hybrid). Then the temperatures and emissivities are
     ! aggregated from the finer (heat equation) mesh boundary elements
     ! and the irradiation is returned to them.
     ! The radiation mesh may also be an independent mesh (GeneralMesh), given by
     ! "Radiation Mesh". Then the mapping between the meshes is computed geometrically.
     LOGICAL :: Hybrid, GeneralMesh
     INTEGER :: RadLevel, nFine
     TYPE(Mesh_t), POINTER :: GeneralRadMesh => NULL()
     CHARACTER(:), ALLOCATABLE :: RadMeshName

     ! The independent radiation mesh follows the motion of the heat equation mesh:
     ! reference (undeformed) coordinates of both meshes, and the fine element with
     ! local coordinates where each coarse node lies in the reference geometry.
     REAL(KIND=dp), ALLOCATABLE :: FineRef(:,:), CoarseRef(:,:), AnchorUV(:,:)
     INTEGER, ALLOCATABLE :: AnchorElem(:)
     LOGICAL :: MotionReported = .FALSE.
     INTEGER, ALLOCATABLE :: FineParent(:), FineElem(:)
     CHARACTER(:), ALLOCATABLE :: MeshDirName, RadDirName, VFSuffix, RadSuffix
     REAL(KIND=dp), ALLOCATABLE :: FineEmis(:), FineAbs(:), FineT(:), CoarseT(:)
     REAL(KIND=dp), ALLOCATABLE :: FineArea(:), FineGrad(:), CoarseA(:), CoarseS(:), FineW(:)

     ! Mapping between the fine and coarse boundary elements (the matrix W): fine element
     ! MapFine(p) covers the coarse element MapCoarse(p) by fraction MapW(p) of its area.
     ! For the current radiation body the same is given in terms of the fine data index
     ! (PcFine) and coarse radiation surface index (PcRad).
     INTEGER :: nMap = 0, nPc
     INTEGER, ALLOCATABLE :: MapFine(:), MapCoarse(:), PcFine(:), PcRad(:)
     REAL(KIND=dp), ALLOCATABLE :: MapW(:), PcW(:), CoarseScale(:)

     TYPE(Factors_t), POINTER :: ViewFactors(:)
     TYPE(Matrix_t),  POINTER :: G => Null()
     TYPE(Element_t), POINTER :: Element
     TYPE(Solver_t),  POINTER :: Solver => Null()
          
     INTEGER :: i,istat,nBndr,nBulk

     REAL (KIND=dp), ALLOCATABLE :: Reflectivity(:),Emissivity(:),Absorptivity(:), &
         Areas(:), RelAreas(:)
     REAL (KIND=dp) :: at, bt,  st, rt
     REAL (KIND=dp) :: SteadyChange, Tol, Sigma

     INTEGER :: RadiationSurfaces, GeometryFixedAfter, TimesVisited=0, RadiationBody, &
                MaxRadiationBody, FactorsFixedAfter

     CHARACTER(:), ALLOCATABLE :: RadiationFlag, SolverType
     INTEGER, ALLOCATABLE :: ElementNumbers(:), InvElementNumbers(:)

     LOGICAL :: SaveFactors, UpdateViewFactors, UpdateGebhartFactors,         &
         ComputeViewFactors, TopologyTest, TopologyFixed, FullMatrix,         &
         IterSolveFactors, ConstantEmissivity, Found,  UpdateRadiatorFactors, &
         ComputeRadiatorFactors, RadiatorsFound, DiffuseGrayRadiationFound,   &
         UpdateGeometry, Radiosity, Spectral
     LOGICAL :: FirstTime = .TRUE.

     INTEGER, PARAMETER :: VFUnit = 10
     LOGICAL, ALLOCATABLE :: ActiveNodes(:), ActiveMe(:), ActiveTasks(:)
     TYPE(ValueList_t), POINTER :: Params, BC

     LOGICAL :: UseFullMatrix
     REAL(KIND=dp), ALLOCATABLE, TARGET, SAVE :: G_Full(:,:)

     INTEGER(KIND=AddrInt) :: mvProc, AddrFunc
     EXTERNAL AddrFunc
     CHARACTER(*), PARAMETER :: Caller = 'RadiationFactors'
     
     SAVE TimesVisited, FirstTime, nMap, MapFine, MapCoarse, MapW, &
         FineRef, CoarseRef, AnchorUV, AnchorElem

!-------------------------------------------------------------------------------------------
     
     Model  => CurrentModel
     IF (.NOT. ASSOCIATED(Model)) THEN
       CALL Fatal(Caller,'No pointer to model')
     END IF

     Mesh => TSolver % Mesh
     IF (.NOT. ASSOCIATED(Mesh) ) THEN
       CALL Fatal(Caller,'No pointer to mesh')
     END IF
     CALL SetCurrentMesh( Model, Mesh )

     Params => TSolver % Values

     RadiatorsFound = .FALSE.
     DiffuseGrayRadiationFound = .FALSE.
     
     DO i=1,Model % NumberOfBCs
       BC => Model % BCs(i) % Values
       RadiatorsFound = RadiatorsFound .OR. GetLogical(BC,'Radiator BC',Found)
       RadiationFlag = GetString(BC,'Radiation',Found)
       IF (RadiationFlag == 'diffuse gray') DiffuseGrayRadiationFound = .TRUE.
     END DO
     IF(.NOT. DiffuseGrayRadiationFound .AND. .NOT. RadiatorsFound) RETURN

     Radiosity = GetLogical( Params, 'Radiosity Model', Found )
     Spectral = GetLogical( Params, 'Spectral Model' ,Found )
     IF( Spectral ) Radiosity = .TRUE.

     FineMesh => Mesh
     RadLevel = ListGetInteger( Params,'Radiation Relative Mesh Level',Hybrid )
     Hybrid = Hybrid .AND. RadLevel < 0
     RadMeshName = ListGetString( Params,'Radiation Mesh',GeneralMesh )
     IF( GeneralMesh ) THEN
       IF( Hybrid ) CALL Fatal(Caller,&
           'Give either "Radiation Mesh" or "Radiation Relative Mesh Level", not both!')
       Hybrid = .TRUE.
     END IF
     ! Factor files are in the directory of the mesh on disk. In the hybrid case the
     ! parent meshes of the hierarchy are renamed, while the finest keeps the name.
     MeshDirName = Mesh % Name
     RadDirName = Mesh % Name
     VFSuffix = ''
     RadSuffix = ''
     IF( Hybrid ) CALL SetHybridRadiationMesh()

     IF(.NOT. TopoCall) TimesVisited = TimesVisited + 1

     UpdateViewFactors = GetLogical( Params, 'Update View Factors', Found )

     GeometryFixedAfter = GetInteger( Params, &
         'View Factors Fixed After Iterations',Found)
     IF(.NOT. Found) GeometryFixedAfter = HUGE(GeometryFixedAfter)

     IF( UpdateViewFactors ) THEN      
       IF(GeometryFixedAfter < TimesVisited) UpdateViewFactors = .FALSE.
       IF(TimesVisited > 1 ) THEN
         SteadyChange = TSolver % Variable % SteadyChange
         Tol = GetConstReal( Params, 'View Factors Fixed Tolerance',Found)
         IF(Found .AND. SteadyChange < Tol) UpdateViewFactors = .FALSE.
       END IF
     END IF 

     CALL GetGebhartFactorsParameters()
     UpdateRadiatorFactors = GetLogical(Params,'Update Radiator Factors',Found)

     IF(.NOT. (FirstTime .OR. UpdateViewFactors .OR. UpdateGebhartFactors .OR. &
             UpdateRadiatorFactors .OR. Radiosity)) RETURN

!------------------------------------------------------------------------------
!    Go for it
!------------------------------------------------------------------------------
     at = CPUTime(); rt = RealTime()

     CALL Info(Caller,'----------------------------------------------------',Level=5)
     CALL Info(Caller,'Computing radiation factors for heat transfer',       Level=5)
     CALL Info(Caller,'----------------------------------------------------',Level=10)

     FullMatrix = GetLogical( Params, 'Radiation Factors Solver Full',Found) 
     IF(.NOT.Found) &
       FullMatrix = GetLogical( Params, 'Gebhart Factors Solver Full',Found) 
     IF(.NOT.Found) &
       FullMatrix = GetLogical( Params, 'Gebhardt Factors Solver Full',Found) 

     IterSolveFactors = GetLogical( Params, 'Radiation Factors Solver Iterative',Found) 
     IF(.NOT.Found) &
       IterSolveFactors  =  GetLogical( Params, 'Gebhart Factors Solver Iterative',Found) 
     IF(.NOT.Found) &
        IterSolveFactors =  GetLogical( Params, 'Gebhardt Factors Solver Iterative',Found) 
     IF(.NOT. Found) THEN
       SolverType = GetString( Params, 'radiation: Linear System Solver', Found )
       IF( Found ) THEN
         IF( SolverType == 'iterative' ) IterSolveFactors = .TRUE. 
       END IF
     END IF

     IF( FirstTime ) THEN
       IF( FullMatrix ) THEN
         CALL Fatal(Caller, &
             'Using full matrix format for radiation problems not available anymore.')
       ELSE
         CALL Info(Caller,'Using sparse matrix format for factor computations.',Level=6)
       END IF
       
       IF( IterSolveFactors ) THEN
         CALL Info(Caller,'Using iterative solver for radiation factors',Level=6)
       ELSE
         CALL Info(Caller,'Using direct solver for radiation factors',Level=6)
       END IF
     END IF

     ComputeViewFactors = GetLogical( Params, 'Compute View Factors',Found )
     ComputeRadiatorFactors = GetLogical( Params, 'Compute Radiator Factors',Found )

!------------------------------------------------------------------------------
!    Compute the number of elements at the surface and check if the 
!    geometry has really changed.
!------------------------------------------------------------------------------
     RadiationSurfaces = 0
     MaxRadiationBody  = 1

     ALLOCATE(ActiveNodes(Mesh % NumberOfNodes), STAT=istat )
     IF (istat/= 0) CALL Fatal(Caller,'Memory allocation error 1.')
     ActiveNodes = .FALSE.

     nBulk = Mesh % NumberOfBulkElements
     nBndr = Mesh % NumberOfBoundaryElements
     ALLOCATE(ElementNumbers(nBndr), InvElementNumbers(nBndr), &
            RelAreas(nBndr), Areas(nBndr), STAT=istat)
     IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 2.')

     CALL GetMeshRadiationSurfaceInfo()

     ALLOCATE(ActiveMe(0:ParEnv % PEs-1), ActiveTasks(0:ParEnv % PEs-1))
     ActiveMe = .FALSE.
     ActiveMe(ParEnv % myPE) = RadiationSurfaces > 0
     IF ( ParEnv % PEs>1 )THEN
       CALL MPI_ALLREDUCE( ActiveMe, ActiveTasks, ParEnv % PEs, &
          MPI_LOGICAL, MPI_LOR, ELMER_COMM_WORLD, i )
     ELSE
       ActiveTasks = ActiveMe
     END IF

     IF ( RadiationSurfaces == 0 ) THEN
       IF( FirstTime ) THEN
         CALL Info(Caller,'No surfaces participating in radiation',Level=5)
       END IF
       RETURN
     END IF

     ! Check that the geometry has really changed before computing the viewfactors 
     IF(.NOT. FirstTime .AND. (UpdateViewFactors .OR. UpdateRadiatorFactors)) THEN
       IF( .NOT. CheckMeshHasChanged() ) THEN
         UpdateViewFactors = .FALSE.
         UpdateRadiatorFactors = .FALSE.
       END IF         
     END IF

     ! If the geometry has not changed and Gebhart factors are fine return
     IF(.NOT. (FirstTime .OR. UpdateViewFactors .OR. UpdateGebhartFactors .OR. &
         UpdateRadiatorFactors .OR. Radiosity)) THEN
       CALL Info(Caller,'Not first time and no updates needed!',Level=12)
       RETURN
     END IF

     IF( FirstTime .OR. UpdateViewFactors .OR. UpdateRadiatorFactors ) THEN
       ! This stays fixed unless the geometry changes. 
       CALL Info(Caller,'Total number of Radiation Surfaces '//I2S(RadiationSurfaces)// &
           ' out of '//I2S(Model % NumberOfBoundaryElements),Level=5)
     END IF

!-----------------------------------------------------------------------------------
!    Check that the needed files exist if os assumed, if not, recompute
!    view factors and radiator factors
!-----------------------------------------------------------------------------------
     CALL CheckFactorsFilesExist()

!------------------------------------------------------------------------------
!    Rewrite the nodes for view factor computations if they have changed
!    and compute the view factors and/or radiator factors with an external
!    function call.
!------------------------------------------------------------------------------

     UpdateGeometry = ListGetLogical(Params,'Update Factors Geometry',Found )  
     IF(.NOT. Found ) THEN 
       UpdateGeometry = ComputeViewFactors .OR. (ComputeRadiatorFactors.AND.RadiatorsFound) .OR. &
           (.NOT. FirstTime .AND. (UpdateViewFactors .OR. UpdateRadiatorFactors))       
     END IF

     IF(UpdateGeometry) THEN
       IF(GetLogical( Params,'Viewfactor Rigid Mesh Mapping', Found ) .OR. &
           ListGetLogicalAnySolver(Model,'Viewfactor Mapping Solver') ) THEN 
         CALL Info(Caller,'Viewfactor geometry will be changed by its own rigid mesh mapping!',Level=4)
         UpdateGeometry = .FALSE.
       END IF
     END IF

     CALL ComputeViewFactorsAndRadiators()

     IF(RadiatorsFound) THEN
       IF (FirstTime .OR. UpdateRadiatorFactors) THEN
         IF( Hybrid ) THEN
           ! Direct radiator irradiation is resolved on the fine mesh,
           ! only its reflected part is treated on the coarse mesh.
           CALL ReadFineRadiatorFactors()
         ELSE
           CALL ReadRadiatorFactorsFromFile(Mesh,RadiationSurfaces,ElementNumbers,Areas)
         END IF
       END IF
     END IF

     IF( .NOT. DiffuseGrayRadiationFound ) THEN
       CALL Info(Caller,'No diffuse grey radiation found!',Level=12)
       RETURN       
     END IF

!------------------------------------------------------------------------------

     TopologyFixed = GetLogical( Params, 'Matrix Topology Fixed',Found)
     SaveFactors = ListGetLogical( Params, 'Save Gebhart Factors',Found )
     IF(.NOT. Found) &
       SaveFactors = ListGetLogical( Params, 'Save Gebhardt Factors',Found )
   
     TopologyTest = .NOT. TopoCall

!------------------------------------------------------------------------------

     IF (.NOT. ALLOCATED(TSolver % Mesh % VFStore)) THEN
       ALLOCATE(TSolver % Mesh % VFStore(MaxRadiationBody))
     END IF

     DO RadiationBody = 1,MaxRadiationBody
       bt = CPUTime()

       CALL Info(Caller,'Computing area info for set '//I2S(RadiationBody),Level=12)
       CALL GetBodyRadiationSurfaceInfo(RadiationBody)
       IF(RadiationSurfaces == 0)  CYCLE

       IF(FirstTime .OR. UpdateViewFactors) THEN
         IF ( .NOT. ReadViewFactorsFromFile(RadiationBody)) CYCLE
       END IF

       ! and finally, compute the Gebhart factor or radiosities:
       ! -------------------------------------------------------
       ViewFactors => TSolver % Mesh % VFStore(RadiationBody) % VF
       IF(.NOT.CheckForQuickFactors()) THEN
         IF( MaxRadiationBody > 1 ) &
           CALL Info(Caller,'Computing radiation for set '//I2S(RadiationBody),Level=12)
         CALL CalculateRadiation()
       END IF

       IF(MaxRadiationBody > 1) THEN
         bt = CPUTime()-bt
         WRITE (Message,'(A,T35,ES15.4)') 'Radiation body '//I2S(RadiationBody)//' done (s)',bt
         CALL Info(Caller,Message)
       END IF
     END DO ! RadiationBody

!------------------------------------------------------------------------------
     
     IF(.NOT. (TopoCall .OR. TopologyTest .OR. TopologyFixed .OR. Radiosity) ) THEN       
       CALL UpdateMatrixTopologyWithFactors()
     END IF     

     FirstTime = .FALSE.

     IF( Radiosity ) THEN
       WRITE (Message,'(A,T35,ES15.4)') 'Radiosity vector determined (s)',CPUTime()-at
     ELSE
       WRITE (Message,'(A,T35,ES15.4)') 'Gebhart factors determined (s)',CPUTime()-at
     END IF
     CALL Info(Caller,Message,Level=4)
     CALL Info(Caller,'----------------------------------------------------',Level=5)

   CONTAINS

     ! This is just to enable quicker testing. It applies only when emissivity is one everywhere...
     FUNCTION CheckForQuickFactors() RESULT(FoundQuick)
       LOGICAL :: FoundQuick

       IF( ListGetLogical( Params,'Use ViewFactors As Gebhart Factors',FoundQuick ) ) THEN
         CALL Warn(Caller,'Used ViewFactors for RadiationFactors (assumes eps=1)')
         CALL UseViewFactorsAsGebhartFactors()
         IF(SaveFactors) CALL SaveGebhartFactors()       
         RETURN       
       END IF
     END FUNCTION CheckForQuickFactors


     ! Boundary element of the radiation mesh. In the hybrid case this is not
     ! the mesh of the (heat equation) solver.
     FUNCTION RadBoundaryElement(t) RESULT(Element)
       INTEGER :: t
       TYPE(Element_t), POINTER :: Element

       IF( Hybrid ) THEN
         Element => Mesh % Elements(Mesh % NumberOfBulkElements+t)
       ELSE
         Element => GetBoundaryElement(t)
       END IF
     END FUNCTION RadBoundaryElement


     ! Pick the coarser mesh for radiation and create the mapping from the
     ! boundary elements of the finer mesh to those of the coarser one.
     SUBROUTINE SetHybridRadiationMesh()
       TYPE(Mesh_t), POINTER :: pMesh
       INTEGER :: i,j,k,n

       IF(.NOT. Radiosity) CALL Fatal(Caller,&
           '"Radiation Relative Mesh Level" and "Radiation Mesh" require "Radiosity Model"!')

       IF( GeneralMesh ) THEN
         CALL SetGeneralRadiationMesh()
         RETURN
       END IF

       DO j=-1,RadLevel,-1
         IF(.NOT. ASSOCIATED(Mesh % Parent) ) THEN
           CALL Fatal(Caller,'Could not find radiation relative mesh level: '//I2S(RadLevel))
         END IF
         IF(.NOT. ASSOCIATED(Mesh % BoundaryParent) ) THEN
           CALL Fatal(Caller,'No boundary element mapping to parent mesh, use "Mesh Levels"!')
         END IF
         Mesh => Mesh % Parent
       END DO

       n = FineMesh % NumberOfBoundaryElements
       ALLOCATE( FineParent(n) )
       DO k=1,n
         i = k
         pMesh => FineMesh
         DO WHILE( .NOT. ASSOCIATED(pMesh, Mesh) )
           i = pMesh % BoundaryParent(i)
           IF( i == 0 ) EXIT
           pMesh => pMesh % Parent
         END DO
         FineParent(k) = i
       END DO

       ! Each fine element lies within one coarse element
       IF( ALLOCATED(MapFine) ) DEALLOCATE( MapFine, MapCoarse, MapW )
       nMap = COUNT( FineParent > 0 )
       ALLOCATE( MapFine(nMap), MapCoarse(nMap), MapW(nMap) )
       nMap = 0
       DO k=1,n
         IF( FineParent(k) <= 0 ) CYCLE
         nMap = nMap + 1
         MapFine(nMap) = k
         MapCoarse(nMap) = FineParent(k)
         MapW(nMap) = 1.0_dp
       END DO

       ! The nodes of the parent mesh are the first nodes of the split mesh.
       ! Hence the coarse mesh may follow the finer one, if that has moved.
       n = Mesh % NumberOfNodes
       Mesh % Nodes % x(1:n) = FineMesh % Nodes % x(1:n)
       Mesh % Nodes % y(1:n) = FineMesh % Nodes % y(1:n)
       Mesh % Nodes % z(1:n) = FineMesh % Nodes % z(1:n)

       IF( FirstTime ) THEN
         CALL Info(Caller,'Computing radiosity on relative mesh level '//I2S(RadLevel)//&
             ' with '//I2S(Mesh % NumberOfBoundaryElements)//' vs. '//I2S(n)//' boundary elements',Level=5)
       END IF

       ! The factor files are named by the mesh level they are computed on, e.g.
       ! ViewFactorsL1.dat on the coarse and RadiatorFactorsL2.dat on the fine mesh.
       VFSuffix = 'L'//I2S(MeshLevel(Mesh))
       RadSuffix = 'L'//I2S(MeshLevel(FineMesh))

       CALL SetCurrentMesh( Model, Mesh )
     END SUBROUTINE SetHybridRadiationMesh


     ! Radiation on an independent mesh: load it (in full to each partition) and
     ! compute the mapping from the boundary elements of the heat equation mesh.
     !---------------------------------------------------------------------------
     SUBROUTINE SetGeneralRadiationMesh()
       IF( .NOT. ASSOCIATED( GeneralRadMesh ) ) THEN
         CALL Info(Caller,'Loading radiation mesh: '//RadMeshName,Level=5)
         GeneralRadMesh => LoadMesh2( Model, OutputPath, TRIM(OutputPath)//'/'//RadMeshName, &
             .FALSE., 1, 0 )
         IF(.NOT. ASSOCIATED(GeneralRadMesh)) CALL Fatal(Caller,'Could not load radiation mesh: '//RadMeshName)
         GeneralRadMesh % Name = RadMeshName
         GeneralRadMesh % OutputActive = .FALSE.
       END IF
       Mesh => GeneralRadMesh

       ! View factors are in the directory of the radiation mesh, radiator factors
       ! in that of the heat equation mesh.
       MeshDirName = RadMeshName
       RadDirName = FineMesh % Name

       IF( nMap == 0 ) THEN
         CALL SetReferenceCoordinates()
         CALL BuildGeneralMapping()
         CALL BuildNodeAnchors()
       END IF
       CALL FollowHeatMeshMotion()

       CALL SetCurrentMesh( Model, Mesh )
     END SUBROUTINE SetGeneralRadiationMesh


     ! The mapping between the meshes is created in the reference geometry. For the
     ! heat equation mesh these are the original coordinates if stored, otherwise the
     ! current ones. The radiation mesh is as loaded.
     !---------------------------------------------------------------------------
     SUBROUTINE SetReferenceCoordinates()
       INTEGER :: n
       LOGICAL :: HaveOrig, Found

       n = FineMesh % NumberOfNodes
       HaveOrig = ASSOCIATED( FineMesh % NodesOrig )
       IF( HaveOrig ) HaveOrig = .NOT. ASSOCIATED( FineMesh % NodesOrig, FineMesh % Nodes )
       IF( HaveOrig ) HaveOrig = SIZE( FineMesh % NodesOrig % x ) >= n

       ALLOCATE( FineRef(3,n) )
       IF( HaveOrig ) THEN
         CALL Info(Caller,'Using original coordinates of heat equation mesh as reference',Level=6)
         FineRef(1,:) = FineMesh % NodesOrig % x(1:n)
         FineRef(2,:) = FineMesh % NodesOrig % y(1:n)
         FineRef(3,:) = FineMesh % NodesOrig % z(1:n)
       ELSE
         ! The mesh may have been mapped already when the solvers were created
         IF( ListGetLogicalAnySolver( Model,'Viewfactor Mapping Solver' ) .OR. &
             ListGetLogical( Params,'Viewfactor Rigid Mesh Mapping',Found ) ) THEN
           CALL Fatal(Caller,'With a moving geometry and "Radiation Mesh" give '//&
               '"Store Original Coordinates = True" for the mesh mapping solver!')
         END IF
         FineRef(1,:) = FineMesh % Nodes % x(1:n)
         FineRef(2,:) = FineMesh % Nodes % y(1:n)
         FineRef(3,:) = FineMesh % Nodes % z(1:n)
       END IF

       n = Mesh % NumberOfNodes
       ALLOCATE( CoarseRef(3,n) )
       CoarseRef(1,:) = Mesh % Nodes % x(1:n)
       CoarseRef(2,:) = Mesh % Nodes % y(1:n)
       CoarseRef(3,:) = Mesh % Nodes % z(1:n)
     END SUBROUTINE SetReferenceCoordinates


     ! Locate the nodes of the coarse radiation elements on the fine radiation elements
     ! in the reference geometry.
     !---------------------------------------------------------------------------
     SUBROUTINE BuildNodeAnchors()
       TYPE(Element_t), POINTER :: Element, CElement
       TYPE(Nodes_t) :: FNodes
       INTEGER, ALLOCATABLE :: Felem(:), BinPtr(:), BinList(:)
       REAL(KIND=dp), ALLOCATABLE :: FBox(:,:), FSize(:)
       LOGICAL, ALLOCATABLE :: NodeUsed(:)
       REAL(KIND=dp) :: Basis(MAX_ELEMENT_NODES), BoxMin(3), hb, x(3), xc(3), u, v, w, ld, d, best, &
           detJ, RelTol
       INTEGER :: i, j, k, l, n, p, nfe, nbin(3), ib(3), bin, ebest, nn, NoFound
       LOGICAL :: stat, Found

       RelTol = ListGetCReal( Params,'Radiation Mesh Mapping Tolerance',Found )
       IF(.NOT. Found) RelTol = 0.5_dp

       nn = Mesh % NumberOfNodes
       ALLOCATE( AnchorElem(nn), AnchorUV(2,nn), NodeUsed(nn) )
       AnchorElem = 0
       AnchorUV = 0.0_dp
       NodeUsed = .FALSE.
       DO i=1,Mesh % NumberOfBoundaryElements
         CElement => Mesh % Elements(Mesh % NumberOfBulkElements+i)
         IF( IsRadiationBC(CElement) ) NodeUsed(CElement % NodeIndexes) = .TRUE.
       END DO

       ! Local fine radiation elements in the reference geometry
       ALLOCATE( FNodes % x(MAX_ELEMENT_NODES), FNodes % y(MAX_ELEMENT_NODES), FNodes % z(MAX_ELEMENT_NODES) )
       n = FineMesh % NumberOfBoundaryElements
       ALLOCATE( Felem(n), FBox(6,n), FSize(n) )
       nfe = 0
       DO k=1,n
         Element => FineMesh % Elements(FineMesh % NumberOfBulkElements+k)
         IF( .NOT. IsRadiationBC(Element) ) CYCLE
         IF( .NOT. ( ASSOCIATED(Element % BoundaryInfo % Left) .OR. &
             ASSOCIATED(Element % BoundaryInfo % Right) ) ) CYCLE
         nfe = nfe + 1
         Felem(nfe) = k
         CALL RefNodes( Element, FNodes )
         j = Element % TYPE % NumberOfNodes
         FSize(nfe) = MAXVAL( [ MAXVAL(FNodes % x(1:j))-MINVAL(FNodes % x(1:j)), &
             MAXVAL(FNodes % y(1:j))-MINVAL(FNodes % y(1:j)), MAXVAL(FNodes % z(1:j))-MINVAL(FNodes % z(1:j)) ] )
         d = RelTol * FSize(nfe)
         FBox(1:3,nfe) = [ MINVAL(FNodes % x(1:j)), MINVAL(FNodes % y(1:j)), MINVAL(FNodes % z(1:j)) ] - d
         FBox(4:6,nfe) = [ MAXVAL(FNodes % x(1:j)), MAXVAL(FNodes % y(1:j)), MAXVAL(FNodes % z(1:j)) ] + d
       END DO
       IF( nfe == 0 ) RETURN

       CALL BuildBins( nfe, FBox, FSize, BoxMin, hb, nbin, BinPtr, BinList )

       NoFound = 0
       DO p=1,nn
         IF( .NOT. NodeUsed(p) ) CYCLE
         x = CoarseRef(:,p)
         ib = MAX(1, MIN(nbin, 1 + FLOOR( (x-BoxMin) / hb )))
         bin = 1 + (ib(1)-1) + nbin(1)*((ib(2)-1) + nbin(2)*(ib(3)-1))
         ebest = 0
         best = HUGE(best)
         DO l=BinPtr(bin),BinPtr(bin+1)-1
           j = BinList(l)
           IF( ANY( x < FBox(1:3,j) ) .OR. ANY( x > FBox(4:6,j) ) ) CYCLE
           Element => FineMesh % Elements(FineMesh % NumberOfBulkElements+Felem(j))
           CALL RefNodes( Element, FNodes )
           CALL GlobalToLocal( u, v, w, x(1), x(2), x(3), Element, FNodes )
           CALL ClampLocal( Element, u, v, ld )
           stat = ElementInfo( Element, FNodes, u, v, 0.0_dp, detJ, Basis )
           n = Element % TYPE % NumberOfNodes
           xc = [ SUM(Basis(1:n)*FNodes % x(1:n)), SUM(Basis(1:n)*FNodes % y(1:n)), &
               SUM(Basis(1:n)*FNodes % z(1:n)) ]
           d = SQRT( SUM( (x-xc)**2 ) )
           IF( d > RelTol * FSize(j) ) CYCLE
           IF( d < best ) THEN
             best = d
             ebest = Felem(j)
             AnchorUV(:,p) = [u,v]
           END IF
         END DO
         AnchorElem(p) = ebest
         IF( ebest == 0 ) NoFound = NoFound + 1
       END DO
       CALL Info(Caller,'Radiation mesh nodes located on heat equation mesh: '//&
           I2S(COUNT(AnchorElem>0))//' out of '//I2S(COUNT(NodeUsed)),Level=6)
     END SUBROUTINE BuildNodeAnchors


     ! Move the radiation mesh with the displacement of the heat equation mesh
     ! relative to the reference geometry.
     !---------------------------------------------------------------------------
     SUBROUTINE FollowHeatMeshMotion()
       TYPE(Element_t), POINTER :: Element
       TYPE(Nodes_t) :: FNodes
       REAL(KIND=dp), ALLOCATABLE :: Disp(:)
       REAL(KIND=dp) :: Basis(MAX_ELEMENT_NODES), detJ, dmax
       INTEGER :: p, n, nn, i
       INTEGER, POINTER :: Ind(:)
       LOGICAL :: stat

       IF( .NOT. ALLOCATED(AnchorElem) ) RETURN
       nn = Mesh % NumberOfNodes
       ALLOCATE( Disp(4*nn), FNodes % x(MAX_ELEMENT_NODES), FNodes % y(MAX_ELEMENT_NODES), &
           FNodes % z(MAX_ELEMENT_NODES) )
       Disp = 0.0_dp
       DO p=1,nn
         IF( AnchorElem(p) == 0 ) CYCLE
         Element => FineMesh % Elements(FineMesh % NumberOfBulkElements+AnchorElem(p))
         CALL RefNodes( Element, FNodes )
         stat = ElementInfo( Element, FNodes, AnchorUV(1,p), AnchorUV(2,p), 0.0_dp, detJ, Basis )
         n = Element % TYPE % NumberOfNodes
         Ind => Element % NodeIndexes
         Disp(4*p-3) = SUM( Basis(1:n) * ( FineMesh % Nodes % x(Ind) - FineRef(1,Ind) ) )
         Disp(4*p-2) = SUM( Basis(1:n) * ( FineMesh % Nodes % y(Ind) - FineRef(2,Ind) ) )
         Disp(4*p-1) = SUM( Basis(1:n) * ( FineMesh % Nodes % z(Ind) - FineRef(3,Ind) ) )
         Disp(4*p) = 1.0_dp
       END DO
       CALL ParallelSum(Disp)

       dmax = 0.0_dp
       DO p=1,nn
         IF( Disp(4*p) > 0.5_dp ) THEN
           Disp(4*p-3:4*p-1) = Disp(4*p-3:4*p-1) / Disp(4*p)
         ELSE
           Disp(4*p-3:4*p-1) = 0.0_dp
         END IF
         dmax = MAX( dmax, MAXVAL(ABS(Disp(4*p-3:4*p-1))) )
         Mesh % Nodes % x(p) = CoarseRef(1,p) + Disp(4*p-3)
         Mesh % Nodes % y(p) = CoarseRef(2,p) + Disp(4*p-2)
         Mesh % Nodes % z(p) = CoarseRef(3,p) + Disp(4*p-1)
       END DO

       IF( dmax > 0.0_dp .AND. .NOT. MotionReported ) THEN
         WRITE(Message,'(A,ES12.3)') 'Radiation mesh follows the heat equation mesh, max displacement: ',dmax
         CALL Info(Caller,Message,Level=5)
         MotionReported = .TRUE.
       END IF
     END SUBROUTINE FollowHeatMeshMotion


     ! Nodes of a fine element in the reference geometry
     SUBROUTINE RefNodes( Element, FNodes )
       TYPE(Element_t), POINTER :: Element
       TYPE(Nodes_t) :: FNodes
       INTEGER :: n
       n = Element % TYPE % NumberOfNodes
       FNodes % x(1:n) = FineRef(1,Element % NodeIndexes)
       FNodes % y(1:n) = FineRef(2,Element % NodeIndexes)
       FNodes % z(1:n) = FineRef(3,Element % NodeIndexes)
     END SUBROUTINE RefNodes


     ! Bins for finding objects (elements) by their bounding boxes. The bin size is about
     ! twice the mean object size.
     !---------------------------------------------------------------------------
     SUBROUTINE BuildBins( nobj, Box, ObjSize, BoxMin, hb, nbin, BinPtr, BinList )
       INTEGER :: nobj, nbin(3)
       REAL(KIND=dp) :: Box(:,:), ObjSize(:), BoxMin(3), hb
       INTEGER, ALLOCATABLE :: BinPtr(:), BinList(:)
       REAL(KIND=dp) :: BoxMax(3)
       INTEGER :: i, j, l, i1, i2, i3, ib0(3), ib1(3), bin, nbins

       BoxMin = MINVAL(Box(1:3,1:nobj),2)
       BoxMax = MAXVAL(Box(4:6,1:nobj),2)
       hb = 2 * SUM(ObjSize(1:nobj)) / nobj
       DO
         nbin = MAX(1, CEILING((BoxMax-BoxMin)/hb))
         IF( PRODUCT(nbin) <= 8*nobj + 1000 ) EXIT
         hb = 1.5_dp * hb
       END DO
       nbins = PRODUCT(nbin)
       ALLOCATE( BinPtr(nbins+1) )
       BinPtr = 0
       DO l=1,2
         DO j=1,nobj
           ib0 = MAX(1, MIN(nbin, 1 + FLOOR( (Box(1:3,j)-BoxMin) / hb )))
           ib1 = MAX(1, MIN(nbin, 1 + FLOOR( (Box(4:6,j)-BoxMin) / hb )))
           DO i3=ib0(3),ib1(3); DO i2=ib0(2),ib1(2); DO i1=ib0(1),ib1(1)
             bin = 1 + (i1-1) + nbin(1)*((i2-1) + nbin(2)*(i3-1))
             IF( l == 1 ) THEN
               BinPtr(bin+1) = BinPtr(bin+1) + 1
             ELSE
               BinList(BinPtr(bin)) = j
               BinPtr(bin) = BinPtr(bin) + 1
             END IF
           END DO; END DO; END DO
         END DO
         IF( l == 1 ) THEN
           BinPtr(1) = 1
           DO i=1,nbins
             BinPtr(i+1) = BinPtr(i+1) + BinPtr(i)
           END DO
           ALLOCATE( BinList(BinPtr(nbins+1)-1) )
         ELSE
           ! Pointers were advanced to the start of the next bin
           DO i=nbins,1,-1
             BinPtr(i+1) = BinPtr(i)
           END DO
           BinPtr(1) = 1
         END IF
       END DO
     END SUBROUTINE BuildBins


     ! Mapping between independent fine (heat) and coarse (radiation) boundary meshes.
     ! Each fine radiation element is projected to the plane (line in 2D) of nearby
     ! coarse radiation elements having the same boundary condition and an aligned
     ! normal, and clipped with them. The weights are the fractions of the overlaps.
     ! The coarse candidates are found by binning the coarse elements.
     !---------------------------------------------------------------------------------
     SUBROUTINE BuildGeneralMapping()
       TYPE(Element_t), POINTER :: Element, CElement
       TYPE(Nodes_t) :: FNodes, CNodes
       INTEGER, ALLOCATABLE :: Cand(:), BinPtr(:), BinList(:), LocC(:), TmpI(:), Mark(:)
       REAL(KIND=dp), ALLOCATABLE :: CBox(:,:), CSize(:), LocW(:), LocD(:), LocS(:), TmpR(:)
       REAL(KIND=dp) :: x0(3), t1(3), t2(3), nf(3), nc(3), BoxMin(3), FBox(6), &
           hb, d, RelTol, MinDot, sw, tol, Overlap, FArea, Total, Lost, &
           Pf(2,16), Pc(2,4), Poly(2,16)
       INTEGER :: nb, nbc, ncand, i, j, k, l, n, nfc, ncc, npoly, ib0(3), ib1(3), i1, i2, i3, &
           bin, nloc, dim, cap, NoLost, nbin(3)
       LOGICAL :: Found

       CALL Info(Caller,'Creating mapping between heat equation and radiation meshes',Level=6)

       dim = CoordinateSystemDimension()
       RelTol = ListGetCReal( Params,'Radiation Mesh Mapping Tolerance',Found )
       IF(.NOT. Found) RelTol = 0.5_dp
       MinDot = ListGetCReal( Params,'Radiation Mesh Mapping Normal Tolerance',Found )
       IF(.NOT. Found) MinDot = 0.5_dp

       ALLOCATE( FNodes % x(MAX_ELEMENT_NODES), FNodes % y(MAX_ELEMENT_NODES), FNodes % z(MAX_ELEMENT_NODES) )
       ALLOCATE( CNodes % x(MAX_ELEMENT_NODES), CNodes % y(MAX_ELEMENT_NODES), CNodes % z(MAX_ELEMENT_NODES) )

       ! Coarse candidates
       nbc = Mesh % NumberOfBoundaryElements
       ALLOCATE( Cand(nbc), CBox(6,nbc), CSize(nbc) )
       ncand = 0
       DO i=1,nbc
         CElement => Mesh % Elements(Mesh % NumberOfBulkElements+i)
         IF( .NOT. IsRadiationBC(CElement) ) CYCLE
         ncand = ncand + 1
         Cand(ncand) = i
         n = CElement % TYPE % NumberOfNodes
         CSize(ncand) = ElementArea(Mesh,CElement,n)
         IF( dim == 3 ) CSize(ncand) = SQRT(CSize(ncand))
         tol = RelTol * CSize(ncand)
         CBox(1,ncand) = MINVAL(Mesh % Nodes % x(CElement % NodeIndexes)) - tol
         CBox(2,ncand) = MINVAL(Mesh % Nodes % y(CElement % NodeIndexes)) - tol
         CBox(3,ncand) = MINVAL(Mesh % Nodes % z(CElement % NodeIndexes)) - tol
         CBox(4,ncand) = MAXVAL(Mesh % Nodes % x(CElement % NodeIndexes)) + tol
         CBox(5,ncand) = MAXVAL(Mesh % Nodes % y(CElement % NodeIndexes)) + tol
         CBox(6,ncand) = MAXVAL(Mesh % Nodes % z(CElement % NodeIndexes)) + tol
       END DO
       IF( ncand == 0 ) CALL Fatal(Caller,'No radiation elements in radiation mesh!')

       CALL BuildBins( ncand, CBox, CSize, BoxMin, hb, nbin, BinPtr, BinList )

       ! Fine elements
       nb = FineMesh % NumberOfBoundaryElements
       cap = 2*nb + 100
       ALLOCATE( MapFine(cap), MapCoarse(cap), MapW(cap), LocC(ncand), LocW(ncand), &
           LocD(ncand), LocS(ncand), Mark(ncand) )
       Mark = 0
       nMap = 0
       Total = 0.0_dp
       Lost = 0.0_dp
       NoLost = 0

       DO k=1,nb
         Element => FineMesh % Elements(FineMesh % NumberOfBulkElements+k)
         IF( .NOT. IsRadiationBC(Element) ) CYCLE
         ! Radiation elements copied from other partitions have no parents
         IF( .NOT. ( ASSOCIATED(Element % BoundaryInfo % Left) .OR. &
             ASSOCIATED(Element % BoundaryInfo % Right) ) ) CYCLE

         ! The fine element in the reference geometry
         n = Element % TYPE % NumberOfNodes
         nfc = CornerCount(Element)
         CALL RefNodes( Element, FNodes )
         nf = NormalVector( Element, FNodes, Check=.FALSE. )
         FArea = ElementArea(FineMesh,Element,n)
         Total = Total + FArea

         FBox(1) = MINVAL(FNodes % x(1:n)); FBox(4) = MAXVAL(FNodes % x(1:n))
         FBox(2) = MINVAL(FNodes % y(1:n)); FBox(5) = MAXVAL(FNodes % y(1:n))
         FBox(3) = MINVAL(FNodes % z(1:n)); FBox(6) = MAXVAL(FNodes % z(1:n))
         ib0 = MAX(1, MIN(nbin, 1 + FLOOR( (FBox(1:3)-BoxMin) / hb )))
         ib1 = MAX(1, MIN(nbin, 1 + FLOOR( (FBox(4:6)-BoxMin) / hb )))

         nloc = 0
         DO i3=ib0(3),ib1(3); DO i2=ib0(2),ib1(2); DO i1=ib0(1),ib1(1)
           bin = 1 + (i1-1) + nbin(1)*((i2-1) + nbin(2)*(i3-1))
           DO l=BinPtr(bin),BinPtr(bin+1)-1
             j = BinList(l)
             IF( Mark(j) == k ) CYCLE
             Mark(j) = k
             IF( ANY( FBox(4:6) < CBox(1:3,j) ) .OR. ANY( FBox(1:3) > CBox(4:6,j) ) ) CYCLE
             CElement => Mesh % Elements(Mesh % NumberOfBulkElements+Cand(j))
             IF( CElement % BoundaryInfo % Constraint /= Element % BoundaryInfo % Constraint ) CYCLE

             i = CElement % TYPE % NumberOfNodes
             ncc = CornerCount(CElement)
             CNodes % x(1:i) = Mesh % Nodes % x(CElement % NodeIndexes)
             CNodes % y(1:i) = Mesh % Nodes % y(CElement % NodeIndexes)
             CNodes % z(1:i) = Mesh % Nodes % z(CElement % NodeIndexes)
             nc = NormalVector( CElement, CNodes, Check=.FALSE. )
             IF( ABS( SUM(nf*nc) ) < MinDot ) CYCLE

             ! Local frame of the coarse element: origin, tangent(s)
             x0 = [ CNodes % x(1), CNodes % y(1), CNodes % z(1) ]
             t1 = [ CNodes % x(2), CNodes % y(2), CNodes % z(2) ] - x0
             t1 = t1 / SQRT(SUM(t1**2))
             t2 = CrossProduct(nc,t1)

             ! Distance of the fine element centre from the coarse plane
             d = ABS( SUM( nc * ( [ SUM(FNodes % x(1:nfc)), SUM(FNodes % y(1:nfc)), &
                 SUM(FNodes % z(1:nfc)) ] / nfc - x0 ) ) )
             IF( d > RelTol * CSize(j) ) CYCLE

             DO i=1,nfc
               Pf(1,i) = SUM( t1 * ( [FNodes % x(i), FNodes % y(i), FNodes % z(i)] - x0 ) )
               Pf(2,i) = SUM( t2 * ( [FNodes % x(i), FNodes % y(i), FNodes % z(i)] - x0 ) )
             END DO
             DO i=1,ncc
               Pc(1,i) = SUM( t1 * ( [CNodes % x(i), CNodes % y(i), CNodes % z(i)] - x0 ) )
               Pc(2,i) = SUM( t2 * ( [CNodes % x(i), CNodes % y(i), CNodes % z(i)] - x0 ) )
             END DO

             IF( nfc == 2 ) THEN
               ! Segments on a line: overlap of the intervals
               Overlap = MIN(MAXVAL(Pf(1,1:2)),MAXVAL(Pc(1,1:2))) - MAX(MINVAL(Pf(1,1:2)),MINVAL(Pc(1,1:2)))
               Overlap = MAX(0.0_dp, Overlap)
             ELSE
               CALL ClipConvex( Pf, nfc, Pc, ncc, Poly, npoly )
               Overlap = 0.0_dp
               IF( npoly >= 3 ) Overlap = ABS( PolyArea(Poly,npoly) )
             END IF
             IF( Overlap <= 1.0d-12 * CSize(j)**(dim-1) ) CYCLE

             nloc = nloc + 1
             LocC(nloc) = Cand(j)
             LocW(nloc) = Overlap
             LocD(nloc) = d
             LocS(nloc) = CSize(j)
           END DO
         END DO; END DO; END DO

         ! Only the closest surface: e.g. the two sides of a thin plate may have the
         ! same boundary condition and be within the tolerance of each other.
         IF( nloc > 1 ) THEN
           d = MINVAL(LocD(1:nloc))
           i = 0
           DO l=1,nloc
             IF( LocD(l) > 2*d + 1.0d-3 * LocS(l) ) CYCLE
             i = i + 1
             LocC(i) = LocC(l)
             LocW(i) = LocW(l)
           END DO
           nloc = i
         END IF

         IF( nloc == 0 ) THEN
           NoLost = NoLost + 1
           Lost = Lost + FArea
           CYCLE
         END IF

         IF( nMap + nloc > cap ) THEN
           cap = 2*cap + nloc
           ALLOCATE( TmpI(cap) ); TmpI(1:nMap) = MapFine(1:nMap); CALL MOVE_ALLOC(TmpI,MapFine)
           ALLOCATE( TmpI(cap) ); TmpI(1:nMap) = MapCoarse(1:nMap); CALL MOVE_ALLOC(TmpI,MapCoarse)
           ALLOCATE( TmpR(cap) ); TmpR(1:nMap) = MapW(1:nMap); CALL MOVE_ALLOC(TmpR,MapW)
         END IF
         sw = SUM(LocW(1:nloc))
         DO l=1,nloc
           nMap = nMap + 1
           MapFine(nMap) = k
           MapCoarse(nMap) = LocC(l)
           MapW(nMap) = LocW(l) / sw
         END DO
       END DO

       CALL Info(Caller,'Number of mapping entries: '//I2S(nMap),Level=6)
       IF( Total > 0.0_dp ) THEN
         WRITE(Message,'(A,ES12.3)') 'Fraction of fine radiation area not mapped: ',Lost/Total
         CALL Info(Caller,Message,Level=6)
       END IF
       IF( NoLost > 0 ) THEN
         CALL Warn(Caller,'Could not map '//I2S(NoLost)//' fine radiation elements to radiation mesh!')
       END IF
     END SUBROUTINE BuildGeneralMapping


     ! Number of corner nodes of a linear boundary element (line, triangle or quad)
     FUNCTION CornerCount(Element) RESULT(n)
       TYPE(Element_t), POINTER :: Element
       INTEGER :: n
       n = Element % TYPE % ElementCode / 100
       IF( n < 2 .OR. n > 4 ) CALL Fatal(Caller,'Radiation mesh mapping for lines, triangles and quads only!')
     END FUNCTION CornerCount


     ! Clip convex polygon P by convex polygon C (Sutherland-Hodgman) in a plane.
     SUBROUTINE ClipConvex( P, np, C, nc, Q, nq )
       INTEGER :: np, nc, nq
       REAL(KIND=dp) :: P(:,:), C(:,:), Q(:,:)
       REAL(KIND=dp) :: W(2,16), Cc(2,4), a(2), b(2), s(2), e(2)
       INTEGER :: i, j, nw
       LOGICAL :: InE, InS

       ! Clip polygon counterclockwise
       Cc(:,1:nc) = C(:,1:nc)
       IF( PolyArea(Cc,nc) < 0.0_dp ) Cc(:,1:nc) = Cc(:,nc:1:-1)

       nq = np
       Q(:,1:np) = P(:,1:np)
       DO i=1,nc
         a = Cc(:,i)
         b = Cc(:,MOD(i,nc)+1)
         nw = nq
         W(:,1:nw) = Q(:,1:nw)
         nq = 0
         IF( nw == 0 ) EXIT
         s = W(:,nw)
         InS = ClipSide(a,b,s) >= 0.0_dp
         DO j=1,nw
           e = W(:,j)
           InE = ClipSide(a,b,e) >= 0.0_dp
           IF( InE ) THEN
             IF( .NOT. InS ) THEN
               nq = nq + 1; Q(:,nq) = ClipCut(a,b,s,e)
             END IF
             nq = nq + 1; Q(:,nq) = e
           ELSE IF( InS ) THEN
             nq = nq + 1; Q(:,nq) = ClipCut(a,b,s,e)
           END IF
           s = e
           InS = InE
         END DO
       END DO

     END SUBROUTINE ClipConvex


     ! Which side of the line a-b the point x is (positive on the left)
     FUNCTION ClipSide(a,b,x) RESULT(f)
       REAL(KIND=dp) :: a(2), b(2), x(2), f
       f = (b(1)-a(1))*(x(2)-a(2)) - (b(2)-a(2))*(x(1)-a(1))
     END FUNCTION ClipSide


     ! Intersection of segment s-e with the line a-b
     FUNCTION ClipCut(a,b,s,e) RESULT(x)
       REAL(KIND=dp) :: a(2), b(2), s(2), e(2), x(2), fs, fe
       fs = ClipSide(a,b,s)
       fe = ClipSide(a,b,e)
       x = s + (e-s) * fs / (fs-fe)
     END FUNCTION ClipCut


     ! Signed area of a planar polygon
     FUNCTION PolyArea(P,n) RESULT(A)
       REAL(KIND=dp) :: P(:,:), A
       INTEGER :: n, i, j
       A = 0.0_dp
       DO i=1,n
         j = MOD(i,n)+1
         A = A + P(1,i)*P(2,j) - P(1,j)*P(2,i)
       END DO
       A = 0.5_dp * A
     END FUNCTION PolyArea


     ! Is the boundary element participating in radiation
     FUNCTION IsRadiationBC(Element) RESULT(IsRad)
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: IsRad
       TYPE(ValueList_t), POINTER :: BC
       LOGICAL :: Found

       IsRad = .FALSE.
       IF ( GetElementFamily(Element)<=1 ) RETURN
       BC => GetBC(Element)
       IF(.NOT. ASSOCIATED(BC)) RETURN
       IsRad = ( GetString(BC,'Radiation',Found) == 'diffuse gray' ) .OR. &
           GetLogical(BC,'Radiator BC',Found)
     END FUNCTION IsRadiationBC


     ! Move local coordinates to the closest point of the reference element,
     ! returning also the distance moved in local coordinates.
     SUBROUTINE ClampLocal( Element, u, v, ld )
       TYPE(Element_t), POINTER :: Element
       REAL(KIND=dp) :: u, v, ld
       REAL(KIND=dp) :: u0, v0, s

       u0 = u; v0 = v
       SELECT CASE( GetElementFamily(Element) )
       CASE(2)
         u = MAX(-1.0_dp, MIN(1.0_dp, u))
         v = 0.0_dp
       CASE(3)
         u = MAX(0.0_dp, u)
         v = MAX(0.0_dp, v)
         s = u + v
         IF( s > 1.0_dp ) THEN
           u = u / s; v = v / s
         END IF
       CASE(4)
         u = MAX(-1.0_dp, MIN(1.0_dp, u))
         v = MAX(-1.0_dp, MIN(1.0_dp, v))
       END SELECT
       ld = ABS(u-u0) + ABS(v-v0)
     END SUBROUTINE ClampLocal


     ! Level of the mesh in the hierarchy, the mesh on disk being level 1.
     FUNCTION MeshLevel(Mesh0) RESULT(Level)
       TYPE(Mesh_t), POINTER :: Mesh0
       INTEGER :: Level
       TYPE(Mesh_t), POINTER :: pMesh

       Level = 1
       pMesh => Mesh0
       DO WHILE( ASSOCIATED(pMesh % Parent) )
         Level = Level + 1
         pMesh => pMesh % Parent
       END DO
     END FUNCTION MeshLevel


     ! Add suffix to a file name before its extension, "name.dat" -> "nameL1.dat".
     FUNCTION SuffixedName(Name,Suffix) RESULT(NewName)
       CHARACTER(*) :: Name, Suffix
       CHARACTER(:), ALLOCATABLE :: NewName
       INTEGER :: i

       i = INDEX(Name,'.',BACK=.TRUE.)
       IF( i > 0 .AND. INDEX(Name(i:),'/') == 0 ) THEN
         NewName = Name(1:i-1)//Suffix//Name(i:)
       ELSE
         NewName = Name//Suffix
       END IF
     END FUNCTION SuffixedName


     ! Fraction of the boundary element owned by this partition, as in the
     ! assembly of boundary conditions. Used to aggregate fine elements to the
     ! coarse ones such that each fine element is accounted for exactly once.
     FUNCTION OwnFraction(Element) RESULT(w)
       TYPE(Element_t), POINTER :: Element
       REAL(KIND=dp) :: w
       INTEGER :: np, no

       IF( ParEnv % PEs <= 1 ) THEN
         w = 1.0_dp
         RETURN
       END IF
       np = 0; no = 0
       IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
         np = np + 1
         IF( Element % BoundaryInfo % Left % PartIndex == ParEnv % myPE ) no = no + 1
       END IF
       IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
         np = np + 1
         IF( Element % BoundaryInfo % Right % PartIndex == ParEnv % myPE ) no = no + 1
       END IF
       w = 0.0_dp
       IF( np > 0 ) w = 1.0_dp * no / np
     END FUNCTION OwnFraction


     ! Sum of coarse element data over the partitions.
     SUBROUTINE ParallelSum(a)
       REAL(KIND=dp) :: a(:)
       REAL(KIND=dp), ALLOCATABLE :: b(:)
       INTEGER :: ierr

       IF( ParEnv % PEs <= 1 ) RETURN
       ALLOCATE( b(SIZE(a)) )
       CALL MPI_ALLREDUCE( a, b, SIZE(a), MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
       a = b
     END SUBROUTINE ParallelSum


     ! In parallel each partition has all the radiation elements, but those of the
     ! other partitions are copies without parents and without correct field values.
     ! Take the element data from the partitions owning the element instead.
     !---------------------------------------------------------------------------
     SUBROUTINE ParallelSurfaceData( a, b, c )
       REAL(KIND=dp) :: a(:)
       REAL(KIND=dp), OPTIONAL :: b(:), c(:)
       REAL(KIND=dp), ALLOCATABLE :: w(:)
       INTEGER :: i

       IF( ParEnv % PEs <= 1 ) RETURN
       ALLOCATE( w(RadiationSurfaces) )
       DO i=1,RadiationSurfaces
         Element => Mesh % Elements(ElementNumbers(i))
         w(i) = OwnFraction(Element)
       END DO
       a(1:RadiationSurfaces) = w * a(1:RadiationSurfaces)
       CALL ParallelSum( a(1:RadiationSurfaces) )
       IF( PRESENT(b) ) THEN
         b(1:RadiationSurfaces) = w * b(1:RadiationSurfaces)
         CALL ParallelSum( b(1:RadiationSurfaces) )
       END IF
       IF( PRESENT(c) ) THEN
         c(1:RadiationSurfaces) = w * c(1:RadiationSurfaces)
         CALL ParallelSum( c(1:RadiationSurfaces) )
       END IF
       CALL ParallelSum( w )
       IF( ANY( w < 0.5_dp ) ) CALL Fatal(Caller,'Radiation element not owned by any partition!')
       a(1:RadiationSurfaces) = a(1:RadiationSurfaces) / w
       IF( PRESENT(b) ) b(1:RadiationSurfaces) = b(1:RadiationSurfaces) / w
       IF( PRESENT(c) ) c(1:RadiationSurfaces) = c(1:RadiationSurfaces) / w
     END SUBROUTINE ParallelSurfaceData


     ! Sum, minimum or maximum of a scalar over the partitions.
     FUNCTION ParallelScalar(x,oper) RESULT(y)
       REAL(KIND=dp) :: x, y
       CHARACTER(*) :: oper
       INTEGER :: ierr

       y = x
       IF( ParEnv % PEs <= 1 ) RETURN
       SELECT CASE(oper)
       CASE('min')
         CALL MPI_ALLREDUCE( x, y, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr )
       CASE('max')
         CALL MPI_ALLREDUCE( x, y, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, ierr )
       CASE DEFAULT
         CALL MPI_ALLREDUCE( x, y, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
       END SELECT
     END FUNCTION ParallelScalar


     ! Tabulate the coarse surface data by aggregating from the finer mesh:
     ! emitted power and absorptivity are conserved over the coarse element.
     ! Also memorize the fine data needed to give the irradiation back.
     !-----------------------------------------------------------------------
     SUBROUTINE TabulateHybrid()
       TYPE(Element_t), POINTER :: Element
       REAL(KIND=dp), POINTER :: Temperature(:)
       INTEGER, POINTER :: TempPerm(:), Inds(:)
       INTEGER, TARGET :: DGInds(27)
       REAL(KIND=dp), ALLOCATABLE :: SumA(:), SumT4(:), SumET4(:), SumE(:), SumAbs(:)
       INTEGER, ALLOCATABLE :: FineIdx(:)
       REAL(KIND=dp) :: A, T, e, ab, r, T4
       INTEGER :: i,j,k,m,n,nb,p
       LOGICAL :: DG, Found, Dummy

       IF(.NOT. ASSOCIATED(TSolver % Variable)) CALL Fatal(Caller, &
           "Radiosity solution can't be completed without the temperature field.")
       Temperature => TSolver % Variable % Values
       TempPerm => TSolver % Variable % Perm

       DG = ListGetLogical(Params, 'Discontinuous Galerkin',Found ) .OR. &
           ListGetLogical(Params, 'DG Reduced Basis',Found )

       nb = FineMesh % NumberOfBoundaryElements
       ALLOCATE( FineElem(nb), FineEmis(nb), FineAbs(nb), FineT(nb), &
           FineArea(nb), FineGrad(nb), FineW(nb), CoarseA(RadiationSurfaces), CoarseS(RadiationSurfaces), &
           CoarseScale(RadiationSurfaces) )
       FineGrad = 0.0_dp; CoarseS = 0.0_dp
       ALLOCATE( CoarseT(RadiationSurfaces), SumA(RadiationSurfaces), SumT4(RadiationSurfaces), &
           SumET4(RadiationSurfaces), SumE(RadiationSurfaces), SumAbs(RadiationSurfaces) )
       SumA = 0.0_dp; SumT4 = 0.0_dp; SumET4 = 0.0_dp; SumE = 0.0_dp; SumAbs = 0.0_dp

       ! Material parameters may depend on fields that live in the fine mesh.
       CALL SetCurrentMesh( Model, FineMesh )

       ! The fine elements covering some coarse element of this radiation body
       ALLOCATE( FineIdx(nb) )
       FineIdx = 0
       DO p=1,nMap
         IF( InvElementNumbers(MapCoarse(p)) > 0 ) FineIdx(MapFine(p)) = -1
       END DO

       Dummy = .FALSE.
       nFine = 0
       DO k=1,nb
         IF( FineIdx(k) == 0 ) CYCLE
         FineIdx(k) = 0

         Element => FineMesh % Elements(FineMesh % NumberOfBulkElements+k)
         n = GetElementNOFNodes(Element)

         CALL GetElementEmissivity(Element, e, ab, r, 0.5_dp, Dummy)

         IF( DG ) THEN
           CALL DgRadiationIndexes(Element,n,DGInds,.TRUE.)
           Inds => DGInds(1:n)
         ELSE
           Inds => Element % NodeIndexes(1:n)
         END IF
         T = 0.0_dp; m = 0
         DO j=1,n
           IF (TempPerm(Inds(j)) > 0) THEN
             T = T + Temperature(TempPerm(Inds(j)))
             m = m + 1
           END IF
         END DO
         IF( m == 0 ) CYCLE
         T = T / m

         nFine = nFine + 1
         FineIdx(k) = nFine
         FineElem(nFine) = Element % ElementIndex
         FineEmis(nFine) = e
         FineAbs(nFine) = ab
         FineT(nFine) = T
         FineArea(nFine) = ElementArea(FineMesh,Element,n)
         FineW(nFine) = OwnFraction(Element)
       END DO

       ! The pieces of fine elements on coarse elements of this radiation body
       ALLOCATE( PcFine(nMap), PcRad(nMap), PcW(nMap) )
       nPc = 0
       DO p=1,nMap
         k = FineIdx(MapFine(p))
         i = InvElementNumbers(MapCoarse(p))
         IF( k <= 0 .OR. i <= 0 ) CYCLE
         nPc = nPc + 1
         PcFine(nPc) = k
         PcRad(nPc) = i
         PcW(nPc) = MapW(p)
       END DO

       DO p=1,nPc
         k = PcFine(p)
         i = PcRad(p)
         A = PcW(p) * FineW(k) * FineArea(k)
         T4 = FineT(k)**4
         e = FineEmis(k)
         SumA(i) = SumA(i) + A
         SumT4(i) = SumT4(i) + A * T4
         SumET4(i) = SumET4(i) + A * e * T4
         SumE(i) = SumE(i) + A * e
         SumAbs(i) = SumAbs(i) + A * FineAbs(k)
       END DO

       CALL SetCurrentMesh( Model, Mesh )

       CALL ParallelSum(SumA); CALL ParallelSum(SumT4); CALL ParallelSum(SumET4)
       CALL ParallelSum(SumE); CALL ParallelSum(SumAbs)

       DO i=1,RadiationSurfaces
         IF( SumA(i) <= 0.0_dp ) THEN
           CALL Fatal(Caller,'No fine mesh elements found for coarse radiation element: '//I2S(i))
         END IF
         Absorptivity(i) = SumAbs(i) / SumA(i)
         Reflectivity(i) = 1.0_dp - Absorptivity(i)
         IF( SumT4(i) > 0.0_dp ) THEN
           Emissivity(i) = SumET4(i) / SumT4(i)
         ELSE
           Emissivity(i) = SumE(i) / SumA(i)
         END IF
         ! Emitted power is conserved over the coarse element of area Areas(i), and
         ! the irradiation of the coarse element is spread over the fine area SumA(i).
         CoarseT(i) = ( SumT4(i) / Areas(i) )**0.25_dp
         CoarseA(i) = SumA(i)
         CoarseScale(i) = Areas(i) / SumA(i)
       END DO

       IF( InfoActive(20) ) THEN
         PRINT *,'Hybrid radiation: fine elements',nFine,' pieces',nPc,' coarse elements',RadiationSurfaces
         PRINT *,'Hybrid radiation: area ratio',SUM(SumA)/SUM(Areas(1:RadiationSurfaces))
         PRINT *,'Hybrid radiation: coarse scale range',MINVAL(CoarseScale),MAXVAL(CoarseScale)
       END IF
     END SUBROUTINE TabulateHybrid


     ! Give the coarse irradiation G back to the fine elements such that the heat
     ! equation sees q = a*G - e*sigma*T^4 in the form q = Fact(1) - e*sigma*T^4
     ! i.e. Fact(1) = a*G where G includes the direct irradiation of the radiators.
     ! As in the single mesh case, the derivative SOL_d is with respect to a uniform
     ! change of temperature, hence it includes the change of irradiation too.
     !------------------------------------------------------------------------------
     SUBROUTINE UpdateHybridFactors(SOL,SOL_d)
       REAL(KIND=dp) :: SOL(:)
       REAL(KIND=dp), OPTIONAL :: SOL_d(:)

       TYPE(Factors_t), POINTER :: RadiosityFactors
       INTEGER :: i,k,p
       REAL(KIND=dp) :: a
       REAL(KIND=dp), ALLOCATABLE :: Irrad(:), Irrad_d(:), FineG(:), FineG_d(:)
       LOGICAL :: Deriv

       Deriv = Newton .AND. PRESENT(SOL_d)

       ! Irradiation from the surfaces only, the direct radiator part is added per fine element
       ALLOCATE( Irrad(RadiationSurfaces), Irrad_d(RadiationSurfaces) )
       Irrad_d = 0.0_dp
       CALL AbsorbedIrradiation(RadiationSurfaces,SOL,Irrad,AbsG=.FALSE.)
       IF( Deriv ) CALL AbsorbedIrradiation(RadiationSurfaces,SOL_d,Irrad_d,AbsG=.FALSE.)

       ! Irradiation of the fine elements from the coarse ones
       ALLOCATE( FineG(nFine), FineG_d(nFine) )
       FineG = FineGrad(1:nFine)
       FineG_d = 0.0_dp
       DO p=1,nPc
         k = PcFine(p)
         i = PcRad(p)
         FineG(k) = FineG(k) + PcW(p) * CoarseScale(i) * Irrad(i)
         IF( Deriv ) FineG_d(k) = FineG_d(k) + PcW(p) * CoarseScale(i) * Irrad_d(i)
       END DO

       DO k=1,nFine
         RadiosityFactors => FineFactors(k)
         a = FineAbs(k)
         RadiosityFactors % Factors(1) = a * FineG(k)
         IF( Deriv ) RadiosityFactors % Factors(2) = a * FineG_d(k)
       END DO

       IF(InfoActive(30)) THEN
         PRINT *,'Irradiation range:',MINVAL(Irrad),MAXVAL(Irrad),SUM(Irrad)/SIZE(Irrad)
       END IF
     END SUBROUTINE UpdateHybridFactors


     ! The radiosity factors of the k:th fine element of the hybrid scheme.
     !---------------------------------------------------------------------
     FUNCTION FineFactors(k) RESULT ( RadiosityFactors )
       INTEGER :: k
       TYPE(Factors_t), POINTER :: RadiosityFactors
       TYPE(Element_t), POINTER :: Element

       Element => FineMesh % Elements(FineElem(k))
       RadiosityFactors => Element % BoundaryInfo % RadiationFactors
       IF ( .NOT. ASSOCIATED( RadiosityFactors ) ) THEN
         ALLOCATE(RadiosityFactors)
         Element % BoundaryInfo % RadiationFactors => RadiosityFactors
       END IF
       IF (.NOT.ALLOCATED(RadiosityFactors % Elements)) THEN
         ALLOCATE( RadiosityFactors % Elements(1) )
         ALLOCATE( RadiosityFactors % Factors(4) )
         RadiosityFactors % Factors = 0.0_dp
         RadiosityFactors % NumberOfFactors = 1
         RadiosityFactors % Elements(1) = FineElem(k)
       END IF
     END FUNCTION FineFactors


     ! Spectral emissivity and absorptivity of the fine elements for radiation of
     ! temperature Trad. With simple temperature dependence the values depend only on
     ! the keyword list and are hence computed once for each list.
     !-----------------------------------------------------------------------------------
     SUBROUTINE TabulateFineSpectral(Trad,IsRadiator,SimpleTdep,FineKey,Ef,Af)
       REAL(KIND=dp) :: Trad, Ef(:), Af(:)
       LOGICAL :: IsRadiator, SimpleTdep
       INTEGER :: FineKey(:)

       TYPE(Variable_t), POINTER :: TVar
       TYPE(ValueList_t), POINTER :: Vlist
       TYPE(Element_t), POINTER :: Element
       REAL(KIND=dp), ALLOCATABLE :: SaveValues(:), CacheE(:), CacheA(:)
       LOGICAL, ALLOCATABLE :: CacheOk(:)
       INTEGER :: k, key, nkeys

       CALL SetCurrentMesh( Model, FineMesh )

       IF(.NOT. SimpleTdep ) THEN
         TVar => VariableGet(FineMesh % Variables,'Temperature')
         IF(.NOT. ASSOCIATED(TVar)) CALL Fatal(Caller,'Temperature not found in fine mesh!')
         ALLOCATE( SaveValues(SIZE(TVar % Values) ) )
         SaveValues = TVar % Values
         TVar % Values = Trad
       END IF

       nkeys = CurrentModel % NumberOfBCs + CurrentModel % NumberOfMaterials
       ALLOCATE( CacheOk(nkeys), CacheE(nkeys), CacheA(nkeys) )
       CacheOk = .FALSE.

       DO k=1,nFine
         key = FineKey(k)
         IF( SimpleTdep ) THEN
           IF( CacheOk(key) ) THEN
             Ef(k) = CacheE(key); Af(k) = CacheA(key)
             CYCLE
           END IF
         END IF
         Element => FineMesh % Elements(FineElem(k))
         IF( key <= CurrentModel % NumberOfBCs ) THEN
           Vlist => CurrentModel % BCs(key) % Values
         ELSE
           Vlist => CurrentModel % Materials(key-CurrentModel % NumberOfBCs) % Values
         END IF
         CALL GetSpectralEmissivity( Element, Vlist, Trad, IsRadiator, SimpleTdep, Ef(k), Af(k) )
         IF( SimpleTdep ) THEN
           CacheOk(key) = .TRUE.; CacheE(key) = Ef(k); CacheA(key) = Af(k)
         END IF
       END DO

       IF(.NOT. SimpleTdep ) THEN
         TVar % Values = SaveValues
       END IF

       CALL SetCurrentMesh( Model, Mesh )
     END SUBROUTINE TabulateFineSpectral


     ! Spectral radiosity where the view factors are on the coarse mesh while the
     ! temperatures, emissivities and radiators are resolved on the fine mesh.
     ! For each temperature interval (and radiator set) the radiation emitted by the
     ! fine elements is aggregated to the coarse ones, and the irradiation computed
     ! on the coarse mesh is absorbed on the fine elements with their own absorptivity
     ! of that interval. Hence Fact(1) is the absorbed irradiation as in SpectralRadiosity.
     !--------------------------------------------------------------------------------------
     SUBROUTINE SpectralHybrid()
       REAL(KIND=dp) :: Tmin, Tmax, dT, Trad, q, qsum, totsum, c, r, a, T, w, Black
       INTEGER :: i,j,k,p,kmin,kmax,nc
       LOGICAL :: RBC, ApproxNewton, AccurateNewton, SimpleTdep
       INTEGER, ALLOCATABLE :: RadiatorSet(:), FineKey(:)
       TYPE(Element_t), POINTER :: Element
       TYPE(ValueList_t), POINTER :: Vlist
       TYPE(Factors_t), POINTER :: RadiosityFactors
       REAL(KIND=dp), ALLOCATABLE :: RadiatorPowers(:), RadiatorTemps(:), &
           RHS(:), RHS_d(:), SOL(:), SOL_d(:), Diag(:), Ec(:), Ec_d(:), SumAbs(:), &
           Gc(:), Gc_d(:), Ef(:), Af(:), Grad(:), Fact1(:), Fact2(:), WAbs(:), WTemp(:)

       nc = RadiationSurfaces
       ALLOCATE( RHS(nc), SOL(nc), Diag(nc), Ec(nc), SumAbs(nc), Gc(nc) )
       ALLOCATE( Ef(nFine), Af(nFine), Grad(nFine), Fact1(nFine), WAbs(nFine), WTemp(nFine) )
       Fact1 = 0.0_dp; WAbs = 0.0_dp; WTemp = 0.0_dp

       ApproxNewton = .FALSE.
       AccurateNewton = .FALSE.
       IF( Newton ) THEN
         AccurateNewton = ListGetLogical( TSolver % Values,'Accurate Spectral Newton',Found )
         ApproxNewton = .NOT. AccurateNewton
         ALLOCATE( Gc_d(nc), Fact2(nFine) )
         Gc_d = 0.0_dp; Fact2 = 0.0_dp
         IF( AccurateNewton ) ALLOCATE( RHS_d(nc), SOL_d(nc), Ec_d(nc) )
       END IF

       SimpleTdep = ListGetLogical( TSolver % Values,'Radiosity Simple Temperature Dependence',Found)

       ! The keyword list for emissivity of each fine element
       ALLOCATE( FineKey(nFine) )
       CALL SetCurrentMesh( Model, FineMesh )
       DO k=1,nFine
         Element => FineMesh % Elements(FineElem(k))
         Vlist => GetEmissivityList( Element, FineKey(k) )
       END DO
       CALL SetCurrentMesh( Model, Mesh )

       Tmin = HUGE(Tmin); Tmax = -HUGE(Tmax)
       IF( nFine > 0 ) THEN
         Tmin = MINVAL(FineT(1:nFine))
         Tmax = MAXVAL(FineT(1:nFine))
       END IF
       ! All partitions must go through the same intervals
       Tmin = ParallelScalar(Tmin,'min')
       Tmax = ParallelScalar(Tmax,'max')
       IF( Tmin < 0.0_dp ) THEN
         CALL Fatal('SpectralHybrid','Negative temperature not a good starting point!')
       END IF

       dT = ListGetCReal( TSolver % Values,'Spectral dT',UnfoundFatal=.TRUE.)
       kmin = FLOOR( Tmin / dT )
       kmax = CEILING( Tmax / dT )
       CALL Info('SpectralHybrid','Going through discrete intervals: '&
           //I2S(kmin)//'-'//I2S(kmax),Level=6)

       totsum = 0.0_dp
       DO k = kmin, kmax
         qsum = 0.0_dp
         DO j=1,nFine
           q = FineT(j) / dT - k
           IF( ABS(q) < 1 ) qsum = qsum + FineW(j) * (1 - ABS(q))
         END DO
         qsum = ParallelScalar(qsum,'sum')
         IF(qsum < 1.0d-6 ) CYCLE
         totsum = totsum + qsum

         Trad = k*dT
         CALL TabulateFineSpectral(Trad,.FALSE.,SimpleTdep,FineKey,Ef,Af)

         ! Aggregate absorptivity and emitted power of this interval to coarse elements
         SumAbs = 0.0_dp; Ec = 0.0_dp
         IF( AccurateNewton ) Ec_d = 0.0_dp
         DO p=1,nPc
           j = PcFine(p)
           i = PcRad(p)
           SumAbs(i) = SumAbs(i) + PcW(p) * FineW(j) * FineArea(j) * Af(j)
           q = FineT(j) / dT - k
           IF( ABS(q) < 1 ) THEN
             T = FineT(j)
             Black = Sigma * T**4
             w = PcW(p) * FineW(j) * FineArea(j) * (1-ABS(q)) * Ef(j)
             Ec(i) = Ec(i) + w * Black
             IF( AccurateNewton ) Ec_d(i) = Ec_d(i) + w * 4 * Black / T
           END IF
         END DO
         CALL ParallelSum(SumAbs); CALL ParallelSum(Ec)
         IF( AccurateNewton ) CALL ParallelSum(Ec_d)
         Absorptivity(1:nc) = SumAbs / CoarseA
         Ec = Ec / Areas(1:nc)
         IF( AccurateNewton ) Ec_d = Ec_d / Areas(1:nc)

         IF ( UseFullMatrix ) THEN
           G_full = 0.0_dp
         ELSE
           G % Values = 0.0_dp
         END IF
         CALL RadiosityAssembly(nc,G,Diag)
         DO i=1,nc
           a = Absorptivity(i)
           c = RelAreas(i) / a
           RHS(i) = -c * Ec(i)
           IF( AccurateNewton ) RHS_d(i) = -c * Ec_d(i)
         END DO
         CALL BlackRadiosityToRHS(nc,RHS)
         IF( AccurateNewton ) CALL BlackRadiosityToRHS(nc,RHS_d)
         CALL RadiationLinearSolver(nc,G,SOL,RHS,Diag,Solver)
         IF( AccurateNewton ) THEN
           CALL RadiationLinearSolver(nc,G,SOL_d,RHS_d,Diag,Solver,Scaling=.FALSE.)
         END IF

         ! Irradiation of the coarse elements from this interval
         CALL AbsorbedIrradiation(nc,SOL,Gc,AbsG=.FALSE.)
         IF( ApproxNewton ) THEN
           Gc_d = 4 * Gc / Trad
         ELSE IF( AccurateNewton ) THEN
           CALL AbsorbedIrradiation(nc,SOL_d,Gc_d,AbsG=.FALSE.)
         END IF

         ! ... and absorbed by the fine elements
         DO p=1,nPc
           j = PcFine(p)
           i = PcRad(p)
           w = Af(j) * PcW(p) * CoarseScale(i) * Gc(i)
           Fact1(j) = Fact1(j) + w
           WAbs(j) = WAbs(j) + Af(j) * w
           WTemp(j) = WTemp(j) + Trad * w
           IF( Newton ) Fact2(j) = Fact2(j) + Af(j) * PcW(p) * CoarseScale(i) * Gc_d(i)
         END DO
       END DO

       ! This should be exactly one!
       WRITE(Message,'(A,G12.5)') 'Checksum for radiosity sources: ',totsum / ParallelScalar(SUM(FineW(1:nFine)),'sum')
       CALL Info('SpectralHybrid',Message,Level=5)

       ! Radiators: direct irradiation on fine elements, reflected part via coarse mesh.
       RBC = CheckForRadiators(RadiatorPowers,RadiatorTemps)
       IF(RBC) THEN
         ALLOCATE( RadiatorSet(SIZE(RadiatorTemps)) )
         RadiatorSet = 0
         kmax = 0
         DO i=1,SIZE(RadiatorTemps)
           IF(RadiatorSet(i) > 0) CYCLE
           kmax = kmax+1
           RadiatorSet(i) = kmax
           DO j=i+1,SIZE(RadiatorTemps)
             IF(ABS(RadiatorTemps(i)-RadiatorTemps(j)) < 1.0e-6) RadiatorSet(j) = kmax
           END DO
         END DO
         CALL Info('SpectralHybrid','Going through radiators in '//I2S(kmax)//' sets',Level=6)

         DO k = 1, kmax
           DO j=1,SIZE(RadiatorSet)
             IF(RadiatorSet(j) == k) Trad = RadiatorTemps(j)
           END DO
           CALL TabulateFineSpectral(Trad,.TRUE.,SimpleTdep,FineKey,Ef,Af)

           SumAbs = 0.0_dp; Ec = 0.0_dp
           DO j=1,nFine
             Grad(j) = 0.0_dp
             Element => FineMesh % Elements(FineElem(j))
             IF(ALLOCATED(Element % BoundaryInfo % Radiators)) THEN
               DO i=1,SIZE(RadiatorSet)
                 IF(RadiatorSet(i) == k) Grad(j) = Grad(j) + &
                     Element % BoundaryInfo % Radiators(i) * RadiatorPowers(i)
               END DO
             END IF
           END DO
           DO p=1,nPc
             j = PcFine(p)
             i = PcRad(p)
             SumAbs(i) = SumAbs(i) + PcW(p) * FineW(j) * FineArea(j) * Af(j)
             Ec(i) = Ec(i) + PcW(p) * FineW(j) * FineArea(j) * (1-Af(j)) * Grad(j)
           END DO
           CALL ParallelSum(SumAbs); CALL ParallelSum(Ec)
           Absorptivity(1:nc) = SumAbs / CoarseA
           Ec = Ec / Areas(1:nc)

           IF ( UseFullMatrix ) THEN
             G_full = 0.0_dp
           ELSE
             G % Values = 0.0_dp
           END IF
           CALL RadiosityAssembly(nc,G,Diag)
           DO i=1,nc
             a = Absorptivity(i)
             c = RelAreas(i) / a
             RHS(i) = -c * Ec(i)
           END DO
           CALL BlackRadiosityToRHS(nc,RHS)
           CALL RadiationLinearSolver(nc,G,SOL,RHS,Diag,Solver)
           CALL AbsorbedIrradiation(nc,SOL,Gc,AbsG=.FALSE.)

           ! Direct irradiation by the radiators ...
           DO j=1,nFine
             w = Af(j) * Grad(j)
             Fact1(j) = Fact1(j) + w
             WAbs(j) = WAbs(j) + Af(j) * w
             WTemp(j) = WTemp(j) + Trad * w
           END DO
           ! ... and their reflected part via the coarse mesh
           DO p=1,nPc
             j = PcFine(p)
             i = PcRad(p)
             w = Af(j) * PcW(p) * CoarseScale(i) * Gc(i)
             Fact1(j) = Fact1(j) + w
             WAbs(j) = WAbs(j) + Af(j) * w
             WTemp(j) = WTemp(j) + Trad * w
           END DO
         END DO
       END IF

       ! Store the results for access by e.g. heat equation solvers
       DO j=1,nFine
         RadiosityFactors => FineFactors(j)
         RadiosityFactors % Factors(1) = Fact1(j)
         IF( Newton ) RadiosityFactors % Factors(2) = Fact2(j)
         IF( ABS(Fact1(j)) > TINY(w) ) THEN
           RadiosityFactors % Factors(3) = WAbs(j) / Fact1(j)
           RadiosityFactors % Factors(4) = WTemp(j) / Fact1(j)
         END IF
       END DO

       IF(InfoActive(30)) THEN
         PRINT *,'Absorbed range:',MINVAL(Fact1),MAXVAL(Fact1),SUM(Fact1)/nFine
       END IF
     END SUBROUTINE SpectralHybrid


     SUBROUTINE GetGebhartFactorsParameters()
       TYPE(Variable_t), POINTER :: Var
       LOGICAL :: Found
       INTEGER :: k
       REAL(KIND=dp) :: TOL, SteadyChange

       UpdateGebhartFactors = GetLogical( Params, 'Update Gebhart Factors',Found )
       IF(.NOT.Found ) &
           UpdateGebhartFactors = GetLogical( Params, 'Update Gebhardt Factors',Found )
       
       IF( UpdateGebhartFactors ) THEN       
         FactorsFixedAfter = GetInteger( Params, &
           'Gebhart Factors Fixed After Iterations',Found)

         IF(.NOT.Found) &
           FactorsFixedAfter = GetInteger( Params, &
             'Gebhardt Factors Fixed After Iterations',Found)

         IF( Found ) THEN       
           IF(FactorsFixedAfter < TimesVisited) UpdateGebhartFactors = .FALSE.
         END IF
       
         FactorsFixedAfter = GetInteger( Params, &
           'Gebhart Factors Fixed After Nonlinear Iterations',Found)
         IF (.NOT. Found) &
            FactorsFixedAfter = GetInteger( Params, &
               'Gebhardt Factors Fixed After Nonlinear Iterations',Found)

         IF( Found ) THEN                
           Var => VariableGet( Mesh % Variables, 'nonlin iter' )
           IF( ASSOCIATED( Var ) ) THEN
             k = NINT( Var % Values(1) ) 
             IF(FactorsFixedAfter < k ) UpdateGebhartFactors = .FALSE.
           END IF
         END IF
      
         Tol = ListGetConstReal(TSolver % Values, &
             'Gebhart Factors Fixed Tolerance',Found)

         IF (.NOT. Found) &
           Tol = ListGetConstReal(TSolver % Values, &
               'Gebhardt Factors Fixed Tolerance',Found)

         IF( Found ) THEN
           IF(TimesVisited > 1 ) THEN
             SteadyChange = TSolver % Variable % SteadyChange
             IF( SteadyChange < Tol) UpdateGebhartFactors = .FALSE.
           END IF
         END IF
       END IF
     END SUBROUTINE GetGebhartFactorsParameters



     SUBROUTINE FixGeometryAfter(n,Element, BC, BCind)
       INTEGER :: n, BCind
       TYPE(ValueList_t), POINTER :: BC
       TYPE(Element_t), TARGET :: Element

       LOGICAL :: Found
       REAL(KIND=dp) :: x0(1), y0(1), MeshU(n)
          
       x0(1) = 1.0; y0(1) = 1.0
       MeshU(1:n) = GetReal(BC, 'Mesh Update 1',Found, Element)
       IF(.NOT. Found) THEN
         CALL Info(Caller,'Freezing Mesh Update 1 for bc: '//I2S(BCind) )
         CALL ListAddDepReal(BC,'Mesh Update 1', 'Mesh Update 1',1,x0,y0 )
       END IF
       MeshU(1:n) = GetReal(BC,'Mesh Update 2',Found, Element)
       IF(.NOT. Found) THEN
         CALL Info(Caller,'Freezing Mesh Update 2 for bc: '//I2S(BCind) )
         CALL ListAddDepReal( BC,'Mesh Update 2', 'Mesh Update 2',1,x0,y0 )
       END IF
       MeshU(1:n) = GetReal(BC,'Mesh Update 3',Found, Element)
       IF(.NOT. Found) THEN
         CALL Info(Caller,'Freezing Mesh Update 3 for bc: '//I2S(BCind) )
         CALL ListAddDepReal( BC, 'Mesh Update 3', 'Mesh Update 3',1,x0,y0 )
       END IF
     END SUBROUTINE FixGeometryAfter


     FUNCTION CheckMeshHasChanged() RESULT ( HasChanged )

       CHARACTER(:), ALLOCATABLE :: OutputName
       LOGICAL :: HasChanged,Found
       INTEGER :: i,j,k,iostat
       REAL(KIND=dp) :: ds,dx,dy,dz,maxds,refds,maxind,x,y,z
       LOGICAL :: Binary, SinglePrec

       !USE iso_c_binding
       REAL(c_double) :: Coords(3)
       REAL(c_float) :: SCoords(3)

       
       HasChanged = .FALSE.
       
       ! This is a dirty thrick where the input file is tampered
       CALL Info(Caller,'Checking changes in mesh.nodes file!',Level=5)

       OutputName = TRIM(OutputPath) // '/' // TRIM(MeshDirName) // '/mesh.nodes.new'
       Binary = .FALSE.
       SinglePrec = .FALSE.
       
       INQUIRE(FILE=OutputName,EXIST=Found)
       IF(.NOT. Found) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(MeshDirName) // '/mesh.nodes'
       END IF
       OPEN( VFUnit,File=OutputName,STATUS='old',ACTION='read',IOSTAT=iostat)
       IF(iostat /= 0) THEN
         Binary = .TRUE.
         OPEN( VFUnit,File=TRIM(OutputName)//'.bin',FORM='unformatted',ACCESS='stream',&
             STATUS='old',ACTION='read',IOSTAT=iostat)
         IF(iostat /= 0) THEN
           SinglePrec = .TRUE.
           OPEN( VFUnit,File=TRIM(OutputName)//'.sbin',FORM='unformatted',ACCESS='stream',&
               STATUS='old',ACTION='read',IOSTAT=iostat)
         END IF
       END IF
       
       dx = MAXVAL(Mesh % Nodes % x) - MINVAL(Mesh % Nodes % x)
       dy = MAXVAL(Mesh % Nodes % y) - MINVAL(Mesh % Nodes % y)
       dz = MAXVAL(Mesh % Nodes % z) - MINVAL(Mesh % Nodes % z)
       refds = SQRT(dx*dx+dy*dy+dz*dz)

       maxds = 0.0       
       maxind = 0
       Found = .FALSE.

       DO i=1,Mesh % NumberOfNodes
         IF( Binary ) THEN          
           IF(SinglePrec) THEN
             READ(VFUnit,ERR=10,END=10) j,SCoords
             Coords = SCoords
           ELSE
             READ(VFUnit,ERR=10,END=10) j,Coords
           END IF
         ELSE
           READ(VFUnit,*,ERR=10,END=10) j,k,Coords
         END IF

         IF(i == Mesh % NumberOfNodes) Found = .TRUE.
         IF(ActiveNodes(i)) THEN
           dx = Mesh % Nodes % x(i) - Coords(1)
           dy = Mesh % Nodes % y(i) - Coords(2)
           dz = Mesh % Nodes % z(i) - Coords(3)
           ds = SQRT(dx*dx+dy*dy+dz*dz)
           IF(ds > maxds) THEN
             maxds = ds
             maxind = i
           END IF
         END IF
       END DO

10     CONTINUE

       CLOSE(VFUnit)

       IF(.NOT. Found) THEN
         HasChanged = .TRUE.
         CALL Info(Caller,'Mismatch in coordinates compared to file: '//OutputName)
       ELSE
         WRITE(Message,'(A,E15.5)') 'Maximum geometry alteration on radiation BCs:',maxds
         CALL Info(Caller,Message)

         x = ListGetConstReal(TSolver % Values,'View Factors Geometry Tolerance',Found)
         IF(.NOT. Found) x = 1.0d-8

         IF(maxds <= refds * x) THEN
           CALL Info(Caller,'Geometry change is neglected and old view factors are used')
         ELSE
           HasChanged = .TRUE.
           CALL Info(Caller,'Geometry change requires recomputation of view factors')
         END IF
       END IF
     END FUNCTION CheckMeshHasChanged


     SUBROUTINE GetMeshRadiationSurfaceInfo
       TYPE(ValueList_t), POINTER :: BC
       INTEGER :: i,j,k,n,t,gmax
       LOGICAL :: Found, Parallel
       TYPE(Element_t), POINTER :: Element
       CHARACTER(:), ALLOCATABLE :: RadiationFlag
       INTEGER, ALLOCATABLE :: gPerm(:), iPerm(:)

       Parallel = ParEnv % PEs > 1 .AND. .NOT. GeneralMesh
       IF (Parallel) THEN
         ALLOCATE(gPerm(GetNOFBoundaryElements()), iPerm(GetNOFBoundaryElements()))
         gPerm = 0; iPerm = 0
         DO i=1,GetNOFBoundaryElements()
           Element => RadBoundaryElement(i)
           IF ( GetElementFamily(Element)<=1 ) CYCLE
           gPerm(i) = Element % GElementIndex
           iPerm(i) = i
         END DO
         CALL Sorti(GetNOFBoundaryElements(),gPerm,iPerm)
       END IF


       ElementNUmbers = 0
       DO i=1,GetNOFBoundaryElements()
         j = i
         IF(Parallel) j = iPerm(i)
         Element => RadBoundaryElement(j)
         IF ( GetElementFamily(Element)<=1 ) CYCLE
         BC => GetBC(Element)
         t = GetBCId(Element)
         IF(t<=0) CYCLE
         RadiationFlag = GetString( BC, 'Radiation', Found )
         IF (RadiationFlag=='diffuse gray' .OR. GetLogical(BC,'Radiator BC',Found)) THEN
           RadiationSurfaces = RadiationSurfaces + 1
           n = GetElementNOFNodes(Element)
           ElementNumbers(RadiationSurfaces) = j + nBulk
           Areas(RadiationSurfaces) = ElementArea(Mesh,Element,n)

           k=MAX(1,GetInteger(BC,'Radiation Boundary',Found) )
           MaxRadiationBody = MAX(k, MaxRadiationBody)

           ActiveNodes(Element % NodeIndexes) = .TRUE.
           IF(GeometryFixedAfter == TimesVisited) THEN
             CALL FixGeometryAfter(n,Element,BC,t)
           END IF
         END IF
       END DO
     END SUBROUTINE GetMeshRadiationSurfaceInfo


     SUBROUTINE GetBodyRadiationSurfaceInfo(RadiationBody)
       INTEGER :: RadiationBody

       TYPE(ValueList_t), POINTER :: BC
       INTEGER :: i,j,k,t,n, gmax
       LOGICAL :: Found, Parallel
       INTEGER, ALLOCATABLE :: gPerm(:), iPerm(:)
       TYPE(Element_t), POINTER :: Element
       CHARACTER(:), ALLOCATABLE :: RadiationFlag

       Parallel = ParEnv % PEs > 1 .AND. .NOT. GeneralMesh
       IF (Parallel) THEN
         ALLOCATE(gPerm(GetNOFBoundaryElements()), iPerm(GetNOFBoundaryElements()))
         gPerm = 0; iPerm = 0
         DO i=1,GetNOFBoundaryElements()
           Element => RadBoundaryElement(i)
           IF ( GetElementFamily(Element)<=1 ) CYCLE
           gPerm(i) = Element % GElementIndex
           iPerm(i) = i
         END DO
         CALL Sorti(GetNOFBoundaryElements(),gPerm,iPerm)
       END IF

       ElementNumbers    = 0
       RadiationSurfaces = 0
       DO i=1,GetNOFBoundaryElements()
         j = i
         IF(Parallel) j = iPerm(i)
         Element => RadBoundaryElement(j)
         IF ( GetElementFamily(Element)<=1 ) CYCLE
         BC => GetBC(Element)
         IF(.NOT.ASSOCIATED(BC)) CYCLE

         RadiationFlag = GetString( BC, 'Radiation',Found )
         IF (RadiationFlag == 'diffuse gray' .OR. GetLogical(BC,'Radiator BC',Found)) THEN
           k = MAX(1,GetInteger(BC,'Radiation Boundary',Found) )

           IF(k == RadiationBody) THEN
             RadiationSurfaces = RadiationSurfaces + 1
             ElementNumbers(RadiationSurfaces) = j + nBulk
             n = GetElementNOFNodes(Element)
             Areas(RadiationSurfaces) = ElementArea(Mesh,Element,n)
           END IF
         END IF
       END DO
       n = RadiationSurfaces
       RelAreas(1:n) = Areas(1:n) / MAXVAL(Areas(1:n))

       IF(MaxRadiationBody > 1) THEN
         CALL Info(Caller,'Number of Radiation Surfaces '//I2S(RadiationSurfaces)// &
             ' for boundary '//I2S(RadiationBody),Level=5)
       END IF

       ! Make the inverse of the list of element numbers of boundaries
       InvElementNumbers = 0
       DO i=1,RadiationSurfaces
         IF ( ElementNumbers(i) <= 0) CYCLE
         InvElementNumbers(ElementNumbers(i)-nBulk) = i
       END DO
     END SUBROUTINE GetBodyRadiationSurfaceInfo


     SUBROUTINE ComputeViewFactorsAndRadiators()

       CHARACTER(:), ALLOCATABLE :: cmd, OutputName, OutputName2
       INTEGER :: cmdStatus
       LOGICAL :: DoScale
       INTEGER :: i,j
       REAL(KIND=dp), POINTER :: Wrk(:,:)
       REAL(KIND=dp) :: BackScale(3), Coord(3)

       ! This is a dirty thrick where the input mesh is scaled after loading.
       ! We need to perform scaling and backscaling then here too. 
       IF( UpdateGeometry ) THEN
         CALL Info(Caller,'Temporarily updating the mesh.nodes file!',Level=5)
         
         OutputName  = TRIM(OutputPath) // '/' // TRIM(MeshDirName) // '/mesh.nodes'         
         OutputName2 = TRIM(OutputPath) // '/' // TRIM(MeshDirName) // '/mesh.nodes.orig'         
         CALL RenameF(OutputName, OutputName2)

         DoScale = ListCheckPresent( Model % Simulation,'Coordinate Scaling')
       
         IF( DoScale ) THEN
           Wrk => ListGetConstRealArray( Model % Simulation,'Coordinate Scaling',Found )    
           BackScale = 1.0_dp
           DO i=1,Mesh % MeshDim 
             j = MIN( i, SIZE(Wrk,1) )
             BackScale(i) = 1.0_dp / Wrk(j,1)
           END DO
         END IF

         OPEN( VFUnit,FILE=OutputName, STATUS='unknown' )                 
         DO i=1,Mesh % NumberOfNodes
           Coord(1) = Mesh % Nodes % x(i)
           Coord(2) = Mesh % Nodes % y(i)
           Coord(3) = Mesh % Nodes % z(i)
           IF( DoScale ) Coord = BackScale * Coord
           WRITE( VFUnit,'(i7,i3,3f20.12)' ) i,-1, Coord
         END DO
         CLOSE(VFUnit)
       END IF

       ! Compute the factors using an external program call.
       ! Only rank 0 spawns the subprocess; other ranks wait at a barrier.
       ! Set ELMER_NO_MPI in the process environment so the child ViewFactors
       ! does not try to join the parent's MPI job via inherited PMI/PMIx vars.
       IF (ComputeViewFactors .OR.  .NOT.FirstTime .AND. UpdateViewFactors ) THEN
         IF ( ParEnv % MyPE == 0 ) THEN
           cmd = SpawnCommand('ViewFactors')//' '//TRIM(GetSifName())
           CALL Info('ComputeViewFactorsAndRadiators','Using system call: '//TRIM(cmd),Level=15)
           CALL ElmerSetNoMPI( 1 )
           CALL SystemCommand( cmd, cmdStatus )
           CALL ElmerSetNoMPI( 0 )
           IF ( cmdStatus /= 0 ) CALL Fatal('ComputeViewFactorsAndRadiators', &
               'View factor computation failed, command was: '//TRIM(cmd))
         END IF
         IF ( ParEnv % PEs > 1 ) CALL MPI_Barrier( ELMER_COMM_WORLD, i )
       END IF

       IF( RadiatorsFound ) THEN
         IF (ComputeRadiatorFactors .OR. .NOT.FirstTime .AND. UpdateRadiatorFactors ) THEN
           IF ( ParEnv % MyPE == 0 ) THEN
             cmd = SpawnCommand('Radiators')//' '//TRIM(GetSifName())
             CALL Info('ComputeViewFactorsAndRadiators','Using system call: '//TRIM(cmd),Level=15)
             CALL ElmerSetNoMPI( 1 )
             CALL SystemCommand( cmd, cmdStatus )
             CALL ElmerSetNoMPI( 0 )
             IF ( cmdStatus /= 0 ) CALL Fatal('ComputeViewFactorsAndRadiators', &
                 'Radiator factor computation failed, command was: '//TRIM(cmd))
           END IF
           IF ( ParEnv % PEs > 1 ) CALL MPI_Barrier( ELMER_COMM_WORLD, i )
         END IF
       END IF
     
       ! Set back the original node coordinates to prevent unwanted user errors
       IF( UpdateGeometry ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(MeshDirName) // '/mesh.nodes'         
         OutputName2 = TRIM(OutputPath) // '/' // TRIM(MeshDirName) // '/mesh.nodes.new'         
         CALL RenameF(OutputName, OutputName2)

         OutputName = TRIM(OutputPath) // '/' // TRIM(MeshDirName) // '/mesh.nodes.orig'         
         OutputName2 = TRIM(OutputPath) // '/' // TRIM(MeshDirName) // '/mesh.nodes'         
         CALL RenameF(OutputName, OutputName2)
       END IF
     END SUBROUTINE ComputeViewfactorsAndRadiators


     SUBROUTINE CheckFactorsFilesExist()
       LOGICAL :: FilesExist, Found
       INTEGER :: RadiationBody
       CHARACTER(:), ALLOCATABLE :: ViewFactorsFile, RadiatorFactorsFile, &
            OutputName

       IF( DiffuseGrayRadiationFound ) THEN 
         FilesExist = .TRUE.
         DO RadiationBody = 1, MaxRadiationBody 
           ViewFactorsFile = GetString(Model % Simulation,'View Factors',Found)
           IF ( .NOT.Found ) ViewFactorsFile = 'ViewFactors.dat'
           ViewFactorsFile = SuffixedName(ViewFactorsFile,VFSuffix)
       
           IF ( LEN_TRIM(MeshDirName) > 0 ) THEN
             OutputName = TRIM(OutputPath) // '/' // TRIM(MeshDirName) &
                     // '/' // ViewFactorsFile
           ELSE
             OutputName = ViewFactorsFile
           END IF
           IF(RadiationBody > 1) OutputName = OutputName//I2S(RadiationBody)
           INQUIRE(FILE=OutputName,EXIST=Found)
           IF(.NOT. Found) FilesExist = .FALSE.
         END DO
         IF(.NOT. FilesExist) ComputeViewFactors = .TRUE.
       END IF

       IF( RadiatorsFound ) THEN
         FilesExist = .TRUE.
         DO RadiationBody = 1, MaxRadiationBody 
           RadiatorFactorsFile = GetString(Model % Simulation,'Radiator Factors',Found)
           IF ( .NOT.Found ) RadiatorFactorsFile = 'RadiatorFactors.dat'
           RadiatorFactorsFile = SuffixedName(RadiatorFactorsFile,RadSuffix)
         
           IF ( LEN_TRIM(RadDirName) > 0 ) THEN
             OutputName = TRIM(OutputPath) // '/' // TRIM(RadDirName) &
                     // '/' // RadiatorFactorsFile
           ELSE
             OutputName = RadiatorFactorsFile
           END IF
           IF(RadiationBody > 1) OutputName = OutputName//I2S(RadiationBody)
           INQUIRE(FILE=OutputName,EXIST=Found)
           IF(.NOT. Found) FilesExist = .FALSE.
         END DO
         IF(.NOT. FilesExist) ComputeRadiatorFactors = .TRUE.
       END IF
     END SUBROUTINE CheckFactorsFilesExist


     ! This is an add-on for including point like radiators into the system.
     ! They act as heat sources with known total power that is distributed among the
     ! surface elements that the point sees.
     !----------------------------------------------------------------------------------     
     SUBROUTINE ReadRadiatorFactorsFromFile(RMesh,nSurf,ElemNums,ElemAreas)
       TYPE(Mesh_t), POINTER :: RMesh
       INTEGER :: nSurf, ElemNums(:)
       REAL(KIND=dp) :: ElemAreas(:)

       LOGICAL :: Success 
       
       TYPE(BoundaryInfo_t), POINTER :: BoundaryInfo
       REAL(KIND=dp), ALLOCATABLE :: Vals(:)
       INTEGER, ALLOCATABLE ::  Cols(:)
       INTEGER :: NofRadiators, i,j,t,n
       LOGICAL :: BinaryMode, SinglePrec, Found
       REAL(KIND=dp), POINTER :: Radiators(:,:)
       REAL :: sval
       TYPE(ValueList_t), POINTER :: RadList
       CHARACTER(:), ALLOCATABLE :: RadiatorFactorsFile, OutputName

       CALL Info(Caller,'Loading radiator factors!',Level=7)
       Success = .TRUE.
       
       RadiatorFactorsFile = GetString(Model % Simulation,'Radiator Factors',Found)       
       IF ( .NOT.Found ) RadiatorFactorsFile = 'RadiatorFactors.dat'
       RadiatorFactorsFile = SuffixedName(RadiatorFactorsFile,RadSuffix)
       
       IF ( LEN_TRIM(RadDirName) > 0 ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(RadDirName) // &
             '/' // RadiatorFactorsFile
       ELSE
         OutputName = RadiatorFactorsFile
       END IF

       INQUIRE(FILE=OutputName,EXIST=Found)
       IF(.NOT. Found) THEN
         CALL Warn(Caller,'Radiator Factors File does NOT exist: '//TRIM(OutputName))
         Success = .FALSE.
         RETURN
       END IF

       BinaryMode = ListGetLogical( Params,'Radiatorfactor Binary Output',Found ) 
       IF(.NOT. Found) BinaryMode = ListGetLogical( Params,'Viewfactor Binary Output',Found ) 

       IF(BinaryMode) THEN
         SinglePrec = ListGetLogical( Params,'Viewfactor Single Precision',Found ) 
       ELSE
         SinglePrec = .FALSE.
       END IF
         
       IF( BinaryMode ) THEN
         CALL Info(Caller,'Loading radiator factors from binary file: '//OutputName,Level=5)
         OPEN( UNIT=VFUnit, FILE=OutputName, FORM = 'unformatted', &
             ACCESS = 'stream', STATUS='old', ACTION='read' )         
         READ( VFUnit ) n
         IF( n /= nSurf ) THEN
           CALL Fatal(Caller,'Mismatch in radiation factor file size: '&
               //I2S(n)//' vs. '//I2S(nSurf))
         END IF
       ELSE
         CALL Info(Caller,'Loading radiator factors from ascii file: '//OutputName,Level=5)
         OPEN( VFUnit,File=OutputName )
       END IF

       IF( .NOT. ListCheckPresentAnyBodyForce( Model,'Radiator Coordinates',RadList ) ) &
           RadList => Params
       
       CALL GetConstRealArray( RadList, Radiators, 'Radiator Coordinates', Found )
       IF(.NOT. Found ) CALL Fatal( Caller, 'No radiators present, quitting' )

       NofRadiators = SIZE(Radiators,1)

       ! Read in the RadiatorFactors
       DO i=1,NofRadiators
         IF( BinaryMode ) THEN
           READ( VFUnit ) n
         ELSE
           READ( VFUnit,* ) n
         END IF

         ALLOCATE( Vals(n), Cols(n) )
         Vals = 0; Cols = 0

         DO j=1,n
           IF( SinglePrec ) THEN
             READ(VFUnit) Cols(j),sval
             Vals(j) = sval
           ELSE IF( BinaryMode ) THEN
             READ(VFUnit) Cols(j),Vals(j)         
           ELSE
             READ(VFUnit,*) t,Cols(j),Vals(j)         
           END IF
           Vals(j) = Vals(j) / ElemAreas(Cols(j))
           Cols(j) = ElemNums(Cols(j))
         END DO

         DO j=1,n
           IF( Cols(j) <= 0 ) CYCLE
           BoundaryInfo => RMesh % Elements(Cols(j)) % BoundaryInfo
           IF ( .NOT.ALLOCATED( BoundaryInfo % Radiators ) ) THEN
             ALLOCATE( BoundaryInfo % Radiators(NofRadiators) )
             BoundaryInfo % Radiators = 0
           END IF
           BoundaryInfo % Radiators(i) = Vals(j)
         END DO
         DEALLOCATE( Cols, Vals )
       END DO
       CLOSE(VFUnit)
              
     END SUBROUTINE ReadRadiatorFactorsFromFile


     ! Radiator factors for the fine mesh boundary elements, in the element order
     ! of the Radiators program.
     !--------------------------------------------------------------------------
     SUBROUTINE ReadFineRadiatorFactors()
       TYPE(Element_t), POINTER :: Element
       TYPE(ValueList_t), POINTER :: BC
       TYPE(Mesh_t), POINTER :: pMesh
       INTEGER, ALLOCATABLE :: FineNums(:), cKey(:), cPerm(:), cOffset(:), First(:)
       REAL(KIND=dp), ALLOCATABLE :: FineAreas(:)
       INTEGER :: i,j,k,b,n,nb,nf,nrc,nl,pos,mult
       LOGICAL :: Found

       nb = FineMesh % NumberOfBoundaryElements

       IF( ParEnv % PEs > 1 .AND. GeneralMesh ) THEN
         ! All radiator elements of the heat equation mesh are present in each
         ! partition, the serial order is that of their global indexes.
         ALLOCATE( FineNums(nb), FineAreas(nb), cKey(nb), cPerm(nb) )
         nf = 0
         DO j=1,nb
           Element => FineMesh % Elements(FineMesh % NumberOfBulkElements+j)
           IF ( .NOT. IsRadiatorBC(Element) ) CYCLE
           nf = nf + 1
           cKey(nf) = Element % GElementIndex
           cPerm(nf) = j
         END DO
         CALL Sorti(nf,cKey,cPerm)
         DO i=1,nf
           Element => FineMesh % Elements(FineMesh % NumberOfBulkElements+cPerm(i))
           n = GetElementNOFNodes(Element)
           FineNums(i) = Element % ElementIndex
           FineAreas(i) = ElementArea(FineMesh,Element,n)
         END DO
       ELSE IF( ParEnv % PEs <= 1 ) THEN
         ALLOCATE( FineNums(nb), FineAreas(nb) )
         nf = 0
         DO j=1,nb
           Element => FineMesh % Elements(FineMesh % NumberOfBulkElements+j)
           IF ( .NOT. IsRadiatorBC(Element) ) CYCLE
           nf = nf + 1
           n = GetElementNOFNodes(Element)
           FineNums(nf) = Element % ElementIndex
           FineAreas(nf) = ElementArea(FineMesh,Element,n)
         END DO
       ELSE
         ! The fine elements are numbered as in the serial split of the whole mesh,
         ! done by the Radiators program: coarse elements in the order of their global
         ! index, followed by their children. All coarse radiation elements are present
         ! in every partition (RadiationParallelMeshDistribute), hence the numbering
         ! can be computed locally.
         n = Mesh % NumberOfBoundaryElements
         ALLOCATE( cKey(n), cPerm(n), cOffset(n) )
         cOffset = -1
         nrc = 0
         DO i=1,n
           Element => Mesh % Elements(Mesh % NumberOfBulkElements+i)
           IF ( .NOT. IsRadiatorBC(Element) ) CYCLE
           nrc = nrc + 1
           cKey(nrc) = Element % GElementIndex
           cPerm(nrc) = i
         END DO
         CALL Sorti(nrc,cKey,cPerm)

         nl = -RadLevel
         nf = 0
         DO j=1,nrc
           Element => Mesh % Elements(Mesh % NumberOfBulkElements+cPerm(j))
           cOffset(cPerm(j)) = nf
           nf = nf + NumberOfChildren(Element)**nl
         END DO

         ALLOCATE( FineNums(nf), FineAreas(nf) )
         FineNums = 0
         FineAreas = 1.0_dp

         DO k=1,nb
           Element => FineMesh % Elements(FineMesh % NumberOfBulkElements+k)
           IF ( .NOT. IsRadiatorBC(Element) ) CYCLE
           IF ( FineParent(k) <= 0 ) CYCLE
           IF ( cOffset(FineParent(k)) < 0 ) CYCLE
           b = NumberOfChildren(Element)

           ! Position among the descendants of the coarse element, the last split
           ! being the least significant.
           pos = 0
           mult = 1
           i = k
           pMesh => FineMesh
           DO WHILE( .NOT. ASSOCIATED(pMesh, Mesh) )
             j = pMesh % BoundaryParent(i)
             First = FirstChildren(pMesh)
             pos = pos + (i - First(j)) * mult
             mult = mult * b
             i = j
             pMesh => pMesh % Parent
           END DO

           j = cOffset(FineParent(k)) + pos + 1
           n = GetElementNOFNodes(Element)
           FineNums(j) = Element % ElementIndex
           FineAreas(j) = ElementArea(FineMesh,Element,n)
         END DO
       END IF

       CALL Info(Caller,'Reading radiator factors for '//I2S(nf)//' fine mesh elements',Level=7)
       CALL ReadRadiatorFactorsFromFile(FineMesh,nf,FineNums,FineAreas)

     END SUBROUTINE ReadFineRadiatorFactors


     ! Same selection as in ExtractSurfaces of the Radiators program
     FUNCTION IsRadiatorBC(Element) RESULT(IsRad)
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: IsRad
       TYPE(ValueList_t), POINTER :: BC

       IsRad = .FALSE.
       IF ( GetElementFamily(Element)<=1 ) RETURN
       BC => GetBC(Element)
       IF(.NOT. ASSOCIATED(BC)) RETURN
       IsRad = GetLogical(BC,'Radiator BC',Found)
     END FUNCTION IsRadiatorBC

     ! Number of children of a boundary element in SplitMeshEqual
     FUNCTION NumberOfChildren(Element) RESULT(nc)
       TYPE(Element_t), POINTER :: Element
       INTEGER :: nc
       SELECT CASE( GetElementFamily(Element) )
       CASE(2)
         nc = 2
       CASE(3,4)
         nc = 4
       CASE DEFAULT
         nc = 1
       END SELECT
     END FUNCTION NumberOfChildren

     ! The first child of each parent boundary element
     FUNCTION FirstChildren(pMesh) RESULT(First)
       TYPE(Mesh_t), POINTER :: pMesh
       INTEGER, ALLOCATABLE :: First(:)
       INTEGER :: i,j
       ALLOCATE( First(pMesh % Parent % NumberOfBoundaryElements) )
       First = 0
       DO i=pMesh % NumberOfBoundaryElements,1,-1
         j = pMesh % BoundaryParent(i)
         IF( j > 0 ) First(j) = i
       END DO
     END FUNCTION FirstChildren


     FUNCTION ReadViewFactorsFromFile(RadiationBody) RESULT (Success)
       LOGICAL :: Success
       INTEGER :: RadiationBody

       LOGICAL :: Found, BinaryMode, SinglePrec
       INTEGER :: i,j,t,n,n2
       INTEGER, POINTER :: Cols(:)
       REAL(KIND=dp), POINTER :: Vals(:)
       REAL :: sval
       CHARACTER(:), ALLOCATABLE :: ViewFactorsFile, OutputName

       Success = .TRUE.
       
       ViewFactors => TSolver % Mesh % VFStore(RadiationBody) % VF
       IF ( .NOT.ASSOCIATED(ViewFactors) ) THEN
         ALLOCATE( ViewFactors(RadiationSurfaces), STAT=istat )
         IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 4.')
       ELSE
         IF (RadiationSurfaces /= SIZE(ViewFactors)) THEN
           DO i=1,SIZE(ViewFactors)
             IF (ALLOCATED(ViewFactors(i) % Factors)) DEALLOCATE(ViewFactors(i) % Factors)
             IF (ALLOCATED(ViewFactors(i) % Elements)) DEALLOCATE(ViewFactors(i) % Elements)
           END DO
           DEALLOCATE( ViewFactors )
           ALLOCATE( ViewFactors(RadiationSurfaces), STAT=istat )
           IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 5.')
         END IF
       END IF
       TSolver % Mesh % VFStore(RadiationBody) % VF => ViewFactors

       ViewFactorsFile = GetString(Model % Simulation,'View Factors',Found)
       IF ( .NOT.Found ) ViewFactorsFile = 'ViewFactors.dat'
       ViewFactorsFile = SuffixedName(ViewFactorsFile,VFSuffix)

       IF ( LEN_TRIM(MeshDirName) > 0 ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(MeshDirName) // &
             '/' // ViewFactorsFile
       ELSE
         OutputName = ViewFactorsFile
       END IF

       IF(RadiationBody > 1) OutputName = OutputName//I2S(RadiationBody)

       INQUIRE(FILE=OutputName,EXIST=Found)
       IF(.NOT. Found) THEN
         CALL Warn(Caller,'View Factors File does NOT exist: '//TRIM(OutputName))
         Success = .FALSE.
         RETURN
       END IF

       BinaryMode = ListGetLogical( Params,'Viewfactor Binary Output',Found ) 
       IF(BinaryMode ) THEN
         SinglePrec = ListGetLogical( Params,'Viewfactor Single Precision',Found ) 
       ELSE
         SinglePrec = .FALSE.
       END IF
         
       IF( BinaryMode ) THEN
         CALL Info(Caller,'Loading view factors from binary file: '//OutputName,Level=5)

         OPEN( UNIT=VFUnit, FILE=OutputName, FORM = 'unformatted', &
             ACCESS = 'stream', STATUS='old', ACTION='read' )         
         READ( VFUnit ) n
         IF( n /= RadiationSurfaces ) THEN
           CALL Fatal(Caller,'Mismatch in viewfactor file size: '&
               //I2S(n)//' vs. '//I2S(RadiationSurfaces))
           Success = .FALSE.
           RETURN
         END IF
       ELSE
         CALL Info(Caller,'Loading view factors from ascii file: '//OutputName,Level=5)
         OPEN( VFUnit,File=OutputName )
       END IF

       ! Read in the ViewFactors
       DO i=1,RadiationSurfaces
         IF( BinaryMode ) THEN
           READ( VFUnit ) n
         ELSE
           READ( VFUnit,* ) n
         END IF

         IF(.NOT.ALLOCATED(ViewFactors(i) % Factors)) THEN
           ViewFactors(i) % NumberOfFactors = n
           ALLOCATE( ViewFactors(i) % Elements(n), ViewFactors(i) % Factors(n), STAT=istat )
           IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 6.')
         ELSE 
           n2 = SIZE( ViewFactors(i) % Factors) 
           IF(n /= n2) THEN
             DEALLOCATE(ViewFactors(i) % Factors, ViewFactors(i) % Elements)
             ALLOCATE( ViewFactors(i) % Factors(n), ViewFactors(i) % Elements(n), STAT=istat )
             IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 7.')
           END IF
           ViewFactors(i) % NumberOfFactors = n
         END IF

         Vals => ViewFactors(i) % Factors
         Cols => ViewFactors(i) % Elements

         Vals = 0; Cols = 0

         DO j=1,n
           IF( SinglePrec ) THEN
             READ(VFUnit) Cols(j),sval
             Vals(j) = sval        
           ELSE IF( BinaryMode ) THEN
             READ(VFUnit) Cols(j),Vals(j)         
           ELSE
             READ(VFUnit,*) t,Cols(j),Vals(j)         
           END IF
           Vals(j) = RelAreas(i) * Vals(j)  ! Scale by area to make symmetric
         END DO
       END DO
       CLOSE(VFUnit)
       
     END FUNCTION ReadViewFactorsFromFile



     ! This uses the view factors mainly to provide a simple test case where
     ! we assumes that emissivity is everywhere 1.
     !------------------------------------------------------------------------
     SUBROUTINE UseViewFactorsAsGebhartFactors()
       INTEGER :: i,j,n
       REAL(KIND=dp) :: s
       TYPE(Factors_t), POINTER :: GebhartFactors

       DO i=1,RadiationSurfaces
         j = ElementNumbers(i)
         Element => Model % Mesh % Elements(j)

         GebhartFactors => Element % BoundaryInfo % RadiationFactors
         IF ( .NOT. ASSOCIATED( GebhartFactors ) ) THEN
           ALLOCATE( GebhartFactors )
           Element % BoundaryInfo % RadiationFactors => GebhartFactors
         END IF

         n = ViewFactors(i) % NumberOfFactors
         GebhartFactors % NumberOfFactors = n
         GebhartFactors % NumberOfImplicitFactors = n
         ALLOCATE( GebhartFactors % Elements(n), GebhartFactors % Factors(n), STAT=istat)
         IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 20.')

         s = SUM(Viewfactors(i) % Factors)
         GebhartFactors % Factors = ViewFactors(i) % Factors/s
         GebhartFactors % Elements = ElementNumbers(Viewfactors(i) % Elements)
       END DO       
     END SUBROUTINE UseViewFactorsAsGebhartFactors


     SUBROUTINE CreateRadiationMatrix(n)
       INTEGER :: n

       INTEGER, ALLOCATABLE :: RowSpace(:),Reorder(:)
       INTEGER, POINTER :: Cols(:)
       LOGICAL :: AllocDone
       INTEGER :: i,j,t,previ,MatrixEntries

       ALLOCATE(RowSpace(n), Reorder(n))

       CALL Info('CreateRadiationMatrix','Creating matrix for Gebhart computation of size '//TRIM(I2S(n)),Level=12)

       UseFullMatrix = .FALSE.
       IF (IterSolveFactors) THEN
         UseFullMatrix = ListGetLogical( Params, 'Use Full Matrix for Radiation', Found )
       END IF
       
       ! Check whether the element already sees itself, it will when 
       ! gebhardt factors are computed. Also compute matrix size.
       RowSpace = 0
       DO i=1,n
         j = ViewFactors(i) % NumberOfFactors
         IF(j==0) CYCLE
         RowSpace(i) = j
         Cols => ViewFactors(i) % Elements
         IF (ALL(Cols/=i)) RowSpace(i) = RowSpace(i)+1
       END DO
       MatrixEntries = SUM(RowSpace(1:n))
       
       AllocDone = ASSOCIATED(G) .AND. .NOT. UseFullMatrix
       AllocDone = AllocDone .OR. ALLOCATED(G_full) .AND. UseFullMatrix

       IF(.NOT. AllocDone .OR. UpdateViewFactors) THEN
         IF (UseFullMatrix) THEN
           IF (ALLOCATED(G_full)) DEALLOCATE(G_full)
           ALLOCATE(G_full(n,n)); G_full = 0.0_dp
         ELSE
           IF(ASSOCIATED(G)) CALL FreeMatrix(G)

           ! Create matrix structures
           Reorder = [(i, i=1,n)]
           G => CRS_CreateMatrix(n,MatrixEntries,RowSpace,1,Reorder,.TRUE. )
 
           ! Create matrix entries
           DO t=1,n
             Cols => ViewFactors(t) % Elements         
             previ = G % Rows(t)-1
             DO j=1,ViewFactors(t) % NumberOfFactors
               CALL CRS_MakeMatrixIndex(G,t,Cols(j),previ)
             END DO
             CALL CRS_MakeMatrixIndex(G,t,t)
           END DO

           CALL CRS_SortMatrix(G)
           CALL CRS_ZeroMatrix(G)
           MatrixEntries = SIZE(G % Cols)
           WRITE(Message,'(A,T35,ES15.4)') 'View factors filling (%)',(100.0*MatrixEntries)/(n**2)
           CALL Info('RadiationFactors',Message,Level=5)
         END IF
         CALL Info('CreateRadiationMatrix','Number of entries in radiation matrix: '//I2S(MatrixEntries),Level=5)
       ELSE IF (UseFullMatrix) THEN
           G_full = 0._dp
       ELSE
          CALL CRS_ZeroMatrix(G)
       END IF
     END SUBROUTINE CreateRadiationMatrix
     

     SUBROUTINE TabulateSurfaceTemperatures(SurfT,T,Tperm)
       REAL(KIND=dp) :: SurfT(:), T(:)
       INTEGER :: Tperm(:)

       TYPE(Element_t), POINTER :: Element
       LOGICAL :: DG, Found
       INTEGER :: i,j,n,m
       REAL(KIND=dp) :: s
       INTEGER, POINTER :: Inds(:)
       INTEGER, TARGET :: DGInds(27)

       DG = ListGetLogical(Params, 'Discontinuous Galerkin',Found ) .OR. &
            ListGetLogical(Params, 'DG Reduced Basis',Found )

       DO i=1,RadiationSurfaces
         Element => Mesh % Elements(ElementNumbers(i))
         n = GetElementNOFNodes(Element)
         IF( DG ) THEN
           CALL DgRadiationIndexes(Element,n,DGInds,.TRUE.)
           Inds => DGInds(1:n)
         ELSE
           Inds => Element % NodeIndexes(1:n)
         END IF
         s = 0.0_dp; m = 0
         DO j=1,n
           IF (Tperm(Inds(j)) > 0) THEN
             s = s + T(Tperm(Inds(j)))
             m = m + 1
           END IF
         END DO
         IF (m > 0) SurfT(i) = s / m
       END DO
     END SUBROUTINE TabulateSurfaceTemperatures


     SUBROUTINE TabulateEmissivity()
       LOGICAL :: Found, SomeEmissivity0
       INTEGER :: i
       TYPE(Element_t), POINTER :: Element
       REAL(KIND=dp) :: Emissivity0

       
       CALL Info('TabulateEmissivity','Setting emissivities for radiation computation',Level=25)

       SomeEmissivity0 = .FALSE.
       IF(TopoCall ) THEN
         Emissivity0 = GetConstReal( Params,'Constant Emissivity',Found )
         IF(.NOT. Found) Emissivity0 = 0.5_dp
       END IF
       
       DO i=1,RadiationSurfaces         
         Element => Mesh % Elements(ElementNumbers(i))
         CALL GetElementEmissivity( Element, Emissivity(i), Absorptivity(i), &
             Reflectivity(i), Emissivity0, SomeEmissivity0 )
       END DO
       
       IF(SomeEmissivity0) THEN
         IF(FirstTime) THEN
           CALL Info('TabulateEmissivity','Using constant emissivity for some elements!',Level=6)
         END IF
         IF(.NOT. UpdateGebhartFactors ) THEN
           CALL Warn('TabulateEmissivity','Gebhart factors should be updated for non-constant emissivities!')
         END IF
       ELSE
         IF(FirstTime) THEN
           CALL Info('TabulateEmissivity','Using real emissivity for all elements!',Level=6)
         END IF
       END IF

     END SUBROUTINE TabulateEmissivity


     ! Emissivity, absorptivity and reflectivity of one boundary element.
     !-------------------------------------------------------------------
     SUBROUTINE GetElementEmissivity( Element, Emis, Abso, Refl, Emissivity0, SomeEmissivity0 )
       TYPE(Element_t), POINTER :: Element
       REAL(KIND=dp) :: Emis, Abso, Refl, Emissivity0
       LOGICAL :: SomeEmissivity0

       REAL(KIND=dp) :: Transmissivity
       LOGICAL :: Found, UseEmissivity0
       INTEGER :: n
       TYPE(ValueList_t), POINTER :: Vlist

       Vlist => GetEmissivityList( Element )

       UseEmissivity0 = .FALSE.
       IF( TopoCall ) THEN
         UseEmissivity0 = .NOT. ListCheckIsConstant( Vlist,'Emissivity' )
       END IF
              
       IF( UseEmissivity0 ) THEN
         Emis = ListGetConstReal( Vlist,'Initial Emissivity', Found )
         IF(.NOT. Found ) Emis = Emissivity0 
         Abso = Emis
         Refl = 1.0_dp - Abso 
         SomeEmissivity0 = .TRUE.
       ELSE          
         n = Element % TYPE % NumberOfNodes          
         CurrentModel % CurrentElement => Element
         Emis = SUM( ListGetReal( Vlist,'Emissivity',n,Element % NodeIndexes) ) / n
         Transmissivity= SUM( ListGetReal( Vlist,'Transmissivity',n,Element % NodeIndexes, Found) ) / n
         IF(.NOT. Found ) Abso = Emis
         Abso = SUM( ListGetReal( Vlist,'Absorptivity',n,Element % NodeIndexes, Found) ) / n
         IF(.NOT. Found ) Abso = Emis
         Refl = 1.0_dp - Abso - Transmissivity
       END IF
     END SUBROUTINE GetElementEmissivity


     ! The list where the emissivity of the boundary element is given: the BC or
     ! the material of a parent. Key is a unique index for the list.
     !------------------------------------------------------------------------------
     FUNCTION GetEmissivityList( Element, Key ) RESULT ( Vlist )
       TYPE(Element_t), POINTER :: Element
       INTEGER, OPTIONAL :: Key
       TYPE(ValueList_t), POINTER :: Vlist

       LOGICAL :: Found
       INTEGER :: k,bc_id,mat_id,i2,j2
       TYPE(Element_t), POINTER :: Parent

       DO bc_id=1,CurrentModel % NumberOfBCs
         IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
       END DO
       IF ( bc_id > CurrentModel % NumberOfBCs ) CALL Fatal('TabulateEmissivity','Could not find BC!')
         
       Vlist => CurrentModel % BCs(bc_id) % Values         
       IF( .NOT. ListCheckPresent(Vlist,'Emissivity') ) THEN
         DO k=1,2
           IF(k==1) THEN
             Parent => Element % BoundaryInfo % Left
           ELSE
             Parent => Element % BoundaryInfo % Right
           END IF
           IF(ASSOCIATED(Parent) ) THEN
             IF( Parent % BodyId > 0 .AND. Parent % BodyId <= CurrentModel % NumberOfBodies ) THEN
               mat_id = ListGetInteger( CurrentModel % Bodies(Parent % BodyId) % Values,'Material',Found)
               IF(Found) THEN
                 Vlist => CurrentModel % Materials(mat_id) % Values

                 IF(ListCheckPresent(Vlist,'Emissivity') ) THEN
                   IF( ASSOCIATED(Parent % DGIndexes) ) THEN
                     IF(.NOT. ASSOCIATED( Element % DGIndexes ) ) THEN
                       ALLOCATE( Element % DGIndexes(Element % TYPE % NumberOfNodes))
                       Element % DGIndexes = 0
                     END IF
                     DO i2 = 1, Element % TYPE % NumberOfNodes
                       DO j2 = 1, Parent % TYPE % NumberOfNodes
                         IF( Element % NodeIndexes(i2) == Parent % NodeIndexes(j2) ) THEN
                           Element % DGIndexes(i2) = Parent % DGIndexes(j2)
                           EXIT
                         END IF
                       END DO
                     END DO
                   END IF
                   EXIT
                 END IF
               END IF
             END IF
           END IF
         END DO
       END IF
       
       ! Radiation elements copied from other partitions have no parents but
       ! the body giving the emissivity.
       IF( .NOT. ( ASSOCIATED(Element % BoundaryInfo % Left) .OR. &
           ASSOCIATED(Element % BoundaryInfo % Right) ) ) THEN
         k = Element % BoundaryInfo % EmissivityBody
         IF( k > 0 .AND. .NOT. ListCheckPresent(Vlist,'Emissivity') ) THEN
           mat_id = ListGetInteger( CurrentModel % Bodies(k) % Values,'Material',Found)
           IF( Found ) Vlist => CurrentModel % Materials(mat_id) % Values
         END IF
       END IF

       IF(.NOT. ASSOCIATED(Vlist) ) CALL Fatal('TabulateEmissivity','Emissivity list not associated!')

       IF( PRESENT(Key) ) THEN
         IF( ASSOCIATED(Vlist, CurrentModel % BCs(bc_id) % Values) ) THEN
           Key = bc_id
         ELSE
           Key = CurrentModel % NumberOfBCs + mat_id
         END IF
       END IF
     END FUNCTION GetEmissivityList


     !------------------------------------------------------------------------------
     ! To save some time tabulate the spectral emissivity data for each temperature.
     !------------------------------------------------------------------------------
     SUBROUTINE TabulateSpectralEmissivity(Emissivity,Absorptivity,Trad,IsRadiator,SimpleTDep)
       REAL(KIND=dp) :: Trad
       REAL(KIND=dp) :: Emissivity(:)
       REAL(KIND=dp) :: Absorptivity(:)
       LOGICAL :: IsRadiator, SimpleTdep
              
       REAL(KIND=dp), ALLOCATABLE :: SaveValues(:)
       TYPE(Variable_t), POINTER :: TVar
       TYPE(ValueList_t), POINTER :: Vlist
       TYPE(Element_t), POINTER :: Element
       INTEGER :: i

       CALL Info('TabulateSpectralEmissivity','Precomputing emissivities for faster radiosity computation',Level=5)       

       ! If we have simple dependence only (dependence just on temperature) we can call it through
       ! a simplefied function call. Otherwise we overwrite the current temperature and use the generic
       ! ListGetReal function, and then rewert back to original temperature.
       IF(.NOT. SimpleTdep ) THEN       
         TVar => VariableGet(Mesh % Variables,'Temperature')
         ALLOCATE( SaveValues(SIZE(TVar % Values) ) )
         SaveValues = TVar % Values
         TVar % Values = Trad
       END IF
                
       DO i=1,RadiationSurfaces         
         Element => Mesh % Elements(ElementNumbers(i))
         Vlist => GetEmissivityList( Element )
         CALL GetSpectralEmissivity( Element, Vlist, Trad, IsRadiator, SimpleTdep, &
             Emissivity(i), Absorptivity(i) )
       END DO

       IF(.NOT. SimpleTdep ) THEN
         TVar % Values = SaveValues
         DEALLOCATE(SaveValues)
       END IF
                
     END SUBROUTINE TabulateSpectralEmissivity


     ! Emissivity and absorptivity of a boundary element for radiation of temperature Trad.
     ! Unless SimpleTdep the temperature field must have been temporarily set to Trad.
     !--------------------------------------------------------------------------------------
     SUBROUTINE GetSpectralEmissivity( Element, Vlist, Trad, IsRadiator, SimpleTdep, Emis, Abso )
       TYPE(Element_t), POINTER :: Element
       TYPE(ValueList_t), POINTER :: Vlist
       REAL(KIND=dp) :: Trad, Emis, Abso
       LOGICAL :: IsRadiator, SimpleTdep

       LOGICAL :: Found
       INTEGER :: n

       IF( SimpleTdep ) THEN
         Emis = ListGetFun( Vlist,'Emissivity',Trad,minv=0.0_dp,maxv=1.0_dp)
         Found = .FALSE.
         IF(IsRadiator) THEN
           Abso = ListGetFun( VList,'Radiator Absorptivity',Trad,Found,minv=0.0_dp,maxv=1.0_dp)
         END IF
         IF(.NOT. Found ) Abso = ListGetFun( VList,'Absorptivity',Trad,Found,minv=0.0_dp,maxv=1.0_dp)         
         IF(.NOT. Found ) Abso = Emis
       ELSE          
         n = Element % TYPE % NumberOfNodes          
         CurrentModel % CurrentElement => Element
         Emis = SUM( ListGetReal( Vlist,'Emissivity',n,Element % NodeIndexes) ) / n
         Found = .FALSE.
         IF(IsRadiator) THEN
           Abso = SUM( ListGetReal( Vlist,'Radiator Absorptivity',n,Element % NodeIndexes, Found) ) / n
         END IF
         IF(.NOT. Found ) Abso = SUM( ListGetReal( Vlist,'Absorptivity',n,Element % NodeIndexes, Found) ) / n
         IF(.NOT. Found ) Abso = Emis
       END IF
     END SUBROUTINE GetSpectralEmissivity
            

     SUBROUTINE CalculateRadiation()

       INTEGER :: istat

       !IF(Radiosity .AND. FirstTime) RETURN

       !CALL Info(Caller,'Computing factors...',Level=5)

       IF(FirstTime) CALL InitRadiationSolver(TSolver,Solver)
       CALL CreateRadiationMatrix(RadiationSurfaces)

       ALLOCATE(Emissivity(RadiationSurfaces), Reflectivity(RadiationSurfaces), &
           Absorptivity(RadiationSurfaces), STAT=istat)
       IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 10.')
       IF( Hybrid ) THEN
         CALL TabulateHybrid()
       ELSE
         CALL TabulateEmissivity()
         CALL ParallelSurfaceData( Emissivity, Absorptivity, Reflectivity )
       END IF

       IF( Radiosity ) THEN
         CALL CalculateRadiosity()
       ELSE
         ! Fill the matrix for gebhardt factors
         CALL CalculateGebhartFactors()
         IF (UseFullMatrix) THEN
           DEALLOCATE(G_full)
         ELSE
           CALL FreeMatrix(G);G => NULL()
         ENDIF
       END IF

       DEALLOCATE(Emissivity,Reflectivity,Absorptivity)
       IF( Hybrid ) DEALLOCATE( FineElem, FineEmis, FineAbs, FineT, CoarseT, &
           FineArea, FineGrad, FineW, CoarseA, CoarseS, CoarseScale, PcFine, PcRad, PcW )
       
     END SUBROUTINE CalculateRadiation


     ! Calculate Gebhart factors for radiation. These result to good convergence with the cost of
     ! computing many smaller linear systems.
     !--------------------------------------------------------------------------------------------
     SUBROUTINE CalculateGebhartFactors()

       REAL(KIND=dp) :: MinFactor, MaxOmittedFactor, ConsideredSum
       INTEGER :: Colj,i,j,k,n,t, ImplicitEntries, MatrixEntries, previ
       LOGICAL :: gTriv, ImplicitLimitIs
       LOGICAL, ALLOCATABLE :: Gray(:)
       REAL(KIND=dp) :: r,s,st,PrevSelf, MinSum,MaxSum,SolSum,FactorSum,ImplicitSum,&
           ImplicitLimit, NeglectLimit

       REAL(KIND=dp), POINTER :: Vals(:)
       INTEGER, POINTER :: Cols(:)
       INTEGER, ALLOCATABLE :: FacPerm(:)
       REAL(KIND=dp), ALLOCATABLE :: Fac(:), RowSums(:), RHS(:), SOL(:), Diag(:)

       TYPE(Factors_t), POINTER :: GebhartFactors

       n = RadiationSurfaces
       ALLOCATE(FacPerm(n), Fac(n), RHS(n), SOL(n), Diag(n))
       RHS = 0.0_dp

       MaxOmittedFactor = 0._dp
       MatrixEntries = 0
       ImplicitEntries = 0
       
       MinFactor = GetConstReal( Params, 'Minimum Gebhart Factor',Found )
       IF (.NOT. Found) &
           MinFactor = GetConstReal( Params, 'Minimum Gebhardt Factor',Found )
       IF(.NOT. Found) MinFactor = 1.0d-20

       ! If all surfaces are black the Gebhart factors are directly given by the view factors
       gTriv = ALL(ABS(Reflectivity)<=AEPS)
       ALLOCATE(Gray(RadiationSurfaces))
       Gray = ABS(Reflectivity) > AEPS


       ImplicitLimit = GetConstReal( Params, 'Implicit Gebhart Factor Fraction', ImplicitLimitIs) 
       IF  (.NOT. ImplicitLimitIs) &
           ImplicitLimit = GetConstReal( Params, 'Implicit Gebhardt Factor Fraction', ImplicitLimitIs) 

       NeglectLimit  = GetConstReal( Params, 'Neglected Gebhart Factor Fraction', Found) 
       IF(.NOT.Found) &
           NeglectLimit  = GetConstReal( Params, 'Neglected Gebhardt Factor Fraction', Found) 
       IF(.NOT. Found) NeglectLimit = 1.0d-6

       ! The equation for the Gebhart factors from surface t is (A-R*AF)x = e_t, with
       ! R=diag(Reflectivity), AF_ij = A_i*F_ij symmetric and Fac = e_t*E*AF*x. For
       ! black surfaces (R_i=0) the row reduces to A_i x_i = e_t(i) and the column is
       ! zero in R*AF. So black unknowns are known a priori and are moved to the RHS.
       ! Dividing the gray rows by R_i the system becomes (A/R-AF)x = e_t/R for gray,
       ! and A x = e_t for black surfaces, which is symmetric & diagonally dominant.
       IF(.NOT. gTriv) THEN
         DO i=1,RadiationSurfaces
           Vals => ViewFactors(i) % Factors
           Cols => ViewFactors(i) % Elements

           IF( .NOT. UseFullMatrix ) previ = G % Rows(i)-1
           IF( Gray(i) ) THEN
             DO j=1,ViewFactors(i) % NumberOfFactors
               IF( .NOT. Gray(Cols(j)) ) CYCLE
               IF (UseFullMatrix) THEN
                 G_full(i,Cols(j)) = G_full(i,Cols(j)) - Vals(j)
               ELSE
                 CALL CRS_AddToMatrixElement(G,i,Cols(j),-Vals(j),previ)
               END  IF
             END DO
             Diag(i) = RelAreas(i) / Reflectivity(i)
           ELSE
             Diag(i) = RelAreas(i)
           END IF
           IF (UseFullMatrix) THEN
             G_full(i,i) = Diag(i)
           ELSE
             CALL CRS_AddToMatrixElement( G,i,i,Diag(i) )
           END IF
         END DO
       END IF
       
       ! Scale matrix to unit diagonals
       Diag = SQRT(1._dp/MAX(ABS(Diag),1.0d-12))
       DO i=1,RadiationSurfaces
         IF ( UseFullMatrix ) THEN
           DO j=1,RadiationSurfaces
             G_full(i,j) = G_full(i,j)*Diag(i)*Diag(j)
           END DO
         ELSE
           DO j=G % Rows(i),G % Rows(i+1)-1
             G % Values(j) = G % Values(j)*Diag(i)*Diag(G % Cols(j))
           END DO
         END IF
       END DO
       
       SOL = 1.0d-4
       st = RealTime()

       n = 0
       ALLOCATE(RowSums(RadiationSurfaces), STAT=istat)
       IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 12.')       
       RowSums=0

       DO t=1,RadiationSurfaces
           
         i = ElementNumbers(t)
         Element => Mesh % Elements(i)

         IF ( gTriv ) THEN

           Vals => ViewFactors(t) % Factors
           Cols => ViewFactors(t) % Elements
           Fac = 0._dp
           DO k=1,ViewFactors(t) % NumberOfFactors
             Fac(Cols(k)) = Vals(k) / RelAreas(t)
           END DO

           DO i=1,RadiationSurfaces
             s = Fac(i)*RelAreas(t)/(RelAreas(i)*Emissivity(i))
             IF (s>MinFactor ) n=n+1
             RowSums(i) = RowSums(i) + s
           END DO

         ELSE
           ! RHS of the unit diagonal scaled system, see above
           RHS = 0.0_dp
           IF( Gray(t) ) THEN
             RHS(t) = Diag(t) / Reflectivity(t)
           ELSE
             RHS(t) = Diag(t)
             Vals => ViewFactors(t) % Factors
             Cols => ViewFactors(t) % Elements
             DO k=1,ViewFactors(t) % NumberOfFactors
               j = Cols(k)
               IF( Gray(j) ) RHS(j) = RHS(j) + Diag(j) * Vals(k) / RelAreas(t)
             END DO
           END IF

           ! It may be a good initial start that the Gii is
           ! the same as previously
           IF (t>1) THEN
             PrevSelf = SOL(t)
             SOL(t) = SOL(t-1)
             SOL(t-1) = PrevSelf
           END IF

           SOL = SOL/Diag

           IF ( Element % PartIndex==ParEnv % MyPE )THEN
             IF(IterSolveFactors) THEN
               Solver % Matrix => G
               IF(UseFullMatrix) THEN
                 BLOCK
                   TYPE(Matrix_t), POINTER :: Gm
                   Gm => AllocateMatrix()
                   Gm % NumberOfRows = RadiationSurfaces
                   fm_G => G_full
                   CALL IterSolver( Gm, SOL, RHS, Solver, MatVecF=AddrFunc(fm_Matvec) )
                   DEALLOCATE(Gm)
                 END BLOCK
               ELSE
                 CALL IterSolver( G, SOL, RHS, Solver )
               END IF
             !------------------------------------------------------------------------------
             ELSE           
               IF (t==1) THEN
                 CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .TRUE. )
                 CALL ListAddLogical( Solver % Values, 'Linear System Free Factorization', .FALSE. )
               ELSE IF(t==2) THEN
                 CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .FALSE. )
               END IF

               CALL DirectSolver( G, SOL, RHS, Solver )
             END IF

             SOL = SOL*Diag
             CALL ListRemove(Solver % Values,'Linear System Free Factorization')
           ELSE
             SOL = 0
           END IF

           n = 0
           DO i=1,RadiationSurfaces
             Vals => ViewFactors(i) % Factors
             Cols => ViewFactors(i) % Elements

             s = 0.0_dp
             DO k=1,ViewFactors(i) % NumberOfFactors
               s = s + Vals(k)*SOL(Cols(k))
             END DO
             Fac(i) = s*Emissivity(t)*Emissivity(i)

             ! rowsums should add up to 1
             s = Fac(i)*RelAreas(t)/(RelAreas(i)*Emissivity(i))
             IF (s>MinFactor ) n=n+1
             RowSums(i) = RowSums(i) + s
           END  DO
         END IF

         FactorSum = SUM(Fac)
         ConsideredSum = 0.0_dp

         IF(ImplicitLimitIs) THEN
           DO i=1,RadiationSurfaces
             FacPerm(i) = i
           END DO
           ! Ensure that the self vision is always implicit to avoid trouble in the future!
           Fac(t) = Fac(t) + FactorSum
           CALL SortR( RadiationSurfaces, FacPerm, Fac) 
           Fac(1) = Fac(1) - FactorSum        

           ConsideredSum = 0.0_dp
           n = 0
           DO i=1,RadiationSurfaces
             IF(ConsideredSum < (1.0_dp-NeglectLimit) * FactorSum) THEN
               ConsideredSum = ConsideredSum + Fac(i)
               n = i
             END IF
           END DO
         ELSE
           n = 0
           DO i=1,RadiationSurfaces
             IF ( Fac(i) > MinFactor ) n = n + 1
           END DO
         END IF

         MatrixEntries = MatrixEntries + n
         Element => Mesh % Elements(ElementNumbers(t))
         GebhartFactors => Element % BoundaryInfo % RadiationFactors
         IF ( .NOT. ASSOCIATED( GebhartFactors ) ) THEN
           ALLOCATE( GebhartFactors )
           Element % BoundaryInfo % RadiationFactors => GebhartFactors
         END IF

         IF(FirstTime) THEN
           GebhartFactors % NumberOfFactors = n
           GebhartFactors % NumberOfImplicitFactors = n
           ALLOCATE( GebhartFactors % Elements(n), GebhartFactors % Factors(n), STAT=istat)
           IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 13.')
         ELSE IF(ImplicitLimitIs) THEN 
           IF( TopologyFixed ) THEN
             CALL Warn(Caller,'Matrix topology cannot be fixed with implicit Gebhart factors')
           END IF
           TopologyFixed = .FALSE.
           TopologyTest = .FALSE.
           DEALLOCATE( GebhartFactors % Elements, GebhartFactors % Factors )
           GebhartFactors % NumberOfFactors = n
           ALLOCATE( GebhartFactors % Elements(n), GebhartFactors % Factors(n), STAT=istat )
           IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 14.')
           GebhartFactors % NumberOfImplicitFactors = 0
         ELSE IF(GebhartFactors % NumberOfFactors /= n .AND. .NOT. TopologyFixed) THEN         
           TopologyTest = .FALSE.
           DEALLOCATE( GebhartFactors % Elements, GebhartFactors % Factors )
           GebhartFactors % NumberOfFactors = n         
           ALLOCATE( GebhartFactors % Elements(n), GebhartFactors % Factors(n), STAT=istat )
           IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 15.')

           GebhartFactors % NumberOfImplicitFactors = n
         END IF

         Vals => GebhartFactors % Factors
         Cols => GebhartFactors % Elements

         IF( ImplicitLimitIs ) THEN

           ImplicitSum = 0.0d0
           DO i=1,n
             Cols(i) = ElementNumbers(FacPerm(i)) 
             Vals(i) = Fac(i)

             IF(ImplicitSum < ImplicitLimit * FactorSum) THEN
               ImplicitSum = ImplicitSum + Fac(i)
               GebhartFactors % NumberOfImplicitFactors = i
             END IF
           END DO
           Vals(2:n) = Vals(2:n) * (FactorSum - Vals(1)) / (ConsideredSum - Vals(1))

           IF(ImplicitLimit < TINY(ImplicitLimit)) GebhartFactors % NumberOfImplicitFactors = 0       
           ImplicitEntries = ImplicitEntries + GebhartFactors % NumberOfImplicitFactors

         ELSE IF(FirstTime .OR. .NOT. TopologyFixed) THEN
           n = 0
           DO i=1,RadiationSurfaces
             IF ( Fac(i) > MinFactor ) THEN
               n = n + 1
               IF(TopologyTest .AND. Cols(n) /= ElementNumbers(i)) TopologyTest = .FALSE.
               Cols(n) = ElementNumbers(i) 
               Vals(n) = Fac(i)
               ConsideredSum = ConsideredSum + Fac(i)
             END IF
           END DO
         ELSE
           ! If the topology is fixed the values are put only according to the existing structure
           ! and others are neglected
           n = GebhartFactors % NumberOfFactors         
           Vals => GebhartFactors % Factors
           Cols => GebhartFactors % Elements

           DO i=1,n
             j = InvElementNumbers(Cols(i)-nBulk)
             Vals(i) = Fac(j)
             ConsideredSum = ConsideredSum + Fac(j)
           END DO
         END IF

         MaxOmittedFactor = MAX(MaxOmittedFactor,(FactorSum-ConsideredSum)/FactorSum) 

         IF ( RealTime() - st > 10.0 ) THEN
           WRITE(Message,'(A,I3,A)' ) '   Solution: ', &
               INT((100.0*t)/RadiationSurfaces),' % done'
           CALL Info( Caller, Message, Level=5 )
           st = RealTime()
         END IF
       END  DO

       MinSum = MINVAL(RowSums)
       MaxSum = MAXVAL(RowSums)

       WRITE(Message,'(A,T35,2ES15.6)') 'Minimum Gebhart factors sum',MINVAL(RowSums)
       CALL Info(Caller,Message,Level=5)
       WRITE(Message,'(A,T35,2ES15.6)') 'Maximum Gebhart factors sum',MAXVAL(RowSums)
       CALL Info(Caller,Message,Level=5)
       WRITE(Message,'(A,T35,ES15.6)') 'Maximum share of omitted factors',MaxOmittedFactor
       CALL Info(Caller,Message,Level=5)
       WRITE(Message,'(A,T35,ES15.6)') 'Gebhart factors filling (%)',(100.0 * MatrixEntries) / &
           (RadiationSurfaces**2)
       CALL Info(Caller,Message,Level=5)
       WRITE(Message,'(A,T38,I0)') 'Gebhart factors count',MatrixEntries
       CALL Info(Caller,Message,Level=5)
       IF(ImplicitEntries > 0) THEN
         WRITE(Message,'(A,T38,I0)') 'Implicit factors count',ImplicitEntries
         CALL Info(Caller,Message,Level=5)
         WRITE(Message,'(A,T35,ES15.6)') 'Implicit factors filling (%)',(100.0 * ImplicitEntries) / &
             (RadiationSurfaces**2)
         CALL Info(Caller,Message,Level=5)
       END IF

       IF(SaveFactors) THEN
         CALL SaveGebhartFactors()       
       END IF
       DEALLOCATE(RowSums)
     END SUBROUTINE CalculateGebhartFactors

     
     ! When Gebhart factors may have changed also modify the matrix topology so that
     ! when we assemble the matrices we are not hitting non-existing entries.
     !--------------------------------------------------------------------------------
     SUBROUTINE UpdateMatrixTopologyWithFactors()

       TYPE(Matrix_t), POINTER :: AMatrix
       LOGICAL :: OptimizeBW, UseGiven, Found
       INTEGER :: j,n,MatrixFormat
       INTEGER, POINTER :: NewPerm(:), TempPerm(:)
     
       CALL Info(Caller,'Recreating the matrix structure for radiation',Level=5)

       MatrixFormat = Tsolver % Matrix % FORMAT

       ! We have different default here!
       OptimizeBW = ListGetLogical(TSolver % Values,'Optimize Bandwidth',Found) 
       IF(.NOT. Found) OptimizeBW = .FALSE.

       ! If we do not use the optimized, we use the previous Perm (which could be optimized as well)
       UseGiven = .NOT. OptimizeBW

       CALL FreeMatrix( TSolver % Matrix)         

       IF ( OptimizeBW ) THEN
         CALL Info(Caller,'Creating new matrix topology')
         ALLOCATE( NewPerm( SIZE(Tsolver % Variable % Perm)), STAT=istat)
         IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 15.')
         TempPerm => Tsolver % Variable % Perm         
       ELSE
         CALL Info(Caller,'Using existing matrix topology')
         NewPerm => Tsolver % Variable % Perm
       END IF
       
       AMatrix => CreateMatrix( CurrentModel,TSolver,TSolver % Mesh, &
           NewPerm, 1, MatrixFormat, OptimizeBW,  &
           ListGetString( TSolver % Values, 'Equation', Found ), UseGivenPerm = UseGiven )       
              
       ! Reorder the primary variable for bandwidth optimization:
       ! --------------------------------------------------------
       IF ( OptimizeBW ) THEN
         WHERE( NewPerm > 0 )
           TSolver % Variable % Values( NewPerm ) = &
               TSolver % Variable % Values( TempPerm )
         END WHERE

         IF ( ASSOCIATED( TSolver % Variable % PrevValues ) ) THEN
           DO j=1,SIZE( TSolver % Variable % PrevValues,2 )
             WHERE( NewPerm > 0 )
               TSolver % Variable % PrevValues( NewPerm,j) = &
                   TSolver % Variable % PrevValues(TempPerm,j)
             END WHERE
           END DO
         END IF

         BLOCK
           TYPE(Variable_t), POINTER :: ExpVar
           CHARACTER(LEN=MAX_NAME_LEN) :: str
           INTEGER, POINTER :: ExpPerm(:)
           INTEGER :: k
           NULLIFY(ExpPerm)         
           DO j=1,10
             str = ListGetString(TSolver % Values,'exported variable '//I2S(j),Found)
             IF(.NOT. Found) EXIT
             ExpVar => VariableGet(TSolver % Mesh % Variables, str, ThisOnly = .TRUE. )             
             IF(ASSOCIATED(ExpVar)) THEN
               IF(ASSOCIATED(ExpVar % Perm, TSolver % Variable % Perm ) ) THEN
                 DO k=1,ExpVar % Dofs
                   WHERE( NewPerm > 0 )
                     ExpVar % Values( ExpVar % Dofs*(NewPerm-1)+k) = &
                         ExpVar % Values( ExpVar % Dofs*(TempPerm-1)+k)
                   END WHERE
                 END DO
               END IF
             END IF
           END DO
         END BLOCK
                  
         Tsolver % Variable % Perm = NewPerm
         DEALLOCATE( NewPerm )
       END IF

       ! TODO: CreateMatrix should do these:
       ! -----------------------------------
       AMatrix % Lumped = GetLogical( Params, 'Lumped Mass Matrix', Found )
       AMatrix % Symmetric = ListGetLogical( Params, 'Linear System Symmetric', Found )       

       n = AMatrix % NumberOFRows
       ALLOCATE( AMatrix % RHS(n), STAT=istat)
       IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 16.')
       
       ! Transient case additional allocations:
       ! --------------------------------------
       IF ( ListGetString( CurrentModel % Simulation,'Simulation Type' ) == 'transient' ) THEN
         ALLOCATE( Amatrix % Force(n, TSolver % TimeOrder+1), STAT=istat )
         IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error 17.')         
         Amatrix % Force = 0.0d0
       END IF

       TSolver % Matrix => Amatrix
       CALL ParallelInitMatrix( TSolver, AMatrix )
       
     END SUBROUTINE UpdateMatrixTopologyWithFactors


     ! Compute radiosities, i.e.
     !
     ! Solve (r*F-1)J = -e*Sigma*T^4 - r*R*P for J
     !  J: Radiosity
     !  e: Emissivity
     !  r: Reflectivity (=1-Emissivity)
     !  F: The "viewfactor" matrix
     !  Sigma: Stefan-Boltzmann constant
     !  T: Temperature
     !  P: "Radiator" power
     !  R: Elementwise visibility of the "radiators"
     SUBROUTINE CalculateRadiosity()

       REAL(KIND=dp), POINTER :: Temperature(:)
       INTEGER, POINTER  :: TempPerm(:)
       REAL(KIND=dp), ALLOCATABLE :: SurfaceTemperature(:)

       Temperature => Null()
       IF(ASSOCIATED(TSolver % Variable))  THEN
         TempPerm => TSolver % Variable % Perm
         Temperature => TSolver % Variable % Values
       END IF

       IF(.NOT.ASSOCIATED(Temperature)) &
         CALL Fatal(Caller, &
              "Radiosity solution can't be completed without the temperature field.")

       Sigma = ListGetConstReal( Model % Constants,&
         'Stefan Boltzmann',UnfoundFatal=.TRUE. )

       ALLOCATE(SurfaceTemperature(RadiationSurfaces))
       IF( Hybrid ) THEN
         SurfaceTemperature = CoarseT
       ELSE
         SurfaceTemperature = 0.0_dp
         CALL TabulateSurfaceTemperatures(SurfaceTemperature,Temperature,TempPerm)
         CALL ParallelSurfaceData( SurfaceTemperature )
       END IF

       IF( InfoActive(30) ) THEN
         PRINT *,'Temp range:',MINVAL(SurfaceTemperature),MAXVAL(SurfaceTemperature)
         PRINT *,'Emis range:',MINVAL(Emissivity),MAXVAL(Emissivity)
         PRINT *,'Abs range:',MINVAL(Absorptivity),MAXVAL(Absorptivity)
       END IF

       IF( Spectral ) THEN
         IF( Hybrid ) THEN
           CALL SpectralHybrid()
         ELSE
           CALL SpectralRadiosity(SurfaceTemperature)
         END IF
       ELSE
         CALL ConstantRadiosity(SurfaceTemperature)
       END IF
     END SUBROUTINE CalculateRadiosity
       

     ! Compute radiosity vector in the case that the emissivity is constant with temperature
     !--------------------------------------------------------------------------------------
     SUBROUTINE ConstantRadiosity(SurfaceTemperature)
       REAL(KIND=dp) :: SurfaceTemperature(:)
 
       LOGICAL :: RBC
       INTEGER :: i,j,k
       REAL(KIND=dp) :: r, e, a, c, Temp, Black
       REAL(KIND=dp), ALLOCATABLE :: RadiatorPowers(:), &
            RHS(:),RHS_d(:),SOL(:),SOL_d(:), Diag(:), AG(:), AG_d(:), Rdir(:)

       ALLOCATE(RHS(RadiationSurfaces),SOL(RadiationSurfaces),Diag(RadiationSurfaces), &
           AG(RadiationSurfaces),Rdir(RadiationSurfaces))
       RHS = 0.0_dp
       Rdir = 0.0_dp

       IF (Newton) THEN
         ALLOCATE( RHS_d(RadiationSurfaces), SOL_d(RadiationSurfaces) )
         ALLOCATE( AG_d(RadiationSurfaces) )
         RHS_d = 0.0_dp
       END IF

       ! Assemble the equations, first coefficient matrix:
       ! -------------------------------------------------
       CALL RadiosityAssembly(RadiationSurfaces,G,Diag)

       ! ... and then the RHS:
       ! ---------------------
       DO i=1,RadiationSurfaces
         e = Emissivity(i)
         a = Absorptivity(i)
         r = 1-a  ! 1-e
         c = RelAreas(i) / a
         Temp = SurfaceTemperature(i)
         Black = Sigma*Temp**4
         RHS(i) = -c*e*Black
         IF(Newton) RHS_d(i) = RHS(i)*(4/Temp)
       END DO

       ! Check for radiation sources:
       RBC = CheckForRadiators(RadiatorPowers)
       IF( RBC .AND. Hybrid ) THEN
         ! Direct irradiation of fine elements; its reflected part is the source of
         ! the coarse element: S = sum(A*r*Grad)/A.
         DO k=1,nFine
           Element => FineMesh % Elements(FineElem(k))
           IF ( .NOT. ALLOCATED(Element % BoundaryInfo % Radiators)) CYCLE
           FineGrad(k) = SUM(Element % BoundaryInfo % Radiators*RadiatorPowers)
         END DO
         DO k=1,nPc
           j = PcFine(k)
           i = PcRad(k)
           CoarseS(i) = CoarseS(i) + PcW(k) * FineW(j) * FineArea(j) * (1-FineAbs(j)) * FineGrad(j)
         END DO
         CALL ParallelSum(CoarseS)
         DO i=1,RadiationSurfaces
           CoarseS(i) = CoarseS(i) / Areas(i)
           a = Absorptivity(i)
           r = 1-a
           c = RelAreas(i) / a
           RHS(i) = RHS(i) - c * CoarseS(i)
         END DO
       ELSE IF( RBC) THEN
         DO i=1,RadiationSurfaces
           Element => Mesh % Elements(ElementNumbers(i))
           IF ( ALLOCATED(Element % BoundaryInfo % Radiators)) THEN
             e = Emissivity(i)
             a = Absorptivity(i)
             !r = Reflectivity(i)
             r = 1-a  ! e
             c = RelAreas(i) / a
             Rdir(i) = SUM(Element % BoundaryInfo % Radiators*RadiatorPowers)
             RHS(i) = RHS(i) - c*r*Rdir(i)
           END IF
         END DO
       END IF

       CALL BlackRadiosityToRHS(RadiationSurfaces,RHS)
       IF( Newton ) CALL BlackRadiosityToRHS(RadiationSurfaces,RHS_d)

       ! Solve for the radiosities and their derivatives with respect
       ! to the temperature
       !-------------------------------------------------------------
       CALL RadiationLinearSolver(RadiationSurfaces,G,SOL,RHS,Diag,Solver)
       IF( Newton ) THEN
         CALL RadiationLinearSolver(RadiationSurfaces,G,SOL_d,RHS_d, &
                      Diag, Solver, Scaling=.FALSE.)
       END IF

       ! Store the results for access by e.g. heat equation solvers:
       !------------------------------------------------------------
       IF( Hybrid ) THEN
         IF(Newton) THEN
           CALL UpdateHybridFactors(SOL,SOL_d)
         ELSE
           CALL UpdateHybridFactors(SOL)
         END IF
       ELSE
         ! The heat equation needs the absorbed irradiation from the surfaces
         ! and directly from the radiators.
         CALL AbsorbedIrradiation(RadiationSurfaces,SOL,AG)
         AG = AG + Absorptivity(1:RadiationSurfaces) * Rdir
         IF(Newton) THEN
           CALL AbsorbedIrradiation(RadiationSurfaces,SOL_d,AG_d)
           CALL UpdateRadiosityFactors(AG,AG_d)
         ELSE
           CALL UpdateRadiosityFactors(AG)
         END IF
       END IF
     END SUBROUTINE ConstantRadiosity
       
     
     ! Divide temperature into intervals.
     ! This is needed in order to compute problems where emissivity depends on temperature.
     !-------------------------------------------------------------------------------------
     SUBROUTINE SpectralRadiosity(SurfaceTemperature)
       REAL(KIND=dp) :: SurfaceTemperature(:)

       REAL(KIND=dp) :: Tmin, Tmax, dT, Trad
       INTEGER :: i,j,k,kmin,kmax
       REAL(KIND=dp) :: q, qsum, totsum, c, r, e, a, s, Temp, Black

       LOGICAL :: RBC, ApproxNewton, AccurateNewton, UsedEdT, SimpleTdep
       INTEGER, ALLOCATABLE :: RadiatorSet(:)
       REAL(KIND=dp), ALLOCATABLE :: RadiatorPowers(:), RadiatorTemps(:), &
            RHS(:),RHS_d(:),SOL(:),SOL_d(:), Diag(:)
       REAL(KIND=dp), ALLOCATABLE :: tmpSOL(:), tmpSOL_d(:), EffAbs(:), EffTemp(:), &
            AG(:), AG_d(:), Rdir(:)

       ALLOCATE(RHS(RadiationSurfaces),SOL(RadiationSurfaces),Diag(RadiationSurfaces))
       RHS = 0.0_dp
       IF (Newton) THEN
         ALLOCATE( RHS_d(RadiationSurfaces), SOL_d(RadiationSurfaces) )
         RHS_d = 0.0_dp
       END IF

       ApproxNewton = .FALSE.
       AccurateNewton = .FALSE.      
       IF( Newton ) THEN
         AccurateNewton = ListGetLogical( TSolver % Values,'Accurate Spectral Newton',Found ) 
         ApproxNewton = .NOT. AccurateNewton
       END IF

       SimpleTdep = ListGetLogical( TSolver % Values,'Radiosity Simple Temperature Dependence',Found)
       
       Tmin = MINVAL(SurfaceTemperature)
       Tmax = MAXVAL(SurfaceTemperature)
       
       WRITE(Message,'(A,ES12.3)') 'Minimum boundary temperature: ',Tmin
       CALL Info('SpectralRadiosity',Message,Level=10)
       WRITE(Message,'(A,ES12.3)') 'Maximum boundary temperature: ',Tmax
       CALL Info('SpectralRadiosity',Message,Level=10)
       
       IF( Tmin < 0.0_dp ) THEN
         CALL Fatal('SpectralRadiosity','Negative temperature not a good starting point!')
       END IF
       
       ! We have a fixed dT instead of having variable one related to Tmin and Tmax since
       ! adaptive intervals could generate funny attractors.
       dT = ListGetCReal( TSolver % Values,'Spectral dT',UnfoundFatal=.TRUE.) 
       
       kmin = FLOOR( Tmin / dT )
       kmax = CEILING( Tmax / dT )

       CALL Info('SpectralRadiosity','Going through discrete intervals: '&
           //I2S(kmin)//'-'//I2S(kmax))

       SOL = 0.0_dp
       IF(Newton) SOL_d = 0.0_dp
       
       ALLOCATE( tmpSOL(RadiationSurfaces), AG(RadiationSurfaces), Rdir(RadiationSurfaces) )
       IF(Newton) ALLOCATE(tmpSOL_d(RadiationSurfaces), AG_d(RadiationSurfaces))

       ALLOCATE(EffAbs(RadiationSurfaces),EffTemp(RadiationSurfaces))
       EffAbs = 0.0_dp
       EffTemp = 0.0_dp
       
       totsum = 0.0_dp

       DO k = kmin, kmax

         qsum = 0.0_dp         
         DO i=1,RadiationSurfaces           
           q = ( SurfaceTemperature(i) / dT - k )
           IF( ABS(q) < 1 ) THEN
             q = 1-ABS(q)
             qsum = qsum + q
           END IF
         END DO

         ! There is nothing to compute here
         ! So no need to resolve equations for this interval.        
         IF(qsum < 1.0d-6 ) THEN
           CALL Info('SpectralRadiosity','Skipping interval '//I2S(k),Level=12)
           CYCLE
         END IF

         WRITE(Message,'(A,G12.5)') 'Spectral radiosity sources '//I2S(k)//': ',qsum
         CALL Info('SpectralRadiosity',Message,Level=10) 
                    
         ! Initialize matrix equation
         Diag = 0.0_dp
         RHS  = 0.0_dp         
         IF(Newton) RHS_d = 0.0_dp

         IF ( UseFullMatrix ) THEN
            G_full = 0.0_dp
         ELSE
            G % Values = 0.0_dp
         END IF
         
         ! This is the temperature under study for which we will get the emissivities for. 
         Trad = k*dT         
         CALL TabulateSpectralEmissivity(Emissivity,Absorptivity,Trad,.FALSE.,SimpleTdep)
         CALL RadiosityAssembly(RadiationSurfaces,G,Diag)
         DO i=1,RadiationSurfaces
           ! The portion of the emissivity to consider for this radiating element
           Temp = SurfaceTemperature(i)
           q = ( Temp / dT - k )
           IF( ABS(q) < 1 ) THEN
             Black = Sigma*Temp**4
             e = Emissivity(i)
             a = Absorptivity(i)
             r = 1-a
             c = RelAreas(i) / a

             ! As a weight we use linear interpolation.
             ! Perfect hit get weight 1 that goes to zero when hitting next temperature interval. 
             q = 1-ABS(q)
             RHS(i) = -q*c*e*Black
             IF (AccurateNewton) RHS_d(i) = 4*RHS(i)/Temp
           END IF
         END DO
         CALL BlackRadiosityToRHS(RadiationSurfaces,RHS)
         IF (AccurateNewton) CALL BlackRadiosityToRHS(RadiationSurfaces,RHS_d)

         ! This is a checksum since integration over all temperature intervals should go through all the
         ! participating surface elements. 
         totsum = totsum + qsum
         !PRINT *,'Trad:',k,Trad,qsum,totsum
         CALL RadiationLinearSolver(RadiationSurfaces,G,tmpSOL,RHS,Diag,Solver)
         ! Irradiation of this interval absorbed by the surfaces
         CALL AbsorbedIrradiation(RadiationSurfaces,tmpSOL,AG)

         ! Newton linearization including only "self"
         IF( ApproxNewton ) THEN
           AG_d = (4.0_dp/Trad) * AG
         ELSE IF( AccurateNewton ) THEN
           CALL RadiationLinearSolver(RadiationSurfaces,G,tmpSOL_d, &
                        RHS_d,Diag,Solver,Scaling=.FALSE.)
           CALL AbsorbedIrradiation(RadiationSurfaces,tmpSOL_d,AG_d)
         END IF
         
         ! Cumulative absorbed irradiation
         SOL = SOL + AG
         IF( Newton ) SOL_d = SOL_d + AG_d

         EffTemp = EffTemp + Trad * AG
         EffAbs = EffAbs + Emissivity(1:RadiationSurfaces) * AG
       END DO

       ! This should be exactly one!
       WRITE(Message,'(A,G12.5)') 'Checksum for radiosity sources: ',totsum / RadiationSurfaces 
       CALL Info('SpectralRadiosity',Message,Level=5) 

       ! Check for radiation sources:
       RBC = CheckForRadiators(RadiatorPowers,RadiatorTemps)
       IF(RBC) THEN
         Tmin = MINVAL(RadiatorTemps)
         Tmax = MAXVAL(RadiatorTemps)

         IF(ABS(Tmin-Tmax) < 1.0e-6 ) THEN           
           WRITE(Message,'(A,ES12.3)') 'Only radiator temperature: ',Tmin
           CALL Info('SpectralRadiosity',Message,Level=10)
         ELSE           
           WRITE(Message,'(A,ES12.3)') 'Minimum radiator temperature: ',Tmin
           CALL Info('SpectralRadiosity',Message,Level=10)
           WRITE(Message,'(A,ES12.3)') 'Maximum radiator temperature: ',Tmax
           CALL Info('SpectralRadiosity',Message,Level=10)
         END IF

         ALLOCATE( RadiatorSet(SIZE(RadiatorTemps)) )
         RadiatorSet = 0
         k = 0
         DO i=1,SIZE(RadiatorTemps)
           IF(RadiatorSet(i) > 0) CYCLE
           k=k+1
           RadiatorSet(i) = k
           DO j=i+1,SIZE(RadiatorTemps)
             IF(ABS(RadiatorTemps(i)-RadiatorTemps(j)) < 1.0e-6) RadiatorSet(j) = RadiatorSet(i)
           END DO
         END DO
         kmax = k
         CALL Info('SpectralRadiosity','Going through radiators in '//I2S(kmax)//' sets')
                           
         DO k = 1, kmax
           DO j=1,SIZE(RadiatorSet)
             IF(RadiatorSet(j) == k) Trad = RadiatorTemps(j)
           END DO
           
           WRITE(Message,'(A,G12.5)') 'Spectral radiosity radiators '//I2S(k)//' at: ',Trad
           CALL Info('SpectralRadiosity',Message,Level=10) 
           
           ! Initialize matrix equation
           Diag = 0.0_dp
           RHS = 0.0_dp
           IF ( UseFullMatrix ) THEN
             G_full = 0.0_dp
           ELSE
             G % Values = 0.0_dp
           END IF

           CALL TabulateSpectralEmissivity(Emissivity,Absorptivity,Trad,.TRUE.,SimpleTdep)
           CALL RadiosityAssembly(RadiationSurfaces,G,Diag)
           Rdir = 0.0_dp
           DO i=1,RadiationSurfaces
             Element => Mesh % Elements(ElementNumbers(i))
             IF(ALLOCATED(Element % BoundaryInfo % Radiators)) THEN
               DO j=1,SIZE(RadiatorSet)
                 IF(RadiatorSet(j) == k) THEN
                   Rdir(i) = Rdir(i) + Element % BoundaryInfo % Radiators(j) * RadiatorPowers(j)
                 END IF
               END DO
               a = Absorptivity(i) 
               r = 1-a
               c = RelAreas(i) / a
               RHS(i) = RHS(i) - c * r * Rdir(i)
             END IF
           END DO
           CALL BlackRadiosityToRHS(RadiationSurfaces,RHS)

           CALL RadiationLinearSolver(RadiationSurfaces,G,tmpSOL,RHS,Diag,Solver)          
           ! Absorbed irradiation from the surfaces and directly from the radiators
           CALL AbsorbedIrradiation(RadiationSurfaces,tmpSOL,AG)
           AG = AG + Absorptivity(1:RadiationSurfaces) * Rdir
           
           ! Cumulative absorbed irradiation
           SOL = SOL + AG
           EffTemp = EffTemp + Trad * AG
           EffAbs = EffAbs + Emissivity(1:RadiationSurfaces) * AG
         END DO
       END IF       
       
       ! Normalize with weight i.e. incoming heat flux
       EffAbs = EffAbs / SOL
       EffTemp = EffTemp / SOL

       ! Store the results for access by e.g. heat equation solvers:
       !------------------------------------------------------------
       IF(Newton) THEN
         CALL UpdateRadiosityFactors(SOL,SOL_d,EffAbs,EffTemp)
       ELSE
         CALL UpdateRadiosityFactors(SOL,EffAbs=EffAbs,EffTemp=EffTemp)
       END IF
     END SUBROUTINE SpectralRadiosity
     

     ! Check whether external radiation sources present:
     ! -------------------------------------------------
     FUNCTION CheckForRadiators(RadiatorPowers,RadiatorTemps) RESULT(RBC)
       LOGICAL :: RBC
       REAL(KIND=dp), ALLOCATABLE :: RadiatorPowers(:)
       REAL(KIND=dp), ALLOCATABLE, OPTIONAL :: RadiatorTemps(:)

       TYPE(ValueList_t), POINTER :: RadList
       INTEGER :: i,t,n
       LOGICAL :: Found
       REAL(KIND=dp), POINTER :: RadiatorCoords(:,:), rWrk(:,:)

       ! If radiator is in body force section then use it:
       ! This will make it easier to make GUIs etc.
       IF( .NOT. ListCheckPresentAnyBodyForce( Model,'Radiator Coordinates',RadList ) ) &
             RadList => TSolver % Values
       CALL GetConstRealArray( RadList, RadiatorCoords, 'Radiator Coordinates',RBC)

       IF(RBC) THEN
         n = SIZE(RadiatorCoords,1)
         ALLOCATE( RadiatorPowers(n))
         CALL GetConstRealArray( RadList, rWrk, 'Radiator Power', Found ) 
         IF( Found ) THEN
           IF(SIZE(rWrk,1)==1) THEN
             RadiatorPowers(1:n) = rWrk(1,1)
           ELSE IF(SIZE(rWrk,1)==n) THEN
             RadiatorPowers(1:n) = rWrk(1:n,1)
           ELSE
             CALL Fatal('ConstantRadiosity','Mismatch between size of "Radiator Coordinates" and "Radiator Power"')
           END IF
         ELSE
           DO t=1,n
             RadiatorPowers(t) = ListGetCReal(RadList, 'Radiator Power '//I2S(t),UnfoundFatal=.TRUE.)
           END DO
         END IF

         IF(PRESENT(RadiatorTemps)) THEN
           ALLOCATE( RadiatorTemps(n))
           CALL GetConstRealArray( RadList, rWrk, 'Radiator Temperature', Found )
           IF( Found ) THEN
             IF(SIZE(rWrk,1)==1) THEN
               RadiatorTemps(1:n) = rWrk(1,1)
             ELSE IF(SIZE(rWrk,1)==n) THEN
                 RadiatorTemps(1:n) = rWrk(1:n,1)
             ELSE
               CALL Fatal('SpectralRadiosity','Mismatch between size of "Radiator Coordinates" and "Radiator Power"')
             END IF
           ELSE
             DO t=1,n
               RadiatorTemps(t) = ListGetCReal(RadList, 'Radiator Temperature '//I2S(t),UnfoundFatal=.TRUE.)
             END DO
           END IF
         END IF
       END IF
     END FUNCTION CheckForRadiators


     ! Assemble the LHS of the radiosity equation (r*F-1)J = -e*sigma*T^4 - r*R*P.
     ! The unknown is x = (a/r)*J and row i is multiplied by A_i/a_i, so that
     ! the matrix A_i*F_ij*(r_i/a_i)*(r_j/a_j) is symmetric since A_i*F_ij = A_j*F_ji,
     ! and the diagonal is -A_i*r_i/a_i**2. The RHS should be multiplied by A_i/a_i
     ! too. For black surfaces (r=0) the row reduces to -(A_i/a_i)*J_i = RHS_i and
     ! the columns vanish, so for them J itself is the unknown with diagonal -A_i.
     ! Their contribution to the gray rows is moved to the RHS by BlackRadiosityToRHS.
     ! The heat equation needs the absorbed irradiation, see AbsorbedIrradiation.
     ! ---------------------------------------------------------------------------
     SUBROUTINE RadiosityAssembly(n,G,Diag)
       TYPE(Matrix_t) :: G
       INTEGER :: n
       REAL(KIND=dp) :: Diag(:)

       REAL(KIND=dp), POINTER :: Vals(:) 
       INTEGER, POINTER :: Cols(:) 
       INTEGER :: i, j, nf, previ
       REAL(KIND=dp) :: s,r,e,a,rj,ej,aj,c

       DO i=1,n
         nf = ViewFactors(i) % NumberOfFactors
         Vals => ViewFactors(i) % Factors
         Cols => ViewFactors(i) % Elements
!        e = Emissivity(i)
         a = Absorptivity(i)
         r = 1-a !e
         IF( IsBlack(i) ) THEN
           c = RelAreas(i)
         ELSE
           c = RelAreas(i) * r / a**2
         END IF
         IF( .NOT. UseFullMatrix ) previ = G % Rows(i)-1
         DO j=1,nf
!          ej = Emissivity(Cols(j))
           aj = Absorptivity(Cols(j))
           rj = 1-aj !ej
           s = Vals(j) * (r/a) * (rj/aj)
           IF ( UseFullMatrix ) THEN
             G_full(i,Cols(j)) = G_full(i,Cols(j)) + s
           ELSE
             CALL CRS_AddToMatrixElement(G,i,Cols(j),s,previ)
           END IF
         END DO
         Diag(i) = -c
         IF( UseFullMatrix ) THEN
           G_full(i,i) = Diag(i)
         ELSE
           CALL CRS_AddToMatrixElement(G,i,i,Diag(i))
         END IF
       END DO
     END SUBROUTINE RadiosityAssembly


     ! Black surfaces have zero reflectivity and hence a priori known radiosity.
     ! ---------------------------------------------------------------------------
     FUNCTION IsBlack(i) RESULT(Black)
       INTEGER :: i
       LOGICAL :: Black
       Black = ( 1-Absorptivity(i) < 1.0e-10_dp )
     END FUNCTION IsBlack


     ! The radiosity of the black surfaces is J_i = -RHS_i/A_i. Move their
     ! contribution (r_i/a_i)*A_i*F_ij*J_j to the RHS of the gray rows.
     ! ---------------------------------------------------------------------------
     SUBROUTINE BlackRadiosityToRHS(n,RHS)
       INTEGER :: n
       REAL(KIND=dp) :: RHS(:)

       REAL(KIND=dp), POINTER :: Vals(:)
       INTEGER, POINTER :: Cols(:)
       INTEGER :: i, j, k
       REAL(KIND=dp) :: a, s

       DO i=1,n
         IF( IsBlack(i) ) CYCLE
         Vals => ViewFactors(i) % Factors
         Cols => ViewFactors(i) % Elements
         s = 0.0_dp
         DO j=1,ViewFactors(i) % NumberOfFactors
           k = Cols(j)
           IF( IsBlack(k) ) s = s - Vals(j) * RHS(k) / RelAreas(k)
         END DO
         a = Absorptivity(i)
         RHS(i) = RHS(i) - (1-a)/a * s
       END DO
     END SUBROUTINE BlackRadiosityToRHS


     ! Absorbed irradiation a_i*G_i = a_i*sum_j F_ij*J_j, where the radiosity of
     ! the gray surfaces is J = (r/a)*x and of the black ones J = x, x being the
     ! solution of the radiosity system. Computed directly from the view factors
     ! instead of from G = (J-e*sigma*T^4)/r which is singular for black surfaces.
     ! Setting AbsG=.FALSE. gives just the irradiation G.
     ! ---------------------------------------------------------------------------
     SUBROUTINE AbsorbedIrradiation(n,SOL,AG,AbsG)
       INTEGER :: n
       REAL(KIND=dp) :: SOL(:), AG(:)
       LOGICAL, OPTIONAL :: AbsG

       REAL(KIND=dp), POINTER :: Vals(:)
       REAL(KIND=dp), ALLOCATABLE :: J(:)
       INTEGER, POINTER :: Cols(:)
       INTEGER :: i, k
       REAL(KIND=dp) :: a, s
       LOGICAL :: Absorbed

       Absorbed = .TRUE.
       IF(PRESENT(AbsG)) Absorbed = AbsG

       ALLOCATE(J(n))
       DO i=1,n
         IF( IsBlack(i) ) THEN
           J(i) = SOL(i)
         ELSE
           a = Absorptivity(i)
           J(i) = (1-a)/a * SOL(i)
         END IF
       END DO

       DO i=1,n
         Vals => ViewFactors(i) % Factors
         Cols => ViewFactors(i) % Elements
         s = 0.0_dp
         DO k=1,ViewFactors(i) % NumberOfFactors
           s = s + Vals(k) * J(Cols(k))
         END DO
         AG(i) = s / RelAreas(i)
         IF( Absorbed ) AG(i) = Absorptivity(i) * AG(i)
       END DO
     END SUBROUTINE AbsorbedIrradiation


     
     ! Scale & solve given linear system Ax=b:
     !----------------------------------------
     SUBROUTINE RadiationLinearSolver(n, A, x, b, Diag,  Solver, Scaling)

       INTEGER :: n
       REAL(KIND=dp), TARGET :: x(n), b(n), Diag(n)
       TYPE(Matrix_t), POINTER :: A
       LOGICAL, OPTIONAL :: Scaling
       TYPE(Solver_t), POINTER :: Solver

       LOGICAL :: Scal,Found
       REAL(KIND=dp) :: bscal, eps
       INTEGER :: i,j, maxiter, FirstActive

       ! Solve serially and distribute the result afterwards, memory bandwidth
       ! destroys the performance otherwise (at least for non-supercomputer systems)
       ! A serially loaded radiation mesh is the same in all partitions:
       ! then each partition solves the (small) system itself.
       IF ( ParEnv % PEs <= 1 .OR. GeneralMesh ) THEN
         FirstActive = ParEnv % myPE
       ELSE
         FirstActive = -1
         DO i=0,ParEnv % PEs-1
           IF (ActiveTasks(i)) THEN
             FirstActive=i; EXIT
           END IF
         END DO
       END IF

       scal = .TRUE.
       IF(PRESENT(Scaling)) scal = Scaling

       IF ( ParEnv % myPE == FirstActive ) THEN
         ! Scale matrix to unit diagonals (if not done already)
         IF(scal) THEN
           Diag = SQRT(1._dp/ABS(Diag))
           IF (UseFullMatrix) THEN
             DO j=1,n
               DO i=1,n
                 G_full(i,j) = G_full(i,j)*Diag(i)*Diag(j)
               END DO
             END DO
           ELSE
             DO i=1,n
               DO j=A % Rows(i),A % Rows(i+1)-1
                 A % Values(j) = A % Values(j)*Diag(i)*Diag(A % Cols(j))
               END DO
             END DO
           END IF
         END IF
         b = b * Diag

         ! Scale rhs to one!
         bscal = SQRT(SUM(b**2))
         x = 0.0_dp
         ! The RHS may be zero, e.g. when radiators light up only black surfaces
         IF( bscal > TINY(bscal) ) THEN
           b = b / bscal

           IF(IterSolveFactors) THEN
             Solver % Matrix => A

             eps = ListGetCReal( Params,'Linear System Convergence Tolerance', Found)
             IF  (.NOT. Found ) eps = 1.0d-8
             maxiter = ListGetInteger( Params,'Linear System Max Iterations', Found)
             IF  (.NOT. Found ) maxiter = 100

             IF (UseFullMatrix) THEN
               BLOCK
                 TYPE(Matrix_t), POINTER :: Gm
                 fm_G => G_full
                 Gm => AllocateMatrix()
                 Gm % NumberOfRows = n
                 mvProc = ADDRFUNC(fm_MatVec)
                 CALL RadiationCG( n, Gm, x, b, eps, maxiter )
                 DEALLOCATE(Gm)
               END BLOCK
             ELSE
               CALL RadiationCG( n, A, x, b, eps, maxiter )
             END IF
           ELSE
             CALL DirectSolver( A, x, b, Solver )
           END IF
           x = x * bscal * Diag
         END IF
       END IF

       IF ( ParEnv % Pes <= 1 .OR. GeneralMesh ) RETURN

       ! Distribute the linear system result
       BLOCK
         INTEGER :: sz, Status(MPI_STATUS_SIZE), ierr
         INTEGER, ALLOCATABLE :: SendInfo(:), RecvInfo(:), RecvPerm(:)
         REAL(KIND=dp), ALLOCATABLE :: y(:)

         IF(ParEnv % myPE==FirstActive ) THEN
           ALLOCATE(SendInfo(n))
           DO i=1,n
             Element => Mesh % Elements(ElementNumbers(i))
             SendInfo(i) = Element % GElementIndex
           END DO

           DO i=0,ParEnv % PEs-1
             IF (i==ParEnv % myPE .OR. .NOT.ActiveTasks(i) ) CYCLE
             CALL MPI_BSEND(SendInfo,n,MPI_INTEGER,i,12006,ELMER_COMM_WORLD,ierr)
             CALL MPI_BSEND(x,n,MPI_DOUBLE_PRECISION,i,12007,ELMER_COMM_WORLD,ierr)
           END DO
         ELSE
           ALLOCATE(RecvInfo(n), y(n))
           CALL MPI_RECV( RecvInfo,n,MPI_INTEGER,FirstActive,12006,ELMER_COMM_WORLD,status,ierr )
           CALL MPI_RECV( y,n,MPI_DOUBLE_PRECISION,FirstActive,12007,ELMER_COMM_WORLD,status,ierr )

           ! create a permutation for finding the correct place for the result ...
           sz  = 0
           DO i=Mesh % NumberOfBulkElements+1,Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
             sz = MAX(sz, Mesh % Elements(i) % GelementIndex)
           END DO

           ALLOCATE(RecvPerm(sz))
           DO i=1,n
             Element => Mesh % Elements(ElementNumbers(i))
             RecvPerm(Element % GelementIndex) = i
           END DO
           ! ... and store the result ...
           DO i=1,n
             x(RecvPerm(RecvInfo(i))) = y(i)
           END DO
         END IF
       END BLOCK
     END SUBROUTINE RadiationLinearSolver

     ! Tailored local CG algo for speed testing (somewhat faster than any of the 
     ! library routines but not so much...)
     !-------------------------------------------------------------------------
     SUBROUTINE  RadiationCG( n, A, x, b, eps, maxiter )
       REAL(KIND=dp) :: x(n),b(n), eps
       INTEGER :: n, maxiter
       TYPE(Matrix_t), POINTER :: A

       REAL(KIND=dp):: alpha, beta, rho, oldrho
       REAL(KIND=dp), ALLOCATABLE :: r(:), p(:), q(:)
       REAL(KIND=dp) :: s
       INTEGER :: iter, i, j, k
       REAL(KIND=dp) :: residual, eps2,st

       ALLOCATE(r(n), p(n), q(n))

       eps2 = eps*eps

       IF ( UseFullMatrix) THEN
         CALL DGEMV('N',n,n,1.0_dp,G_full,n,x,1,0.0_dp,r,1)
       ELSE
         CALL CRS_MatrixVectorMultiply(A,x,r) 
       END IF
       r = b - r
       residual = SUM(r*r)
       IF(residual<eps2) RETURN

       DO iter=1,maxiter
         rho = SUM(r*r)
         IF(rho==0.0_dp) ERROR STOP 'CG, rho=0'
  
         IF ( iter==1 ) THEN
           p = r
         ELSE
           beta = rho / oldrho
           p = r + beta * p
         END IF

         IF ( UseFullMatrix) THEN
           CALL DGEMV('N',n,n,1.0_dp,G_full,n,p,1,0.0_dp,q,1)
         ELSE
           CALL CRS_MatrixVectorMultiply(A,p,q) 
         END IF
         alpha = rho/SUM(p*q)

         x = x + alpha * p
         r = r - alpha * q
         residual = SUM(r*r)
         IF ( residual < eps2) EXIT

         oldrho = rho
       END DO

       IF ( UseFullMatrix) THEN
         CALL DGEMV('N',n,n,1.0_dp,G_full,n,x,1,0.0_dp,r,1)
       ELSE
         CALL CRS_MatrixVectorMultiply(A,x,r) 
       END IF
       r = b - r
       residual = SQRT(SUM(r*r))
       WRITE (*, '(I8, E11.4)') iter, residual
       IF( residual > eps ) THEN
         WRITE(Message,'(A,I0,A,ES11.4)') 'Radiosity CG not converged in ',maxiter,&
             ' iterations, residual: ',residual
         CALL Warn('RadiationCG',Message)
       END IF
       DEALLOCATE(r, p, q)
     END SUBROUTINE RadiationCG


     ! Update the outside (heat equation solver) view of the radiosities:
     ! ------------------------------------------------------------------
     SUBROUTINE UpdateRadiosityFactors(SOL,SOL_d,EffAbs,EffTemp)
       REAL(KIND=dp) :: SOL(:)
       REAL(KIND=dp), OPTIONAL :: SOL_d(:), EffAbs(:), EffTemp(:)
       
       TYPE(Element_t), POINTER :: Element
       INTEGER :: i
       TYPE(Factors_t), POINTER :: RadiosityFactors
                
       DO i=1,RadiationSurfaces
         Element => Mesh % Elements(ElementNumbers(i))

         RadiosityFactors => Element % BoundaryInfo % RadiationFactors       
         IF ( .NOT. ASSOCIATED( RadiosityFactors ) ) THEN
           ALLOCATE(RadiosityFactors)
           Element % BoundaryInfo % RadiationFactors => RadiosityFactors
         END IF

         IF (.NOT.ALLOCATED(RadiosityFactors % Elements)) THEN
           ALLOCATE( RadiosityFactors % Elements(1) )
           ALLOCATE( RadiosityFactors % Factors(4) )
           RadiosityFactors % Factors = 0.0_dp
           RadiosityFactors % NumberOfFactors = 1
           RadiosityFactors % Elements(1) = ElementNumbers(i)
         END IF

         RadiosityFactors % Factors(1) = SOL(i)
         IF(Newton .AND. PRESENT(SOL_d)) RadiosityFactors % Factors(2) = SOL_d(i)
         IF(PRESENT(EffAbs))  RadiosityFactors % Factors(3) = EffAbs(i)
         IF(PRESENT(EffTemp)) RadiosityFactors % Factors(4) = EffTemp(i)
       END DO

       IF(InfoActive(30)) THEN
         PRINT *,'SOL_0 range:',MINVAL(SOL),MAXVAL(SOL),SUM(SOL)/SIZE(SOL)       
         IF(Newton .AND. PRESENT(SOL_d)) PRINT *,'SOL_d range:',MINVAL(SOL_d),MAXVAL(SOL_d),SUM(SOL_d)/SIZE(SOL_d)
       END IF
       
     END SUBROUTINE UpdateRadiosityFactors



     ! Save factors is mainly for debugging purposes
     !-------------------------------------------------------------------
     SUBROUTINE SaveGebhartFactors()

       CHARACTER(:), ALLOCATABLE :: OutputName, GebhartFactorsFile
       INTEGER ::  i,t,n
       INTEGER, POINTER :: Cols(:)
       REAL(KIND=dp), POINTER :: Vals(:)
       TYPE(Factors_t), POINTER :: GebhartFactors

       GebhartFactorsFile = GetString(Model % Simulation, 'Gebhart Factors',Found )
       IF (.NOT.Found) &
         GebhartFactorsFile = GetString(Model % Simulation, 'Gebhardt Factors',Found )

       IF ( .NOT.Found ) THEN
         GebhartFactorsFile = 'GebhartFactors.dat'
       END IF

       IF ( LEN_TRIM(MeshDirName) > 0 ) THEN
         OutputName = TRIM(OutputPath) // '/' // TRIM(MeshDirName) // &
             '/' // GebhartFactorsFile
       ELSE
         OutputName = GebhartFactorsFile
       END IF

       IF(RadiationBody > 1) OutputName = OutputName//I2S(RadiationBody)

       OPEN( VFUnit,File=OutputName )

       WRITE (Message,'(A,A)') 'Writing Gephardt Factors to file: ',OutputName
       CALL Info(Caller,Message,Level=5)

       WRITE( VFUnit,* ) RadiationSurfaces

       DO t=1,RadiationSurfaces
         WRITE(VFUnit,*) t,ElementNumbers(t)
       END DO

       DO t=1,RadiationSurfaces
         Element => Mesh % Elements(ElementNumbers(t))
         GebhartFactors => Element % BoundaryInfo % RadiationFactors

         n = GebhartFactors % NumberOfFactors 
         Vals => GebhartFactors % Factors
         Cols => GebhartFactors % Elements

         WRITE( VFUnit,* ) n
         DO i=1,n
           WRITE(VFUnit,*) t,InvElementNumbers(Cols(i)-nBulk),Vals(i)
         END DO
       END DO

       CLOSE(VFUnit)
     END SUBROUTINE SaveGebhartFactors


     ! Update parameters for the iterative linear equation solver (default+user
     ! selected).
     ! -------------------------------------------------------------------------
     SUBROUTINE InitRadiationSolver(TSolver, Solver)
       TYPE(Solver_t) :: TSolver
       TYPE(Solver_t), POINTER :: Solver

       IF(.NOT. ASSOCIATED(Solver)) ALLOCATE(Solver)
       IF(.NOT. ASSOCIATED(Solver % Values)) Solver % Values => ListAllocate()

       CALL ListCopyPrefixedKeywords( TSolver % Values, Solver % Values, 'radiation:' )

       CALL ListAddNewString(Solver % Values,'Linear System Iterative Method', 'CGS')
       CALL ListAddNewString(Solver % Values,'Linear System Direct Method','Umfpack')
       CALL ListAddNewInteger(Solver % Values,'Linear System Max Iterations',500)
       CALL ListAddNewInteger(Solver % Values,'Linear System Residual Output',10 )
       CALL ListAddNewString(Solver % Values,'Linear System Preconditioning','None' )
       CALL ListAddNewConstReal(Solver % Values,'Linear System Convergence Tolerance',1.0d-9)
     END SUBROUTINE InitRadiationSolver

   END SUBROUTINE RadiationFactorsMesh
