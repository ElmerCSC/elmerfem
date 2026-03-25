! *****************************************************************************/
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
! *  Original Date: 02 Apr 2001
! *
! *****************************************************************************/
  
!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Utilities used for interface conditions of various types including mortar.
!------------------------------------------------------------------------------

MODULE MeshAllocations

  USE Types
  USE Messages
  USE ElementDescription
!  USE BandwidthOptimize
!  USE Interpolation, ONLY : PointInElement
!  USE ParallelUtils
  USE Lists
!  USe ListMatrix
  USE ElementUtils, ONLY : FreeMatrix !Find_Face, Find_Edge, AllocateMesh, FreeMatrix, TangentDirections !mGetBoundaryIndexesFromParent, &
!        NormalDirection, CreateMatrix, TangentDirections, &
!        FreeMatrix
  IMPLICIT NONE

CONTAINS
  
!> Allocate mesh structure and return handle to it.
!------------------------------------------------------------------------------
   FUNCTION AllocateMesh(NumberOfBulkElements, NumberOfBoundaryElements, &
       NumberOfNodes, InitParallel ) RESULT(Mesh)
!------------------------------------------------------------------------------
     INTEGER, OPTIONAL :: NumberOfBulkElements, NumberOfBoundaryElements, NumberOfNodes
     LOGICAL, OPTIONAL :: InitParallel
     TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
     INTEGER :: istat, i, n
     CHARACTER(*), PARAMETER :: Caller = 'AllocateMesh'
     
     ALLOCATE( Mesh, STAT=istat )
     IF ( istat /= 0 ) CALL Fatal( Caller, 'Unable to allocate a few bytes of memory?' )

!    Nothing computed on this mesh yet!
!    ----------------------------------
     Mesh % SavesDone    = 0
     Mesh % OutputActive = .FALSE.

     Mesh % AdaptiveDepth = 0
     Mesh % Changed   = .FALSE. !  TODO: Change this sometime
     Mesh % Stabilize = .FALSE.
     Mesh % MeshTag = 1

     Mesh % Variables => NULL()
     Mesh % Parent => NULL()
     Mesh % Child => NULL()
     Mesh % Next => NULL()
     Mesh % RootQuadrant => NULL()
     Mesh % Edges => NULL()
     Mesh % Faces => NULL()
     Mesh % Projector => NULL()
     Mesh % NumberOfEdges = 0
     Mesh % NumberOfFaces = 0

     Mesh % NumberOfBulkElements = 0
     Mesh % NumberOfBoundaryElements = 0
     Mesh % Elements => NULL()
     
     Mesh % DiscontMesh = .FALSE.
     Mesh % SingleMesh  = .FALSE.
     Mesh % InvPerm => NULL()

     Mesh % MinFaceDOFs = 1000
     Mesh % MinEdgeDOFs = 1000
     Mesh % MaxNDOFs = 0
     Mesh % MaxFaceDOFs = 0
     Mesh % MaxEdgeDOFs = 0
     Mesh % MaxBDOFs = 0
     Mesh % MaxElementDOFs  = 0
     Mesh % MaxElementNodes = 0

     Mesh % ViewFactors => NULL()

     ALLOCATE( Mesh % Nodes, STAT=istat )
     IF ( istat /= 0 ) CALL Fatal( Caller, 'Unable to allocate a few bytes of memory?' )
     
     NULLIFY( Mesh % Nodes % x )
     NULLIFY( Mesh % Nodes % y )
     NULLIFY( Mesh % Nodes % z )
     Mesh % Nodes % NumberOfNodes = 0
     Mesh % NumberOfNodes = 0
       
     Mesh % NodesOrig => Mesh % Nodes
     NULLIFY( Mesh % NodesMapped )

     Mesh % EntityWeightsComputed = .FALSE.
     Mesh % BCWeight => NULL()
     Mesh % BodyForceWeight => NULL()
     Mesh % BodyWeight => NULL()
     Mesh % MaterialWeight => NULL()
    
     Mesh % ParallelInfo % NumberOfIfDOFs =  0        
     NULLIFY( Mesh % ParallelInfo % GlobalDOFs )
     NULLIFY( Mesh % ParallelInfo % GInterface )
     NULLIFY( Mesh % ParallelInfo % NeighbourList )     

     i = 0
     IF( PRESENT( NumberOfBulkElements ) ) THEN       
       Mesh % NumberOfBulkElements = NumberOfBulkElements
       i = i + 1
     END IF
     
     IF( PRESENT( NumberOfBoundaryElements ) ) THEN
       Mesh % NumberOfBoundaryElements = NumberOfBoundaryElements
       i = i + 1
     END IF

     IF( PRESENT( NumberOfNodes ) ) THEN
       Mesh % NumberOfNodes = NumberOfNodes
       i = i + 1
     END IF
     
     IF( i > 0 ) THEN
       IF( i < 3 ) CALL Fatal(Caller,'Either give all or no optional parameters!')
       CALL InitializeMesh( Mesh, InitParallel )         
     END IF       
     
!------------------------------------------------------------------------------
   END FUNCTION AllocateMesh
!------------------------------------------------------------------------------

   ! Initialize mesh structures after the size information has been 
   ! retrieved.
   !----------------------------------------------------------------
   SUBROUTINE InitializeMesh(Mesh, InitParallel)     
     TYPE(Mesh_t), POINTER :: Mesh
     LOGICAL, OPTIONAL :: InitParallel
     
     INTEGER :: i,j,k,NoElems,istat
     TYPE(Element_t), POINTER :: Element
     CHARACTER(*), PARAMETER :: Caller = 'InitializeMesh'
     LOGICAL :: DoParallel
     
     IF( Mesh % NumberOfNodes == 0 ) THEN
       CALL Warn(Caller,'Mesh has zero nodes!')
       RETURN
     ELSE
       CALL Info(Caller,'Number of nodes in mesh: '&
           //I2S(Mesh % NumberOfNodes),Level=8)
     END IF

     CALL Info(Caller,'Number of bulk elements in mesh: '&
         //I2S(Mesh % NumberOfBulkElements),Level=8)        

     CALL Info(Caller,'Number of boundary elements in mesh: '&
         //I2S(Mesh % NumberOfBoundaryElements),Level=8)        

     Mesh % Nodes % NumberOfNodes = Mesh % NumberOfNodes          

     NoElems = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

     IF( NoElems == 0 ) THEN
       CALL Fatal('InitializeMesh','Mesh has zero elements!')
     END IF

     Mesh % MaxElementDOFs  = 0
     Mesh % MinEdgeDOFs     = 1000
     Mesh % MinFaceDOFs     = 1000
     Mesh % MaxEdgeDOFs     = 0
     Mesh % MaxFaceDOFs     = 0
     Mesh % MaxBDOFs        = 0

     Mesh % DisContMesh = .FALSE.
     Mesh % DisContPerm => NULL()
     Mesh % DisContNodes = 0

     CALL Info(Caller,'Initial number of max element nodes: '&
         //I2S(Mesh % MaxElementNodes),Level=10) 

     ! Allocate the elements
     !-------------------------------------------------------------------------
     CALL AllocateVector( Mesh % Elements, NoElems, Caller )

     DO j=1,NoElems        
       Element => Mesh % Elements(j)        

       Element % DGDOFs = 0
       Element % BodyId = 0
       Element % TYPE => NULL()
       Element % BoundaryInfo => NULL()
       Element % PDefs => NULL()
       Element % DGIndexes => NULL()
       Element % EdgeIndexes => NULL()
       Element % FaceIndexes => NULL()
       Element % BubbleIndexes => NULL()
     END DO

     ! Allocate the nodes
     !-------------------------------------------------------------------------
     CALL AllocateVector( Mesh % Nodes % x, Mesh % NumberOfNodes, Caller )
     CALL AllocateVector( Mesh % Nodes % y, Mesh % NumberOfNodes, Caller )
     CALL AllocateVector( Mesh % Nodes % z, Mesh % NumberOfNodes, Caller )
     
     IF( .NOT. PRESENT( InitParallel ) ) RETURN
     IF( .NOT. InitParallel ) RETURN
     
     CALL Info( Caller,'Allocating parallel info',Level=12)
     
     ALLOCATE(Mesh % ParallelInfo % GlobalDOFs(Mesh % NumberOfNodes), STAT=istat )
     IF ( istat /= 0 ) &
         CALL Fatal( Caller, 'Unable to allocate Mesh % ParallelInfo % NeighbourList' )
     ALLOCATE(Mesh % ParallelInfo % GInterface(Mesh % NumberOfNodes), STAT=istat )
     IF ( istat /= 0 ) &
         CALL Fatal( Caller, 'Unable to allocate Mesh % ParallelInfo % NeighbourList' )
     ALLOCATE(Mesh % ParallelInfo % NeighbourList(Mesh % NumberOfNodes), STAT=istat )
     IF ( istat /= 0 ) &
         CALL Fatal( Caller, 'Unable to allocate Mesh % ParallelInfo % NeighbourList' )
     DO i=1,Mesh % NumberOfNodes
       NULLIFY(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
     END DO
     
   END SUBROUTINE InitializeMesh


   !------------------------------------------------------------------------------
  SUBROUTINE ReleaseMesh( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    TYPE(Projector_t), POINTER :: Projector
    TYPE(Projector_t), POINTER :: Projector1
    TYPE(Variable_t), POINTER  :: Var, Var1
    TYPE(BoundaryInfo_t), POINTER :: bInfo
    INTEGER :: i,j,k
    LOGICAL :: GotIt
    REAL(KIND=dp), POINTER :: ptr(:)
!------------------------------------------------------------------------------
 
!    Deallocate mesh variables:
!    --------------------------

    CALL Info('ReleaseMesh','Releasing mesh variables',Level=15)
    CALL ReleaseVariableList( Mesh % Variables )
    Mesh % Variables => NULL()

!    Deallocate mesh geometry (nodes,elements and edges):
!    ----------------------------------------------------
    IF ( ASSOCIATED( Mesh % Nodes ) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh nodes',Level=15)
      IF ( ASSOCIATED( Mesh % Nodes % x ) ) DEALLOCATE( Mesh % Nodes % x )
      IF ( ASSOCIATED( Mesh % Nodes % y ) ) DEALLOCATE( Mesh % Nodes % y )
      IF ( ASSOCIATED( Mesh % Nodes % z ) ) DEALLOCATE( Mesh % Nodes % z )
      DEALLOCATE( Mesh % Nodes )
    END IF
    Mesh % Nodes => NULL()


    IF ( ASSOCIATED( Mesh % ParallelInfo % GlobalDOFs ) ) &
        DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs )

    IF ( ASSOCIATED( Mesh % ParallelInfo % NeighbourList ) ) THEN 
      DO i=1,Mesh % NumberOfNodes
        IF(ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) &
            DEALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
      END DO
      DEALLOCATE( Mesh % ParallelInfo % NeighbourList )
    END IF

    IF ( ASSOCIATED( Mesh % ParallelInfo % GInterface ) ) &
        DEALLOCATE( Mesh % ParallelInfo % GInterface )

    IF ( ASSOCIATED( Mesh % ParallelInfo % EdgeInterface ) ) &
        DEALLOCATE( Mesh % ParallelInfo % EdgeInterface )

    IF ( ASSOCIATED( Mesh % ParallelInfo % EdgeNeighbourList ) ) THEN 
      DO i=1,Mesh % NumberOfNodes
        IF(ASSOCIATED( Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours ) ) &
            DEALLOCATE( Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours )
      END DO
      DEALLOCATE( Mesh % ParallelInfo % EdgeNeighbourList )
    END IF

    IF ( ASSOCIATED( Mesh % ParallelInfo % FaceInterface ) ) &
        DEALLOCATE( Mesh % ParallelInfo % FaceInterface )

    IF ( ASSOCIATED( Mesh % ParallelInfo % FaceNeighbourList ) ) THEN
      DO i=1,Mesh % NumberOfNodes
        IF(ASSOCIATED( Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours ) ) &
            DEALLOCATE( Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours )
      END DO
      DEALLOCATE( Mesh % ParallelInfo % FaceNeighbourList )
    END IF

    IF ( ASSOCIATED( Mesh % ParallelInfo % EdgeNeighbourList ) ) THEN 
      DO i=1,Mesh % NumberOfNodes
        IF(ASSOCIATED( Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours ) ) &
           DEALLOCATE( Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours )
      END DO
      DEALLOCATE( Mesh % ParallelInfo % EdgeNeighbourList )
    END IF

    IF ( ASSOCIATED( Mesh % Edges ) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh edges',Level=15)
      CALL ReleaseMeshEdgeTables( Mesh )
      Mesh % Edges => NULL()
    END IF

    IF ( ASSOCIATED( Mesh % Faces ) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh faces',Level=15)
      CALL ReleaseMeshFaceTables( Mesh )
      Mesh % Faces => NULL()
    END IF

    IF (ASSOCIATED(Mesh % ViewFactors) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh view factors',Level=15)
      CALL ReleaseMeshFactorTables( Mesh % ViewFactors )
      Mesh % ViewFactors => NULL()
    END IF


!    Deallocate mesh to mesh projector structures:
!    ---------------------------------------------
    Projector => Mesh % Projector
    DO WHILE( ASSOCIATED( Projector ) )
      CALL Info('ReleaseMesh','Releasing mesh projector',Level=15)
      CALL FreeMatrix( Projector % Matrix )
      CALL FreeMatrix( Projector % TMatrix )
      Projector1 => Projector
      Projector => Projector % Next
      DEALLOCATE( Projector1 )
    END DO
    Mesh % Projector => NULL()

    IF(ASSOCIATED(Mesh % InvPerm)) DEALLOCATE(Mesh % InvPerm)

!    Deallocate quadrant tree (used in mesh to mesh interpolation):
!    --------------------------------------------------------------
    IF( ASSOCIATED( Mesh % RootQuadrant ) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh quadrant tree',Level=15)
      CALL FreeQuadrantTree( Mesh % RootQuadrant )
      Mesh % RootQuadrant => NULL()
    END IF

    CALL ReleaseMeshElements( Mesh ) 
         
    Mesh % NumberOfNodes = 0
    Mesh % NumberOfBulkElements = 0
    Mesh % NumberOfBoundaryElements = 0
    
    CALL Info('ReleaseMesh','Releasing mesh finished',Level=15)
    
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMesh
!------------------------------------------------------------------------------


  SUBROUTINE ReleaseMeshElements(Mesh)

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(BoundaryInfo_t), POINTER :: bInfo
    INTEGER :: i, n

    IF(.NOT. ASSOCIATED(Mesh % Elements) ) THEN
      CALL Info('ReleaseMeshElements','Elements not associated, nothing to release',Level=30)
      RETURN
    END IF

    n = SIZE( Mesh % Elements )     
    CALL Info('ReleaseMeshElements','Releasing number of elements: '//I2S(n),Level=30)


    DO i=1,n

      ! Boundaryinfo structure for boundary elements
      !---------------------------------------------
      IF ( Mesh % Elements(i) % Copy ) CYCLE

      IF ( i > Mesh % NumberOfBulkElements ) THEN
        bInfo => Mesh % Elements(i) % BoundaryInfo
        IF ( ASSOCIATED(bInfo) ) THEN
          IF (ASSOCIATED(bInfo % RadiationFactors)) THEN
            IF ( ALLOCATED(bInfo % RadiationFactors % Elements ) ) THEN
              DEALLOCATE(bInfo % RadiationFactors % Elements )
              DEALLOCATE(bInfo % RadiationFactors % Factors )
            END IF
            DEALLOCATE(bInfo % RadiationFactors)
          END IF
          DEALLOCATE(bInfo)
        END IF
      END IF

      IF ( ASSOCIATED( Mesh % Elements(i) % NodeIndexes ) ) &
          DEALLOCATE( Mesh % Elements(i) % NodeIndexes )
      Mesh % Elements(i) % NodeIndexes => NULL()

      IF ( ASSOCIATED( Mesh % Elements(i) % DGIndexes ) ) &
          DEALLOCATE( Mesh % Elements(i) % DGIndexes )
      Mesh % Elements(i) % DGIndexes => NULL()

      IF ( ASSOCIATED( Mesh % Elements(i) % BubbleIndexes ) ) &
          DEALLOCATE( Mesh % Elements(i) % BubbleIndexes )
      Mesh % Elements(i) % BubbleIndexes => NULL()

      ! This creates problems later on!!!
      !IF ( ASSOCIATED( Mesh % Elements(i) % PDefs ) ) &
      !   DEALLOCATE( Mesh % Elements(i) % PDefs )

      Mesh % Elements(i) % PDefs => NULL() 
    END DO

    DEALLOCATE( Mesh % Elements )
    Mesh % Elements => NULL()

    
  END SUBROUTINE ReleaseMeshElements


  
!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshEdgeTables( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: i
    TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Mesh % Edges ) ) THEN
      CALL Info('ReleaseMeshEdgeTables','Releasing number of edges: '&
          //I2S(Mesh % NumberOfEdges),Level=30)
      
       DO i=1,Mesh % NumberOfEdges
          Edge => Mesh % Edges(i)
          IF ( ASSOCIATED( Edge % NodeIndexes ) ) THEN
             DEALLOCATE( Edge % NodeIndexes )
          END IF
          IF ( ASSOCIATED( Edge % BoundaryInfo ) ) THEN
             DEALLOCATE( Edge % BoundaryInfo )
          END IF
       END DO
       DEALLOCATE( Mesh % Edges )

       NULLIFY( Mesh % Edges )
       IF( Mesh % NumberOfEdges == 0 ) RETURN
       Mesh % NumberOfEdges = 0
       
       IF( ASSOCIATED( Mesh % Elements ) ) THEN      
         DO i=1,SIZE(Mesh % Elements)
           IF ( ASSOCIATED( Mesh % Elements(i) % EdgeIndexes ) ) THEN
             DEALLOCATE( Mesh % Elements(i) % EdgeIndexes )
             Mesh % Elements(i) % EdgeIndexes => NULL()
           END IF
         END DO
       END IF
     END IF
       
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshEdgeTables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshFaceTables( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: i
    TYPE(Element_t), POINTER :: Face
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Mesh % Faces ) ) THEN
      CALL Info('ReleaseMeshFaceTables','Releasing number of faces: '&
          //I2S(Mesh % NumberOfFaces))

      DO i=1,Mesh % NumberOfFaces
          Face => Mesh % Faces(i)
          IF ( ASSOCIATED( Face % NodeIndexes ) ) THEN
             DEALLOCATE( Face % NodeIndexes )
          END IF
          IF ( ASSOCIATED( Face % BoundaryInfo ) ) THEN
             DEALLOCATE( Face % BoundaryInfo )
          END IF
       END DO

       DEALLOCATE( Mesh % Faces )
       NULLIFY( Mesh % Faces )
       IF( Mesh % NumberOfFaces == 0 ) RETURN
       
       Mesh % NumberOfFaces = 0

       IF( ASSOCIATED( Mesh % Elements ) ) THEN
         DO i=1,SIZE(Mesh % Elements)
           IF ( ASSOCIATED( Mesh % Elements(i) % FaceIndexes ) ) THEN
             DEALLOCATE( Mesh % Elements(i) % FaceIndexes )
             Mesh % Elements(i) % FaceIndexes => NULL()
           END IF
         END DO
       END IF
     END IF
       
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshFaceTables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshFactorTables( Factors )
!------------------------------------------------------------------------------
    TYPE(Factors_t), POINTER :: Factors(:)
!------------------------------------------------------------------------------
    INTEGER :: i
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Factors ) ) THEN
       DO i=1,SIZE( Factors)
          IF (ALLOCATED(Factors(i) % Factors))  DEALLOCATE(Factors(i) % Factors)
          IF (ALLOCATED(Factors(i) % Elements)) DEALLOCATE(Factors(i) % Elements)
       END DO
       DEALLOCATE(  Factors )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshFactorTables
!------------------------------------------------------------------------------


END MODULE MeshAllocations



  
