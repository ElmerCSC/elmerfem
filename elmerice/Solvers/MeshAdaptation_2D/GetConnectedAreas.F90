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
! *  Original Date: 01/02/2018
! *  
! *  Find conneceted areas in a mesh
! *****************************************************************************
!!!  
      SUBROUTINE GetConnectedAreas( Model,Solver,dt,TransientSimulation)
      USE DefUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t), TARGET :: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      TYPE Queue_t
         INTEGER :: maxsize
         INTEGER :: top
         INTEGER,ALLOCATABLE :: items(:)
      END TYPE Queue_t

      TYPE(Variable_t), POINTER :: Var,Area
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Queue_t) :: Queue
      TYPE(Element_t),POINTER :: Element,Parent
      TYPE(Element_t),POINTER :: Faces(:),Face
      TYPE(GaussIntegrationPoints_t) :: IP
      TYPE(Nodes_t),SAVE :: Nodes
      TYPE(ValueList_t), POINTER :: SolverParams
      REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), &
                                      ddBasisddx(:,:,:)
      REAL(KIND=dp) :: s, detJ
      REAL(KIND=dp) :: smallest,largest
      REAL(KIND=dp), ALLOCATABLE :: RegionArea(:)
      INTEGER :: EIndex
      INTEGER :: label
      INTEGER :: region
      INTEGER,ALLOCATABLE :: ElementLabel(:)
      INTEGER :: t,i,n,p
      INTEGER :: nfaces
      INTEGER :: NOFActive

      LOGICAL :: stat,Found
      LOGICAL :: SAVE_REGIONS

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='GetConnectedAreas'

      Mesh => Solver % Mesh
      SolverParams => GetSolverParams()

      n = MAX(Mesh % MaxElementNodes,Mesh % MaxElementDOFs)
      ALLOCATE(ElementLabel(Solver%Mesh%NumberOfBulkElements))
      ALLOCATE( Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3))
      ALLOCATE(Nodes%x(n),Nodes%y(n),Nodes%z(n))


      CALL FindMeshEdges(Mesh,.FALSE.)
      SELECT CASE(Mesh % MeshDim)
       CASE(2)
        Faces => Mesh % Edges
       CASE(3)
        Faces => Mesh % Faces
      END SELECT

      CALL QueueInit(Queue,Mesh%NumberOfBulkElements)

      NOFActive=GetNOFActive()

      ElementLabel=-1
      label=1
      
      DO t=1,NOFActive
         Element => GetActiveElement(t)
         EIndex= Element % ElementIndex 

         IF (CheckPassiveElement(Element)) CYCLE
         IF (ElementLabel(EIndex).GT.0) CYCLE

         ElementLabel(EIndex) = label
         CALL QueuePush(Queue,EIndex)

         DO WHILE (Queue%top.GT.0)
            CALL QueuePop(Queue,EIndex)
            Element => Mesh % Elements(EIndex)

            SELECT CASE(Mesh % MeshDim)
             CASE(2)
              nfaces=Element % TYPE % NumberOfEdges
             CASE(3)
              nfaces=Element % TYPE % NumberOfFaces
            END SELECT

            DO i=1,nfaces
              SELECT CASE(Mesh % MeshDim)
               CASE(2)
                Face =>  Mesh % Edges (Element % EdgeIndexes(i))
               CASE(3)
                Face =>  Mesh % Faces (Element % FaceIndexes(i))
              END SELECT

              IF (.NOT.ASSOCIATED(Face)) &
                CALL FATAL(SolverName,'Face not found')

              Parent => Face % BoundaryInfo % Left
              IF (ASSOCIATED(Parent)) THEN 
                 IF ((.NOT.CheckPassiveElement(Parent))&
                      .AND.(ElementLabel(Parent%ElementIndex).LE.0)) THEN
                   ElementLabel(Parent%ElementIndex)=label      
                   CALL QueuePush(Queue,Parent%ElementIndex)
                 END IF
              END IF
              Parent => Face % BoundaryInfo % Right
              IF (ASSOCIATED(Parent)) THEN 
                 IF ((.NOT.CheckPassiveElement(Parent))&
                      .AND.(ElementLabel(Parent%ElementIndex).LE.0)) THEN
                   ElementLabel(Parent%ElementIndex)=label      
                   CALL QueuePush(Queue,Parent%ElementIndex)
                 END IF
              END IF

            END DO           
            
         END DO

         label = label + 1
      END DO
 
      label = label - 1
      WRITE(Message,'(A,i0,A)') 'There is ',label,' unconnected areas'
      CALL INFO(SolverName,TRIM(Message),level=3)

      Var => VariableGet( Mesh % Variables,'AreaNumber',UnfoundFatal=.True.)
      Area => VariableGet( Mesh % Variables,'RegionArea',UnfoundFatal=.True.)

      ALLOCATE(RegionArea(label))
      
      RegionArea=0._dp
      Var % Values = -1
      Area % Values = 0.

      DO t=1,NOFActive
         Element => GetActiveElement(t)
         EIndex= Element % ElementIndex 
         n  = GetElementNOFNodes(Element)

         region=ElementLabel(EIndex)
         IF (region.LE.0) CYCLE

         Var % Values(Var % Perm(Element % NodeIndexes(1:n)))= &
           region

         Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes(1:n))
         Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes(1:n))
         Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes(1:n))
         IP = GaussPoints( Element )
         DO p = 1, IP % n
           stat = ElementInfo( Element, Nodes, IP % U(p), IP % V(p), &
             IP % W(p), detJ, Basis, dBasisdx, ddBasisddx, .FALSE.) 
           s = detJ * IP % S(p)                           

            RegionArea(region) = RegionArea(region)  + s
        END DO


      END DO

      !$omp parallel do &
      !$omp& private(region)
      DO t=1,Mesh % NumberOfNodes
         region=NINT(Var % Values(Var % Perm(t)))
         IF (region.GT.0) THEN
           Area % Values(Area % Perm(t))=& 
             MAX(RegionArea(region),Area % Values(Area % Perm(t)))
         END IF
      END DO
      !$omp end parallel do

      largest=MAXVAL(RegionArea)
      smallest=MINVAL(RegionArea)
      WRITE(Message,'(A,e15.7)') 'largest area :',largest
      CALL INFO(SolverName,TRIM(Message),level=3)
      WRITE(Message,'(A,e15.7)') 'smallest area :',smallest
      CALL INFO(SolverName,TRIM(Message),level=3)

      SAVE_REGIONS=ListGetLogical(SolverParams,'Save regions labels',Found)
      IF (.NOT.Found) SAVE_REGIONS=.FALSE.
      IF (SAVE_REGIONS) THEN
        Open(12,file='regions.txt')
        DO i=1,label
          write(12,*) i,RegionArea(i)
        END DO
        close(12)
      END IF

      DEALLOCATE(ElementLabel)
      DEALLOCATE(RegionArea)
      DEALLOCATE(Basis,dBasisdx,ddBasisddx)
      DEALLOCATE(Queue%items)
      DEALLOCATE(Nodes%x,Nodes%y,Nodes%z)


      CONTAINS 
        SUBROUTINE QueueInit(Queue,n)
        IMPLICIT NONE
        TYPE(Queue_t) :: Queue
        INTEGER :: n
          ALLOCATE(Queue%items(n))
          Queue%top=0
          Queue%maxsize=n
        END SUBROUTINE QueueInit

        SUBROUTINE QueuePush(Queue,x)
        IMPLICIT NONE
        TYPE(Queue_t) :: Queue
        INTEGER :: x
          IF (Queue%top.EQ.Queue%maxsize) &
            CALL FATAL(SolverName,'Too many elements in the queue?')

          Queue%top=Queue%top+1
          Queue%items(Queue%top)=x

        END SUBROUTINE QueuePush

        SUBROUTINE QueuePop(Queue,x)
        IMPLICIT NONE
        TYPE(Queue_t) :: Queue
        INTEGER :: x
          IF (Queue%top.EQ.0) &
            CALL FATAL(SolverName,'No more element in the queue')

          x=Queue%items(Queue%top)
          Queue%top=Queue%top-1
        END SUBROUTINE QueuePop
      END



