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
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 17 Oct 1996
! *
! *****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!-------------------------------------------------------------------------------
!>  Module for reordering variables for bandwidth and/or gaussian elimination
!>  fillin optimization. Also computes node to element connections (which
!>  implies node to node connections, and thus the global matrix structure).
!-------------------------------------------------------------------------------
MODULE BandwidthOptimize

!-------------------------------------------------------------------------------
   USE ElementDescription
!-------------------------------------------------------------------------------

   IMPLICIT NONE

!-------------------------------------------------------------------------------
   TYPE Label_t
     INTEGER :: Value
     TYPE(Label_t), POINTER :: Next
   END TYPE Label_t

   TYPE LabelPointer_t
     TYPE(Label_t), POINTER :: ListHead
   END TYPE LabelPointer_t

   LOGICAL, PRIVATE :: ForceReorder
!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
!> Subroutine for computing the bandwidth of a sparse matrix.
!-------------------------------------------------------------------------------
   FUNCTION ComputeBandwidth( N, List, Reorder, &
               InvInitialReorder ) RESULT(HalfBandWidth)
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t) :: List(:)
     INTEGER :: n
     INTEGER :: HalfBandWidth
     INTEGER, OPTIONAL, TARGET :: Reorder(:), InvInitialReorder(:)
!-------------------------------------------------------------------------------
     INTEGER :: i,j,k,istat
     TYPE(ListMatrixEntry_t), POINTER :: CList
     INTEGER, POINTER :: pReorder(:), pInvInitialReorder(:)
!-------------------------------------------------------------------------------
     HalfBandWidth = 0

     ! Let's try with pointers as gcc 10 might not like PRESENT within OMP pragmas..
     pReorder => NULL()
     pInvInitialReorder => NULL()

     IF( PRESENT(Reorder) ) pReorder => Reorder
     IF( PRESENT(InvInitialReorder) ) pInvInitialReorder => InvInitialReorder
     
     !$OMP PARALLEL DO &
     !$OMP SHARED(List, pReorder, pInvInitialReorder, N) & 
     !$OMP PRIVATE(Clist, j, k) & 
     !$OMP REDUCTION(max:HalfBandWidth) &
     !$OMP DEFAULT(NONE)
     DO i=1,n
       CList => List(i) % Head
       j = i
       IF ( ASSOCIATED(pInvInitialReorder ) ) j = pInvInitialReorder(j)
       DO WHILE( ASSOCIATED( CList ) )
         k = CList % Index
         IF ( ASSOCIATED(pInvInitialReorder) ) k = pInvInitialReorder(k)
         IF ( ASSOCIATED( pReorder ) ) THEN
           HalfBandwidth = MAX( HalfBandWidth, ABS(pReorder(j)-pReorder(k)) )
         ELSE
           HalfBandwidth = MAX( HalfBandWidth, ABS(j-k) )             
         END IF
         Clist => Clist % Next
       END DO
     END DO
     !$OMP END PARALLEL DO
!-------------------------------------------------------------------------------
   END FUNCTION ComputeBandwidth
!-------------------------------------------------------------------------------


#if 0
   SUBROUTINE OrderPermByMortars(Mesh,Perm)
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: Perm(:)
     
     INTEGER :: SlaveTag, MasterTag, DefaultTag, i,j,k,n
     INTEGER, ALLOCATABLE :: NodeTag(:)
     LOGICAL, ALLOCATABLE :: SlaveBC(:), MasterBC(:)
     TYPE(Element_t), POINTER :: Element
     LOGICAL :: Found
     
     n = CurrentModel % NumberOfBCs
     ALLOCATE(SlaveBC(n), MasterBC(n) )
     SlaveBC = .FALSE.; MasterBC = .FALSE.

     DO i=1,CurrentModel % NumberOfBCs
       j = ListGetInteger( Currentmodel % BCs(i) % Values,'Mortar BC', Found )
       IF(Found ) THEN
         SlaveBC(i) = .TRUE.
         MasterBC(j) = .TRUE.
       END IF
     END DO

     IF(.NOT. ANY(SlaveBC)) RETURN

     ! Tags should have values 1,2,3
     SlaveTag = ListGetInteger( CurrentModel % Solver % Values,'Slave Tag',UnfoundFatal=.TRUE.)
     MasterTag = ListGetInteger( CurrentModel % Solver % Values,'Master Tag',UnfoundFatal=.TRUE.)     
     DefaultTag = 6 - SlaveTag - MasterTag

     ALLOCATE( NodeTag( Mesh % NumberOfNodes ) )
     NodeTag = DefaultTag     
     
     DO i=Mesh % NumberOfBulkElements+1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 
       Element => Mesh % Elements(i)
       IF(.NOT. ASSOCIATED(Element % BoundaryInfo) ) CYCLE
       j = Element % BoundaryInfo % Constraint
       
       IF(SlaveBC(j)) NodeTag(Element % NodeIndexes) = SlaveTag
       IF(MasterBC(j)) NodeTag(Element % NodeIndexes) = MasterTag
     END DO
     
     k = 0
     ! Here we go through cases 1,2,3
     DO j=1,3
       DO i=1, Mesh % NumberOfNodes
         IF(Perm(i)==0) CYCLE
         IF(NodeTag(i) == j) THEN
           k = k + 1
           Perm(i) = k
         END IF
       END DO
     END DO
     
   END SUBROUTINE OrderPermByMortars
#endif

   
!-------------------------------------------------------------------------------
!> Subroutine for reordering variables for bandwidth and/or gaussian elimination
!> fillin optimization. Also computes node to element connections (which
!>  implies node to node connections, and thus the global matrix structure).
!-------------------------------------------------------------------------------
   FUNCTION OptimizeBandwidth( ListMatrix, Perm, InvInitialReorder, LocalNodes, &
             Optimize, UseOptimized, Equation ) RESULT( HalfBandWidth )
!-------------------------------------------------------------------------------
use spariterglobals

     INTEGER, DIMENSION(:) :: Perm, InvInitialReorder
     LOGICAL :: Optimize, UseOptimized
     CHARACTER(LEN=*) :: Equation
     TYPE(ListMatrix_t) :: ListMatrix(:)

     INTEGER :: HalfBandWidth, LocalNodes
!-------------------------------------------------------------------------------

     LOGICAL(KIND=1), ALLOCATABLE :: DoneAlready(:)
     INTEGER, ALLOCATABLE :: PermLocal(:),DoneIndex(:)
     LOGICAL :: Newroot, Finished
     INTEGER :: MinDegree,StartNode,MaxLevel
     INTEGER :: Indx,i,j,k,n,k1,k2,HalfBandWidthBefore,HalfBandWidthAfter,istat

     TYPE(Element_t),POINTER :: Element
     TYPE(ListMatrixEntry_t), POINTER :: p
!-------------------------------------------------------------------------------

     CALL Info( 'OptimizeBandwidth', &
               '---------------------------------------------------------', Level=6 )
     CALL Info( 'OptimizeBandwidth', 'Computing matrix structure for: '// TRIM(Equation) , Level=6)

     HalfBandwidth = ComputeBandWidth( LocalNodes, ListMatrix )+1

     CALL Info( 'OptimizeBandwidth','Initial bandwidth for '&
         //TRIM(Equation)//': '//I2S(HalfBandwidth),Level=4)

     IF ( .NOT.Optimize ) THEN
       CALL Info( 'OptimizeBandwidth', &
               '---------------------------------------------------------', Level=6 )
       RETURN
     END IF

!-------------------------------------------------------------------------------
     HalfBandWidthBefore = HalfBandWidth

!-------------------------------------------------------------------------------
!    Search for node to start
!-------------------------------------------------------------------------------
     StartNode = 1
     MinDegree = ListMatrix(StartNode) % Degree
     DO i=1,LocalNodes
       IF ( ListMatrix(i) % Degree < MinDegree ) THEN
         StartNode = i
         MinDegree = ListMatrix(i) % Degree
       END IF
       ListMatrix(i) % Level = 0
     END DO

     ALLOCATE(DoneAlready(LocalNodes),STAT=istat)
     IF( istat /= 0 ) THEN
       CALL Fatal('OptimizeBandwidth','Allocation error for DoneAlready of size: '&
           //I2S(LocalNodes))
     END IF

     MaxLevel = 0
     DoneAlready = .FALSE.
 
     CALL Levelize( StartNode,0 )
 
     NewRoot = .TRUE.
     DO WHILE( NewRoot )
       NewRoot = .FALSE.
       MinDegree = ListMatrix(StartNode) % Degree
       k = StartNode

       DO i=1,LocalNodes
         IF ( ListMatrix(i) % Level == MaxLevel ) THEN
           IF ( ListMatrix(i) % Degree < MinDegree ) THEN
             k = i
             MinDegree = ListMatrix(i) % Degree
           END IF
         END IF
       END DO

       IF ( k /= StartNode ) THEN
         j = MaxLevel
         MaxLevel = 0
         DoneAlready = .FALSE.

         CALL Levelize( k,0 )

         IF ( j > MaxLevel ) THEN
           NewRoot = .TRUE.
           StartNode = j
         END IF
       END IF
     END DO
!-------------------------------------------------------------------------------
     ALLOCATE( PermLocal(SIZE(Perm)), STAT=istat )
     IF( istat /= 0 ) THEN
       CALL Fatal('OptimizeBandwidth','Allocation error for PermLocal: '//&
           I2S(SIZE(Perm)))
     END IF
     PermLocal = 0

     ALLOCATE( DoneIndex(LocalNodes), STAT=istat )
     IF( istat /= 0 ) THEN
       CALL Fatal('OptimizeBandwidth','Allocation error for DoneIndex: '&
           //I2S(LocalNodes))
     END IF
     DoneIndex = 0
!-------------------------------------------------------------------------------
!    This loop really does the thing
!-------------------------------------------------------------------------------
     Indx = 1
     PermLocal(Indx) = StartNode
     DoneIndex(StartNode) = Indx
     Indx = Indx + 1

     DO i=1,LocalNodes
       IF ( PermLocal(i)==0 ) THEN
         DO j=1,LocalNodes
           IF ( DoneIndex(j)==0 ) THEN
             PermLocal(Indx) = j
             DoneIndex(j) = Indx
             Indx = Indx + 1
             EXIT
           END IF
         END DO
       END IF
       CALL Renumber(ListMatrix(PermLocal(i)) % Head)
     END DO
!-------------------------------------------------------------------------------
!    Store it the other way round for FEM, and reverse order for profile
!    optimization
!-------------------------------------------------------------------------------
     DoneIndex = 0
     DO i=1,LocalNodes
       DoneIndex(PermLocal(i)) = LocalNodes-i+1
     END DO

     PermLocal = Perm
     Perm      = 0
     DO i=1,SIZE(Perm)
       k = PermLocal(i)
       IF (k>0) Perm(i) = DoneIndex(k)
     END DO
     DEALLOCATE( DoneIndex )

     HalfBandWidthAfter = ComputeBandwidth( LocalNodes, &
           ListMatrix,Perm,InvInitialReorder )+1

     CALL Info( 'OptimizeBandwidth','Optimized bandwidth for '&
         //TRIM(Equation)//': '//I2S(HalfBandwidthAfter),Level=4)

     HalfBandWidth = HalfBandWidthAfter

     IF ( HalfBandWidthBefore < HalfBandWidth .AND. .NOT. UseOptimized ) THEN
       CALL Info( 'OptimizeBandwidth',&
             'Bandwidth optimization rejected, using original ordering.',Level=4 )
       HalfBandWidth = HalfBandWidthBefore
       Perm = PermLocal
     END IF
     CALL Info( 'OptimizeBandwidth',&
         '---------------------------------------------------------',Level=6 )

     DEALLOCATE( PermLocal,DoneAlready )
!-------------------------------------------------------------------------------

     CONTAINS

!-------------------------------------------------------------------------------
       SUBROUTINE Renumber( Current )
!-------------------------------------------------------------------------------
         TYPE(ListMatrixEntry_t), POINTER :: Current
!-------------------------------------------------------------------------------
         INTEGER :: k
         TYPE(ListMatrixEntry_t), POINTER :: p
!-------------------------------------------------------------------------------
         p => Current
         DO WHILE( ASSOCIATED(p) )
           k = p % Index
           IF ( k <= LocalNodes ) THEN
             IF ( DoneIndex(k) == 0 ) THEN
               PermLocal(Indx) = k
               DoneIndex(k) = Indx
               Indx = Indx + 1
             END IF
           END IF
           p => p % Next
         END DO
!-------------------------------------------------------------------------------
       END SUBROUTINE Renumber
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
       SUBROUTINE Levelize(nin,Levelin)
!-------------------------------------------------------------------------------
         INTEGER :: nin,Levelin
!-------------------------------------------------------------------------------
         INTEGER :: j,k,n, Level
         TYPE(ListMatrixEntry_t), POINTER :: p=>Null()

         TYPE Stack_t
           TYPE(ListMatrixEntry_t), POINTER :: p
         END TYPE Stack_t

         INTEGER :: stackp
         TYPE(Stack_t), POINTER CONTIG :: stack(:), copystack(:)

         INTEGER :: istat
!-------------------------------------------------------------------------------
         n = nin
         Level=Levelin

         ALLOCATE(stack(512),STAT=istat)
         IF( istat /= 0 ) THEN
           CALL Fatal('Levelize','Allocation error for stack of size 512')
         END IF

         stackp = 0

         p=>ListMatrix(n) % Head
         DO WHILE(ASSOCIATED(p))
           IF ( stackp>=SIZE(stack) ) THEN
             ALLOCATE( copystack(stackp*2) )
             DO j=1,stackp
               copystack(j) = stack(j)
             END DO
             DEALLOCATE(stack)
             stack => copystack
           END IF
           stackp = stackp+1
           stack(stackp) % p => p

           ListMatrix(n) % Level = Level
           DoneAlready(n) = .TRUE.
           MaxLevel = MAX( MaxLevel,Level )

           p => ListMatrix(n) % Head

           DO WHILE(.TRUE.)
             IF ( ASSOCIATED(p) ) THEN
               n = p % Index
               IF ( n <= LocalNodes ) THEN
                 IF ( .NOT.DoneAlready(n) ) THEN
                   Level = Level+1; EXIT
                 END IF
               END IF
             ELSE IF ( stackp>=1 ) THEN
               p => stack(stackp) % p
               Level  = Level-1
               Stackp = Stackp-1
             ELSE
               EXIT
             END IF
             p => p % Next
           END DO
         END DO

         DEALLOCATE(Stack)
!-------------------------------------------------------------------------------
       END SUBROUTINE Levelize
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   END FUNCTION OptimizeBandwidth
!-------------------------------------------------------------------------------

END MODULE BandwidthOptimize

!> \{ ElmerLib
