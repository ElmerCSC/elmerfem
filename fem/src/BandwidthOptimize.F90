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
     INTEGER, OPTIONAL :: Reorder(:), InvInitialReorder(:)
!-------------------------------------------------------------------------------
     INTEGER :: i,j,k,istat
     TYPE(ListMatrixEntry_t), POINTER :: CList
!-------------------------------------------------------------------------------
     HalfBandWidth = 0
     DO i=1,n
        CList => List(i) % Head
        j = i
        IF ( PRESENT( InvInitialReorder ) ) j = InvInitialReorder(j)
        DO WHILE( ASSOCIATED( CList ) )
           k = CList % Index
           IF ( PRESENT(InvInitialReorder) ) k = InvInitialReorder(k)
           IF ( .NOT. PRESENT( Reorder ) ) THEN
              HalfBandwidth = MAX( HalfBandWidth, ABS(j-k) )
           ELSE
              HalfBandwidth = MAX( HalfBandWidth, ABS(Reorder(j)-Reorder(k)) )
           END IF
           Clist => Clist % Next
        END DO
     END DO
!-------------------------------------------------------------------------------
   END FUNCTION ComputeBandwidth
!-------------------------------------------------------------------------------


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
               '---------------------------------------------------------', Level=4 )
     CALL Info( 'OptimizeBandwidth', 'Computing matrix structure for: ' &
                 // TRIM(Equation) //  '...', .TRUE., Level=4)

     HalfBandwidth = ComputeBandWidth( LocalNodes, ListMatrix )+1

     CALL Info( 'OptimizeBandwidth', 'done.', Level=4 )
     WRITE( Message,'(A,I0)' ) 'Half bandwidth without optimization: ', HalfBandwidth
     CALL Info( 'OptimizeBandwidth', Message, Level=4 )

     IF ( .NOT.Optimize ) THEN
       CALL Info( 'OptimizeBandwidth', &
               '---------------------------------------------------------', Level=4 )
       RETURN
     END IF

!-------------------------------------------------------------------------------
     HalfBandWidthBefore = HalfBandWidth

     CALL Info( 'OptimizeBandwidth', ' ', Level=4 )
     CALL Info( 'OptimizeBandwidth', 'Bandwidth Optimization ...', .TRUE.,Level=4 )
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
           //TRIM(I2S(LocalNodes)))
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
           TRIM(I2S(SIZE(Perm))))
     END IF
     PermLocal = 0

     ALLOCATE( DoneIndex(LocalNodes), STAT=istat )
     IF( istat /= 0 ) THEN
       CALL Fatal('OptimizeBandwidth','Allocation error for DoneIndex: '&
           //TRIM(I2S(LocalNodes)))
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
     CALL Info( 'OptimizeBandwidth', 'done.', Level=4 )

     WRITE( Message,'(A,I0)') 'Half bandwidth after optimization: ', HalfBandwidthAfter
     CALL Info( 'OptimizeBandwidth', Message, Level=4 )
     HalfBandWidth = HalfBandWidthAfter

     IF ( HalfBandWidthBefore < HalfBandWidth .AND. .NOT. UseOptimized ) THEN
       CALL Info( 'OptimizeBandwidth',&
             'Bandwidth optimization rejected, using original ordering.',Level=4 )
       HalfBandWidth = HalfBandWidthBefore
       Perm = PermLocal
     END IF
     CALL Info( 'OptimizeBandwidth', &
             '---------------------------------------------------------',Level=4 )

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

!> \deprecated The LexiographicSearch for bandwidth optimization is not in use. Can it be removed? 
#if 0
NOT CURRENT AT THE MOMENT...
   FUNCTION LexiographicSearch(Model,Perm) RESULT(HalfBandWidth)
     TYPE(Model_t), POINTER :: Model
     INTEGER :: Perm(:)

     INTEGER :: HalfBandWidth

     TYPE(LabelPointer_t), POINTER :: Label(:)
     TYPE(Label_t), POINTER :: Lptr1,Lptr2, CurrentLabel

     TYPE(Element_t), POINTER :: Element,Element1
     TYPE(ElementList_t), POINTER :: p,q
     TYPE(ElementListPointer_t), POINTER :: List(:)

     LOGICAL(KIND=1), ALLOCATABLE :: DoneAlready(:)

     INTEGER :: i,j,k,l,m,n,k1,k2,CurrentVertex

     PRINT*,' '
     PRINT*,'Lexiographic search for fillin optimization: '
     PRINT*,'---------------------------------------------------------'

     HalfBandWidth = 0
     DO i=1,Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
       Element => Model % Elements(i)
       DO j=1,Element % Type % NumberOfNodes
         k1 = Element % NodeIndexes(j)
         DO k=j+1,Element % Type % NumberOfNodes
           k2 = Element % NodeIndexes(k)
           HalfBandWidth = MAX(HalfBandWidth,ABS(k1-k2+1))
         END DO
       END DO
     END DO

     PRINT*,'Half bandwidth before optimizing: ',HalfBandWidth

!    CALL NodeToElementList(Model,List)

     n = Model % NumberOfNodes

     ALLOCATE( Label(n) )
     DO i=1,N
       NULLIFY( Label(i) % ListHead )
     END DO


     ALLOCATE( DoneAlready(Model % NumberOfNodes) )
     DoneAlready = .FALSE.
 
     Perm = 0

     DO i=N,1,-1
       NULLIFY( CurrentLabel )
       CurrentVertex = 1

       DO j=1,N
         IF ( .NOT.DoneAlready(j) ) THEN
           IF ( LabelCompare(Label(j) % ListHead, CurrentLabel) > 0 ) THEN
             CurrentLabel => Label(j) % ListHead
             CurrentVertex = j 
           END IF
         END IF
       END DO

       Perm(CurrentVertex) = i
       DoneAlready(CurrentVertex) = .TRUE.

       Lptr1 => Label(CurrentVertex) % ListHead
       DO WHILE( ASSOCIATED(Lptr1) )
         Lptr2 => Lptr1
         DEALLOCATE(Lptr1)
         Lptr1 => Lptr2 % Next
       END DO
      
       p => List(CurrentVertex) % ListHead
       DO WHILE( ASSOCIATED(p) )
         Element => p % Element

         DO j=1,Element % Type % NumberOfNodes
           k = Element % NodeIndexes(j)
           IF ( .NOT.DoneAlready(k) ) THEN
             CALL LabelAdd( Label(k), i )
           END IF
           CALL AddPath( k,i,0 )
         END DO

         p => p % Next
       END DO
     END DO

!
!  WARNING:: Node to element list deallocated here
!
     CALL FreeNodeToElementList( Model % NumberOfNodes,List )
!
!
!
     HalfBandWidth = 0
     DO i=1,Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
       Element => Model % Elements(i)
       DO j=1,Element % Type % NumberOfNodes
         k1 = Perm(Element % NodeIndexes(j))
         DO k=j+1,Element % Type % NumberOfNodes
           k2 = Perm(Element % NodeIndexes(k))
           HalfBandWidth = MAX(HalfBandWidth,ABS(k1-k2)+1)
         END DO
       END DO
     END DO

     PRINT*,'Half bandwidth is:                ',HalfBandWidth
     PRINT*,'---------------------------------------------------------'

     DEALLOCATE( DoneAllready,Label )

     CONTAINS

       RECURSIVE SUBROUTINE AddPath(n,i,Level)
         INTEGER :: n,i,Level

         INTEGER :: j,k
         TYPE(ElementList_t), POINTER :: p

         p => List(n) % ListHead
         DO WHILE( ASSOCIATED(p) ) 
           DO j=1,Element % Type % NumberOfNodes
             k = Element % NodeIndexes(j)

             IF ( .NOT.DoneAlready(k) ) THEN
               IF ( LabelCompare( Label(k) % ListHead, Label(n) % ListHead ) > 0 ) THEN
                 CALL AddPath( k,i,Level+1 )
               ELSE IF ( Level /= 0 )  THEN
                 CALL LabelAdd( Label(n),i )
               END IF
             END IF

           END DO
           p => p % Next
         END DO
       END SUBROUTINE AddPath


       SUBROUTINE LabelAdd(Label,Value)
         TYPE(LabelPointer_t) :: Label
         INTEGER :: Value

         TYPE(Label_t), POINTER :: tmp,prev,p

         NULLIFY(prev)
         p => Label % ListHead
 
         NULLIFY(prev)
         DO WHILE( ASSOCIATED(p) )
           IF ( p % Value  < Value ) EXIT
           IF ( p % Value == Value ) RETURN
           prev => p
           p => p % Next
         END DO
 
         ALLOCATE( tmp )
         tmp % Value = Value
 
         NULLIFY( tmp % Next )
         IF ( ASSOCIATED(p) ) tmp % Next => p
 
         IF ( ASSOCIATED(prev) ) THEN
           prev % Next => tmp
         ELSE
           Label % ListHead => tmp
         END IF
       END SUBROUTINE LabelAdd


       FUNCTION LabelCompare(Label1,Label2) RESULT(res)
         TYPE(Label_t), POINTER :: Label1,Label2
         INTEGER :: Res

         TYPE(Label_t), POINTER :: p,q

         p => Label1
         q => Label2

         DO WHILE( ASSOCIATED(p) .AND. ASSOCIATED(q) )
           IF ( p % Value == q % Value ) THEN
             p => p % Next
             q => q % Next
           ELSE
             IF ( p % Value > q % Value ) THEN
               res =  1
             ELSE
               res = -1
             END IF
             RETURN
           END IF
         END DO

         IF ( .NOT.ASSOCIATED(q) .AND. .NOT.ASSOCIATED(p) ) THEN
           res =  0
         ELSE IF ( ASSOCIATED(p) ) THEN
           res =  1 
         ELSE IF ( ASSOCIATED(q) ) THEN
           res = -1
         END IF
       END FUNCTION LabelCompare

   END FUNCTION LexiographicSearch



   FUNCTION DepthFirstSearch(Model,Perm) RESULT(HalfBandWidth)
     TYPE(Model_t), POINTER :: Model
     INTEGER, DIMENSION(:) :: Perm

     INTEGER :: HalfBandWidth

     LOGICAL(KIND=1), ALLOCATABLE :: DoneAlready(:)
     INTEGER :: MaxDegree,StartNode
     INTEGER :: Indx,i,j,k,n,k1,k2

     TYPE(Element_t),POINTER :: Element
     TYPE(ElementListPointer_t), DIMENSION(:), POINTER :: List

     N = Model % NumberOfNodes

!    CALL NodeToElementList(Model,List)

     MaxDegree = List(1) % Degree
     StartNode = 1
     DO i=1,N
       IF ( List(i) % Degree >= MaxDegree ) THEN
         StartNode = i
         MaxDegree = List(i) % Degree
       END IF
     END DO

     ALLOCATE( DoneAlready(Model % NumberOfNodes) )

     DoneAlready = .FALSE.
     DoneAlready(StartNode) = .TRUE.

     Perm = 0
     Perm(StartNode) = N

     Indx = N-1
     CALL Renumber(List(StartNode) % ListHead)

!
!  WARNING:: Node to element list deallocated here
!
     CALL FreeNodeToElementList( Model % NumberOfNodes,List )
!
!
!
     HalfBandWidth = 0
     DO i=1,Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
       Element => Model % Elements(i)
       DO j=1,Element % Type % NumberOfNodes
         k1 = Perm(Element % NodeIndexes(j))
         DO k=j+1,Element % Type % NumberOfNodes
           k2 = Perm(Element % NodeIndexes(k))
           HalfBandWidth = MAX(HalfBandWidth,ABS(k1-k2))
         END DO
       END DO
     END DO

     DEALLOCATE( DoneAlready )

     CONTAINS

       RECURSIVE SUBROUTINE Renumber(Current)
         TYPE(ElementList_t), POINTER :: Current

         TYPE(ElementList_t), POINTER :: p
         TYPE(Element_t),POINTER :: Element
         INTEGER :: i,j,k,l,n,Ind(500),Deg(500)

         p => Current
         n = 0
         DO WHILE( ASSOCIATED(p) )
           Element => p % Element
           DO i=1,Element % Type % NumberOfNodes
             k = Element % NodeIndexes(i)
             IF ( .NOT.DoneAlready(k) ) THEN
               n = n + 1
               Ind(n) = k
               Deg(n) = List(k) % Degree
             END IF
           END DO
           p => p % Next
         END DO

         CALL SortI(n,Deg,Ind) 

         DO i=n,1,-1
           k = Ind(i)
           IF ( .NOT.DoneAlready(k) ) THEN
             DO j=i-1,1,-1
               IF ( .NOT.DoneAlready(Ind(j)) ) THEN
                 IF ( List(k) % Degree == List(Ind(j)) % Degree ) THEN
                   l = NumberVisited(List(k) % ListHead)
                   IF ( NumberVisited(List(Ind(j)) % ListHead) < l ) THEN
                     l = Ind(j)
                     Ind(j) = k
                     k = l
                   END IF
                 ELSE
                   EXIT
                 END IF
               END IF
             END DO

             Perm(k) = Indx
             Indx = Indx - 1
             DoneAlready(k) = .TRUE.
             CALL Renumber( List(k) % ListHead )
           END IF
         END DO

       END SUBROUTINE Renumber


       FUNCTION NumberVisited(Current) RESULT(n)

         TYPE(ElementList_t),POINTER :: Current
         INTEGER :: n

         TYPE(ElementList_t),POINTER :: p

          n = 0
          p => Current
          DO WHILE( ASSOCIATED(p) )
            Element => p % Element
            DO i=1,Element % Type % NumberOfNodes
              k = Element % NodeIndexes(i)
              IF ( .NOT.DoneAlready(k) ) n = n + 1
            END DO
            p => p % Next
          END DO

       END FUNCTION NumberVisited
 
   END FUNCTION DepthFirstSearch
#endif

END MODULE BandwidthOptimize

!> \{ ElmerLib
