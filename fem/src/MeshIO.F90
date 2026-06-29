!/******************************************************************************
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
! *  Authors: Juha Ruokolainen, Peter Raback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Mesh I/O utilities: element allocation helpers, mesh file reading and writing.
!>  Extracted from MeshUtils.
!------------------------------------------------------------------------------

MODULE MeshIO

    USE ElementDescription
    USE ParallelUtils
    USE Lists
    USE ListMatrix
    USE MeshAllocations
    IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
!> Allocated one single element. 
!------------------------------------------------------------------------------
   FUNCTION AllocateElement() RESULT( Element )
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

     ALLOCATE( Element, STAT=istat )
     IF ( istat /= 0 ) &
        CALL Fatal( 'AllocateElement', 'Unable to allocate a few bytes of memory?' )
     Element % BDOFs    =  0
     Element % NDOFs    =  0
     Element % BodyId   = -1
     Element % Splitted =  0
     Element % hK = 0
     Element % ElementIndex = 0
     Element % StabilizationMk = 0
     NULLIFY( Element % TYPE )
     NULLIFY( Element % PDefs )
     NULLIFY( Element % BubbleIndexes )
     NULLIFY( Element % DGIndexes )
     NULLIFY( Element % NodeIndexes )
     NULLIFY( Element % EdgeIndexes )
     NULLIFY( Element % FaceIndexes )
     NULLIFY( Element % BoundaryInfo )
!------------------------------------------------------------------------------
   END FUNCTION AllocateElement
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
   SUBROUTINE AllocatePDefinitions(Element)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: istat,n

     LOGICAL :: Found
     TYPE(Element_t) :: Element

     ! Sanity check to avoid memory leaks
     IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
        ALLOCATE(Element % PDefs, STAT=istat)
        IF ( istat /= 0) CALL Fatal('AllocatePDefinitions','Unable to allocate memory')
     ELSE
       CALL Info('AllocatePDefinitions','P element definitions already allocated',Level=32)
     END IF

     ! Initialize fields
     Element % PDefs % P = 0 
     Element % PDefs % TetraType = 0
     Element % PDefs % isEdge = .FALSE.
     Element % PDefs % localNumber = 0
     Element % PDefs % GaussPoints = 0

     Element % PDefs % Serendipity = ListGetLogical( CurrentModel % Simulation, &
           'Serendipity P elements', Found )
     IF(.NOT.Found) Element % PDefs % Serendipity = .TRUE.
!------------------------------------------------------------------------------
   END SUBROUTINE AllocatePDefinitions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE AllocateBoundaryInfo(Element)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: istat,n

     TYPE(Element_t) :: Element

     ALLOCATE(Element % BoundaryInfo, STAT=istat)
     IF ( istat /= 0) CALL Fatal('AllocateBoundaryInfo','Unable to allocate memory')

     Element % BoundaryInfo % Left => NULL()
     Element % BoundaryInfo % Right => NULL()
     Element % BoundaryInfo % Constraint =  0
     Element % BoundaryInfo % RadiationFactors => NULL()

!------------------------------------------------------------------------------
   END SUBROUTINE AllocateBoundaryInfo
!------------------------------------------------------------------------------


 !> Fortran reader for Elmer ascii and binary mesh file format.
 !> The ascii format is tried out first, if not success, binary is followed. 
 !> This is a Fortran replacement for the old C++ eio library. 
 !------------------------------------------------------------------------
 SUBROUTINE ElmerMeshReader(Step, PMesh, MeshNamePar, ThisPe, NumPEs, &
                 IsParallel, BoundariesOnly )

   IMPLICIT NONE

   INTEGER :: Step
   CHARACTER(LEN=*), OPTIONAL :: MeshNamePar
   TYPE(Mesh_t), POINTER, OPTIONAL :: PMesh
   INTEGER, OPTIONAL :: ThisPe, NumPEs
   LOGICAL, OPTIONAL :: IsParallel
   LOGICAL, OPTIONAL :: BoundariesOnly

   TYPE(Mesh_t), POINTER :: Mesh
   INTEGER :: PrevStep=0, iostat
   INTEGER, PARAMETER :: FileUnit = 10
   INTEGER :: i,j,k,n,BaseNameLen, SharedNodes = 0, mype = 0, numprocs = 0
   INTEGER, POINTER :: NodeTags(:), ElementTags(:), LocalPerm(:)
   INTEGER :: MinNodeTag = 0, MaxNodeTag = 0, istat
   LOGICAL :: ElementPermutation=.FALSE., NodePermutation=.FALSE., Parallel, &
       PseudoParallel, Found
   CHARACTER(:), ALLOCATABLE :: BaseName, FileName


   SAVE PrevStep, BaseName, BaseNameLen, Mesh, mype, Parallel, &
       NodeTags, ElementTags, LocalPerm, PseudoParallel

   CALL Info('ElmerMeshReader','Performing step: '//I2S(Step),Level=8)

   IF( Step - PrevStep /= 1 ) THEN
     CALL Fatal('ElmerMeshReader','The routine should be called in sequence: '// &
         I2S(PrevStep)//' : '//I2S(Step) )
   END IF
   PrevStep = Step
   IF( PrevStep == 6 ) PrevStep = 0 

   IF( Step == 1 ) THEN
     IF(.NOT. PRESENT( MeshNamePar ) ) THEN
       CALL Fatal('ElmerMeshReader','When calling in mode one give MeshNamePar!')
     END IF
     BaseName = TRIM( MeshNamePar ) 
     IF(.NOT. PRESENT( PMesh ) ) THEN
       CALL Fatal('ElmerMeshReader','When calling in mode one give PMesh!')
     END IF
     Mesh => PMesh
     IF(.NOT. PRESENT( ThisPe ) ) THEN
       CALL Fatal('ElmerMeshReader','When calling in mode one give ThisPe!')
     END IF
     mype = ThisPe 
     IF(.NOT. PRESENT( NumPEs) ) THEN
       CALL Fatal('ElmerMeshReader','When calling in mode one give NumPEs!')
     END IF
     numprocs = NumPEs
     IF(.NOT. PRESENT( IsParallel ) ) THEN
       CALL Fatal('ElmerMeshReader','When calling in mode one give IsParallel!')
     END IF
     Parallel = IsParallel

     PseudoParallel = .FALSE.
     IF(.NOT. Parallel ) THEN
       IF( ParEnv % PEs > 1 ) THEN
         PseudoParallel = ListGetLogical(CurrentModel % Simulation,'Enforce Parallel',Found ) 
         IF(.NOT. Found ) PseudoParallel = ListGetLogicalAnySolver(CurrentModel,'Enforce Parallel')
       END IF
     END IF
     
     i = LEN_TRIM(MeshNamePar)
     DO WHILE(MeshNamePar(i:i) == CHAR(0))
       i=i-1
     END DO
     BaseNameLen = i
     CALL Info('ElmerMeshReader','Base mesh name: '//TRIM(MeshNamePar(1:BaseNameLen)))
   END IF
   

   SELECT CASE( Step ) 

   CASE(1)       
     CALL ReadHeaderFile()

   CASE(2)
     CALL ReadNodesFile()

   CASE(3)
     CALL ReadElementsFile()

   CASE(4)
     CALL ReadBoundaryFile(BoundariesOnly)
     CALL PermuteNodeNumbering()

   CASE(5)
     IF( PseudoParallel ) THEN
       CALL InitPseudoParallel()
     ELSE
       CALL InitParallelInfo()
       CALL ReadSharedFile()
     END IF
       
   CASE(6)
     IF( ASSOCIATED( LocalPerm) ) DEALLOCATE( LocalPerm ) 
     IF( ASSOCIATED( ElementTags) ) DEALLOCATE( ElementTags )

   END SELECT


 CONTAINS


   FUNCTION read_ints(s,j,halo) RESULT(n)
     INTEGER :: j(:)
     CHARACTER(LEN=*) :: s
     LOGICAL :: halo
     
     INTEGER :: i,k,l,m,n,ic
     INTEGER, PARAMETER :: ic0 = ICHAR('0'), ic9 = ICHAR('9'), icm = ICHAR('-'), &
         icd = ICHAR('/'), ics = ICHAR(' ')
     
     k = LEN_TRIM(s)
     l = 1
     n = 0
     halo = .FALSE.
     DO WHILE(l<=k.AND.n<SIZE(j))
       DO WHILE(l<=k)
         ic = ICHAR(s(l:l))
         IF( ic == ics ) THEN
           CONTINUE
         ELSE IF( ic == icd ) THEN
           halo = .TRUE.
         ELSE
           EXIT
         END IF
         l=l+1
       END DO
       IF(l>k) EXIT
       IF(.NOT.(ic==icm .OR. ic>=ic0 .AND. ic<=ic9)) EXIT
       
       m = l+1
       DO WHILE(m<=k)
         ic = ICHAR(s(m:m))
         IF(ic<ic0 .OR. ic>ic9) EXIT
         m=m+1
       END DO
       
       n = n + 1
       j(n) = s2i(s(l:m-1),m-l)
       l = m
     END DO
   END FUNCTION read_ints
   

   !---------------------------------------------------
   ! Read header file and allocate some mesh structures
   !---------------------------------------------------
   SUBROUTINE ReadHeaderFile()

     INTEGER :: TypeCount
     INTEGER :: Types(64),CountByType(64)

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//I2S(numprocs)//&
           '/part.'//I2S(mype+1)//'.header'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.header'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadHeaderFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadHeaderFile','Reading header info from file: '//TRIM(FileName),Level=10)
     END IF

     READ(FileUnit,*,IOSTAT=iostat) Mesh % NumberOfNodes, &
         Mesh % NumberOfBulkElements,&
         Mesh % NumberOfBoundaryElements
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadHeaderFile','Could not read header 1st line in file: '//TRIM(FileName))
     END IF

     Types = 0
     CountByType = 0
     READ(FileUnit,*,IOSTAT=iostat) TypeCount
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadHeaderFile','Could not read the type count in file: '//TRIM(FileName))
     END IF
     DO i=1,TypeCount
       READ(FileUnit,*,IOSTAT=iostat) Types(i),CountByType(i)
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadHeaderFile','Could not read type count '&
             //I2S(i)//'in file: '//TRIM(FileName))
       END IF
     END DO

     IF( Parallel ) THEN
       READ(FileUnit,*,IOSTAT=iostat) SharedNodes
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadHeaderFile','Could not read shared nodes in file: '//TRIM(FileName))
       END IF
     ELSE
       SharedNodes = 0
     END IF

     Mesh % MaxElementNodes = 0
     DO i=1,TypeCount
       Mesh % MaxElementNodes = MAX( &
           Mesh % MaxElementNodes, MODULO( Types(i), 100) )
     END DO

     CLOSE(FileUnit)

   END SUBROUTINE ReadHeaderFile


   !-----------------------------------------------------------------------
   ! Read nodes file and create nodal permutation if needed
   !-----------------------------------------------------------------------
   SUBROUTINE ReadNodesFile()

     !USE iso_c_binding
     REAL(c_double) :: Coords(3)
     REAL(c_float) :: SCoords(3)
     INTEGER :: NodeTag
     LOGICAL :: Binary, singlePrec

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//I2S(numprocs)//&
           '/part.'//I2S(mype+1)//'.nodes'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.nodes'
     END IF

     Binary = .FALSE.
     SinglePrec = .FALSE.
     
     OPEN( Unit=FileUnit, File=FileName, STATUS='old', ACTION='read', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       ! ascii file was not successfull, try with binary.
       Binary = .TRUE.
       OPEN( Unit=FileUnit, File=TRIM(FileName)//".bin", FORM='unformatted', &
           ACCESS = 'stream', STATUS='old', ACTION='read', IOSTAT = iostat )
       IF(iostat /= 0 ) THEN         
         SinglePrec = .TRUE.
         OPEN( Unit=FileUnit, File=TRIM(FileName)//".sbin", FORM='unformatted', &
             ACCESS = 'stream', STATUS='old', ACTION='read', IOSTAT = iostat )
       END IF
     END IF
     
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadNodesFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadNodesFile','Reading nodes from file: '//TRIM(FileName),Level=10)
     END IF

     ALLOCATE( NodeTags(Mesh % NumberOfNodes ) ) 
     NodeTags = 0

     NodePermutation = .FALSE.
     DO j = 1, Mesh % NumberOfNodes
       IF(SinglePrec) THEN
         READ(FileUnit,IOSTAT=iostat) NodeTag, SCoords
         Coords = SCoords
       ELSE IF(Binary) THEN
         READ(FileUnit,IOSTAT=iostat) NodeTag, Coords
       ELSE
         READ(FileUnit,*,IOSTAT=iostat) NodeTag, k, Coords
       END IF
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadNodesFile','Problem load node '//I2S(j)//' in file: '//TRIM(Filename))
       END IF

       IF( NodeTags(j) /= j ) NodePermutation = .TRUE. 
       NodeTags(j) = NodeTag
       
       Mesh % Nodes % x(j) = Coords(1)
       Mesh % Nodes % y(j) = Coords(2)
       Mesh % Nodes % z(j) = Coords(3)
     END DO

     CLOSE(FileUnit)

   END SUBROUTINE ReadNodesFile


   !------------------------------------------------------------------------------
   ! Read elements file and create elemental permutation if needed 
   !------------------------------------------------------------------------------
   SUBROUTINE ReadElementsFile()
     TYPE(Element_t), POINTER :: Element
     INTEGER :: ElemType, Tag, Body, ElemNo, Ivals(64),nread, ioffset, partn
     CHARACTER(256) :: str
     LOGICAL :: halo, Binary 


     CALL AllocateVector( ElementTags, Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements, 'ReadElementsFile')   
     ElementTags = 0
     ElementPermutation = .FALSE.

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)// &
          '/partitioning.'//I2S(numprocs)//&
             '/part.'//I2S(mype+1)//'.elements'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.elements'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='old', iostat=IOSTAT )
     IF( iostat == 0 ) THEN
       Binary = .FALSE.
     ELSE
       ! ascii file was not successfull, try with binary.
       Binary = .TRUE.       
       OPEN( Unit=FileUnit, File=TRIM(FileName)//".bin", FORM='unformatted', &
           ACCESS = 'stream', STATUS='old', ACTION='read', IOSTAT = iostat )
     END IF

     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadElementsFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadElementsFile','Reading bulk elements from file: '//TRIM(FileName),Level=10)
     END IF


     DO j=1,Mesh % NumberOfBulkElements

       Element => Mesh % Elements(j)
       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('ReadElementsFile','Element '//I2S(i)//' not associated!')
       END IF

       IF(Binary) THEN
         READ(FileUnit,IOSTAT=iostat) Tag, PartN, body, elemtype
       ELSE
         READ(FileUnit, '(a)', IOSTAT=iostat) str
         IF( iostat /= 0 ) THEN
           CALL Fatal('ReadElementsFile','Could not read start of element entry: '//I2S(j))
         END IF

         nread = read_ints(str,ivals,halo)         
         tag = ivals(1)

         IF( halo ) THEN
           ioffset = 1
           partn = ivals(2) 
         ELSE
           ioffset = 0
           partn = 0 
         END IF
         body = ivals(ioffset+2)
         ElemType = ivals(ioffset+3)
       END IF
         
       ElementTags(j) = tag
       IF( j /= tag ) ElementPermutation = .TRUE.             
       Element % ElementIndex = j
       Element % BodyId = body

       IF( partn > 0 ) THEN
         Element % PartIndex = partn-1
       ELSE
         Element % PartIndex = mype
       END IF

       Element % TYPE => GetElementType(ElemType)

       IF ( .NOT. ASSOCIATED(Element % TYPE) ) THEN
         CALL Fatal('ReadElementsFile','Element of type '&
             //I2S(ElemType)//' could not be associated!')
       END IF

       n = Element % TYPE % NumberOfNodes
       CALL AllocateVector( Element % NodeIndexes, n )

       IF( Binary ) THEN
         READ(FileUnit,IOSTAT=iostat) Element % NodeIndexes(1:n)
       ELSE
         IF( nread < n + ioffset + 3 ) THEN
           CALL Fatal('ReadElementsFile','Line '//I2S(j)//' does not contain enough entries')
         END IF
         Element % NodeIndexes(1:n) = IVals(4+ioffset:nread)
       END IF
     END DO
     CLOSE( FileUnit ) 

   END SUBROUTINE ReadElementsFile
   !------------------------------------------------------------------------------


   !------------------------------------------------------------------------------
   ! Read boundary elements file and remap the parents if needed.  
   !------------------------------------------------------------------------------
   SUBROUTINE ReadBoundaryFile( BoundariesOnly )
     LOGICAL, OPTIONAL :: BoundariesOnly
     INTEGER, POINTER :: LocalEPerm(:)
     INTEGER :: MinEIndex, MaxEIndex, ElemNodes, i, j
     INTEGER :: Left, Right, bndry, tag, ElemType, IVals(64), nread, ioffset, partn, nswap
     TYPE(Element_t), POINTER :: Element
     CHARACTER(256) :: str
     LOGICAL :: halo, Binary, BOnly

     BOnly = .FALSE.
     IF ( PRESENT(BoundariesOnly) ) BOnly=BoundariesOnly

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//I2S(numprocs)//&
           '/part.'//I2S(mype+1)//'.boundary'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.boundary'
     END IF

     ! Create permutation for the elements. This is needed when the element 
     ! parents are mapped to the new order. This is needed for mapping of the 
     ! parents. Otherwise the element numbering is arbitrary. 
     !------------------------------------------------------------------------------
     IF( ElementPermutation ) THEN
       MinEIndex = MINVAL( ElementTags(1:Mesh % NumberOfBulkElements) )
       MaxEIndex = MAXVAL( ElementTags(1:Mesh % NumberOfBulkElements) )

       LocalEPerm => NULL()
       CALL AllocateVector( LocalEPerm, MaxEIndex - MinEIndex + 1, 'ReadBoundaryFile' )
       LocalEPerm = 0
       DO i=1,Mesh % NumberOfBulkElements
         LocalEPerm( ElementTags(i) - MinEIndex + 1 ) = i
       END DO
     ELSE
       MinEIndex = 1 
       MaxEIndex = Mesh % NumberOfBulkElements
     END IF


     OPEN( Unit=FileUnit, File=FileName, STATUS='old', iostat=IOSTAT )

     IF( iostat == 0 ) THEN
       Binary = .FALSE.
     ELSE
       ! ascii file was not successfull, try with binary.
       Binary = .TRUE.       
       OPEN( Unit=FileUnit, File=TRIM(FileName)//".bin", FORM='unformatted', &
           ACCESS = 'stream', STATUS='old', ACTION='read', IOSTAT = iostat )
     END IF

     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadBoundaryFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadBoundaryFile','Reading boundary elements from file: '//TRIM(FileName),Level=10)
     END IF

     nswap = 0
     
     DO j=Mesh % NumberOfBulkElements+1, &
         Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements

       Element => Mesh % Elements(j)
       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('ReadBoundaryFile','Element '//I2S(i)//' not associated!')
       END IF

       IF(Binary) THEN
         READ(FileUnit,IOSTAT=iostat) Tag, PartN, bndry, left, right, elemtype
       ELSE
         READ(FileUnit, '(a)', IOSTAT=iostat) str
         IF( iostat /= 0 ) THEN
           CALL Fatal('ReadBoundaryFile','Could not read boundary element entry: '//I2S(j))
         END IF
         nread = read_ints(str,ivals,halo)
         
         tag = ivals(1)
         ElementTags(j) = tag
         
         IF( halo ) THEN
           partn = ivals(2)
           ioffset = 1
         ELSE
           partn = 0
           ioffset = 0
         END IF
         
         bndry = ivals(ioffset+2)
         left = ivals(ioffset+3)
         right = ivals(ioffset+4)
         ElemType = ivals(ioffset+5)
       END IF
         
       Element % ElementIndex = j
       Element % TYPE => GetElementType(ElemType)
       IF ( .NOT. ASSOCIATED(Element % TYPE) ) THEN
         CALL Fatal('ReadBoundaryFile','Element of type '//I2S(ElemType)//'could not be associated!')
       END IF

       ElemNodes = Element % TYPE % NumberOfNodes
       Mesh % MaxElementNodes = MAX( Mesh % MaxElementNodes, ElemNodes )

       IF( partn == 0 ) THEN
         Element % PartIndex = mype
       ELSE
         Element % PartIndex = partn-1
       END IF

       CALL AllocateBoundaryInfo( Element ) 

       Element % BoundaryInfo % Constraint = bndry
       Element % BoundaryInfo % Left => NULL()
       Element % BoundaryInfo % Right => NULL()

       IF ( Left >= MinEIndex .AND. Left <= MaxEIndex ) THEN
         IF( ElementPermutation ) THEN
           Left  = LocalEPerm(Left - MinEIndex + 1)
         END IF
       ELSE IF ( Left>0 .AND. .NOT. BOnly ) THEN
         WRITE( Message, * ) mype,'BOUNDARY PARENT out of range: ', Tag, Left
         CALL Error( 'ReadBoundaryFile', Message )
         Left = 0
       END IF

       IF ( Right >= MinEIndex .AND. Right <= MaxEIndex ) THEN
         IF( ElementPermutation ) THEN
           Right = LocalEPerm(Right - MinEIndex + 1)
         END IF
       ELSE IF ( Right > 0 .AND. .NOT. BOnly ) THEN
         WRITE( Message, * ) mype,'BOUNDARY PARENT out of range: ', Tag,Right
         CALL Error( 'ReadBoundaryFile', Message )
         Right = 0
       END IF

       ! We always want to have "Left" present.
       IF( Right >= 1 .AND. Left == 0 ) THEN
         Left = Right
         Right = 0
         nswap = nswap + 1
       END IF
       
       IF ( Left >= 1 ) THEN
         Element % BoundaryInfo % Left => Mesh % Elements(left)
       END IF

       IF ( Right >= 1 ) THEN
         Element % BoundaryInfo % Right => Mesh % Elements(right)
       END IF

       n = Element % TYPE % NumberOfNodes
       CALL AllocateVector( Element % NodeIndexes, n )

       IF( binary ) THEN
         READ(FileUnit,IOSTAT=iostat) Element % NodeIndexes(1:n)
       ELSE
         IF( nread < 5 + n + ioffset ) THEN
           CALL Fatal('ReadBoundaryFile','Line '//I2S(j)//' does not contain enough entries')
         END IF
         Element % NodeIndexes(1:n) = Ivals(6+ioffset:nread)
       END IF
     END DO
     CLOSE( FileUnit )

     IF(nswap > 0) THEN
       CALL Info('ReadBoundaryFile',&
           'Swapped '//I2S(nswap)//' "right" owners to "left" to always have left parent existing!')
     END IF
            
     IF( ElementPermutation ) THEN
       DEALLOCATE( LocalEPerm ) 
     END IF

   END SUBROUTINE ReadBoundaryFile
   !------------------------------------------------------------------------------



   ! Make a permutation for the bulk and boundary element topology if 
   ! the nodes are permuted. This is always the case in parallel.
   ! The initial numbering is needed only when the nodes are loaded and 
   ! hence this is a local subroutine. 
   !----------------------------------------------------------------------
   SUBROUTINE PermuteNodeNumbering()

     TYPE(Element_t), POINTER :: Element

     IF( NodePermutation ) THEN
       CALL Info('PermuteNodeNumbering','Performing node mapping',Level=6)

       MinNodeTag = MINVAL( NodeTags )
       MaxNodeTag = MAXVAL( NodeTags )

       CALL AllocateVector( LocalPerm, MaxNodeTag-MinNodeTag+1, 'PermuteNodeNumbering' )
       LocalPerm = 0
       DO i=1,Mesh % NumberOfNodes
         LocalPerm(NodeTags(i) - MinNodeTag + 1) = i
       END DO

       DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements       
         Element => Mesh % Elements(i)
         n = Element % TYPE % NumberOfNodes

         DO j=1,n
           k = Element % NodeIndexes(j) 
           Element % NodeIndexes(j) = LocalPerm(k - MinNodeTag + 1)
         END DO
       END DO
     ELSE
       CALL Info('PermuteNodeNumbering','Node mapping is continuous',Level=8)
     END IF

     ! Set the for now, if the case is truly parallel we'll have to revisit these
     ! when reading the parallel information. 
     Mesh % ParallelInfo % NumberOfIfDOFs = 0
     Mesh % ParallelInfo % GlobalDOFs => NodeTags

   END SUBROUTINE PermuteNodeNumbering


   ! Initialize some parallel structures once the non-nodal 
   ! element types are known. 
   ! Currently this is here mainly because the 
   ! Elemental and Nodal tags are local
   !-------------------------------------------------------
   SUBROUTINE InitParallelInfo()

     INTEGER, POINTER :: TmpGlobalDofs(:)

     ! These two have already been set, and if the case is serial
     ! case they can be as is.
     !Mesh % ParallelInfo % NumberOfIfDOFs = 0
     !Mesh % ParallelInfo % GlobalDOFs => NodeTags


     ! This also for serial runs ...
     DO i=1,Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
       Mesh % Elements(i) % GElementIndex = ElementTags(i)
     END DO

     IF(.NOT. Parallel ) RETURN

     n = Mesh % NumberOfNodes + &
         Mesh % MaxEdgeDOFs * Mesh % NumberOFEdges + &
         Mesh % MaxFaceDOFs * Mesh % NumberOFFaces + &
         Mesh % MaxBDOFs    * Mesh % NumberOFBulkElements

     ALLOCATE( TmpGlobalDOFs(n) )
     TmpGlobalDOFs = 0
     TmpGlobalDOFs(1:Mesh % NumberOfNodes) = &
         Mesh % ParallelInfo % GlobalDOFs(1:Mesh % NumberOfNodes)
     DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs ) 
     Mesh % ParallelInfo % GlobalDofs => TmpGlobalDofs

     ALLOCATE(Mesh % ParallelInfo % NeighbourList(n), STAT=istat)
     IF (istat /= 0) CALL Fatal('InitParallelInfo', 'Unable to allocate NeighbourList array.')

     DO i=1,n
       NULLIFY( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
     END DO

     CALL AllocateVector( Mesh % ParallelInfo % GInterface, n, 'InitParallelInfo')
     Mesh % ParallelInfo % GInterface = .FALSE.       

   END SUBROUTINE InitParallelInfo


   ! Read the file that shows the shared nodes.
   !------------------------------------------------------------------------
   SUBROUTINE ReadSharedFile()

     INTEGER :: Ivals(64)
     INTEGER :: npart, tag, nread
     CHARACTER(256) :: str
     LOGICAL :: halo

     IF(.NOT. Parallel) RETURN

     FileName = BaseName(1:BaseNameLen)//&
       '/partitioning.'//I2S(numprocs)//&
         '/part.'//I2S(mype+1)//'.shared'

     OPEN( Unit=FileUnit, File=FileName, STATUS='old', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadSharedFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadSharedFile','Reading nodes from file: '//TRIM(FileName),Level=10)
     END IF

     ! This loop could be made more effective, for example
     ! by reading tags and nparts to a temporal vector
     ! The operation using the str takes much more time.
     !-----------------------------------------------------
     DO i=1,SharedNodes          
       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadSharedFile','Could not read shared nodes entry: '//I2S(i))
       END IF
       nread = read_ints(str,ivals,halo)

       tag = ivals(1)
       npart = ivals(2)       

       k = LocalPerm( tag-MinNodeTag+1 )
       Mesh % ParallelInfo % GInterface(k) = .TRUE.
       CALL AllocateVector(Mesh % ParallelInfo % NeighbourList(k) % Neighbours,npart)

       IF( nread < 2 + npart ) THEN
         CALL Fatal('ReadSharedFile','Line '//I2S(j)//' does not contain enough entries')
       END IF
       
       Mesh % ParallelInfo % NeighbourList(k) % Neighbours = ivals(3:nread) - 1

       ! this partition does not own the node
       IF ( ivals(3)-1 /= mype ) THEN
         Mesh % ParallelInfo % NumberOfIfDOFs = &
             Mesh % ParallelInfo % NumberOfIfDOFs + 1
       END IF
     END DO

     CLOSE( FileUnit )

   END SUBROUTINE ReadSharedFile


   ! Initialize parallel info for pseudo parallel meshes
   !-------------------------------------------------------
   SUBROUTINE InitPseudoParallel()

     INTEGER, POINTER :: TmpGlobalDofs(:)

     ! This also for serial runs ...
     n = ParEnv % MyPe * Mesh % NumberOfBulkElements

     DO i=1,Mesh % NumberOfBulkElements
       Mesh % Elements(i) % GElementIndex = ElementTags(i) + n
     END DO

     n = Mesh % NumberOfNodes + &
         Mesh % MaxEdgeDOFs * Mesh % NumberOFEdges + &
         Mesh % MaxFaceDOFs * Mesh % NumberOFFaces + &
         Mesh % MaxBDOFs    * Mesh % NumberOFBulkElements

     ALLOCATE( TmpGlobalDOFs(n) )
     TmpGlobalDOFs = 0
     TmpGlobalDOFs(1:Mesh % NumberOfNodes) = &
         Mesh % ParallelInfo % GlobalDOFs(1:Mesh % NumberOfNodes) + n
     DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs ) 
     Mesh % ParallelInfo % GlobalDofs => TmpGlobalDofs
     
     ALLOCATE(Mesh % ParallelInfo % NeighbourList(n), STAT=istat)
     IF (istat /= 0) CALL Fatal('InitParallelInfo', 'Unable to allocate NeighbourList array.')
     
     DO i=1,n
       ALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) )
       Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) = ParEnv % MyPe
     END DO

     CALL AllocateVector( Mesh % ParallelInfo % GInterface, n, 'InitParallelInfo')
     Mesh % ParallelInfo % GInterface = .FALSE.       

   END SUBROUTINE InitPseudoParallel

   
 END SUBROUTINE ElmerMeshReader

 !> An interface over potential mesh loading strategies. 
 !----------------------------------------------------------------- 
 SUBROUTINE LoadMeshStep( Step, PMesh, MeshNamePar, ThisPe, NumPEs, &
         IsParallel, BoundariesOnly ) 
   
   IMPLICIT NONE

   INTEGER :: Step
   CHARACTER(LEN=*), OPTIONAL :: MeshNamePar
   TYPE(Mesh_t), POINTER, OPTIONAL :: PMesh
   INTEGER, OPTIONAL :: ThisPe, NumPEs
   LOGICAL, OPTIONAL :: IsParallel
   LOGICAL, OPTIONAL :: BoundariesOnly

   ! Currently only one strategy to get the mesh is implemented 
   ! but there could be others.
   !
   ! This has not yet been tested in parallel and for sure
   ! it does not work for halo elements. 
   !-----------------------------------------------------------------
   CALL ElmerMeshReader( Step, PMesh, MeshNamePar, ThisPe, NumPEs, &
           IsParallel, BoundariesOnly ) 

 END SUBROUTINE LoadMeshStep


 ! Set the mesh dimension by studying the coordinate values.
 ! This could be less conservative also...
 !------------------------------------------------------------------------------    
 SUBROUTINE SetMeshDimension( Mesh )
   TYPE(Mesh_t) :: Mesh
   
   REAL(KIND=dp) :: x, y, z
   LOGICAL :: C(3)
   INTEGER :: i
   
   IF( Mesh % NumberOfNodes == 0 ) RETURN

   ! Compare value to some node, why not the 1st one
   x = Mesh % Nodes % x(1)
   y = Mesh % Nodes % y(1)
   z = Mesh % Nodes % z(1)
   
   C(1) = ANY( Mesh % Nodes % x /= x ) 
   C(2) = ANY( Mesh % Nodes % y /= y )  
   C(3) = ANY( Mesh % Nodes % z /= z )  

   ! This version is perhaps too liberal 
   Mesh % MeshDim = COUNT( C )
   Mesh % MaxDim = 0
   DO i=1,3
     IF( C(i) ) Mesh % MaxDim = i
   END DO
      
   CALL Info('SetMeshDimension','Dimension of mesh is: '//I2S(Mesh % MeshDim),Level=8)
   CALL Info('SetMeshDimension','Max dimension of mesh is: '//I2S(Mesh % MaxDim),Level=8)

 END SUBROUTINE SetMeshDimension

 
 SUBROUTINE ReadTargetNames(Model,Filename)
   CHARACTER(LEN=*) :: FileName
   TYPE(Model_t) :: Model
!------------------------------------------------------------------------------
   INTEGER, PARAMETER :: FileUnit = 10
   INTEGER, PARAMETER :: A=ICHAR('A'),Z=ICHAR('Z'),U2L=ICHAR('a')-ICHAR('A')
   INTEGER :: i,j,k,iostat,i1,i2,i3,n
   INTEGER :: ivals(256)
   CHARACTER(LEN=1024) :: str, name0, name1
   TYPE(ValueList_t), POINTER :: Vlist
   LOGICAL :: Found, AlreadySet, DoIt, DoBCs, DoBodies
   INTEGER :: BodyMaps, BCMaps
   CHARACTER(*), PARAMETER :: Caller = 'ReadTargetNames'

   
   DoIt = ListGetLogical( Model % Simulation,'Use Mesh Names',Found )
   IF(DoIt) THEN   
     DoBCs = .TRUE.
     DoBodies = .TRUE.
   ELSE     
     DoBCs = .FALSE.
     DoBodies = .FALSE.
   END IF

   DoIt = ListGetLogical( Model % Simulation,'Use Mesh Body Names',Found )
   IF(Found) DoBodies = DoIt   
   DoIt = ListGetLogical( Model % Simulation,'Use Mesh Boundary Names',Found ) 
   IF(Found) DoBCs = DoIt

   IF(.NOT. (DoBodies .OR. DoBCs )) RETURN
   
   BodyMaps = 0
   BCMaps = 0

   OPEN( Unit=FileUnit, File=FileName, STATUS='old', IOSTAT=iostat )
   IF( iostat /= 0 ) THEN
     CALL Fatal(Caller,'Requested the use of entity names but this file does not exits: '//TRIM(FileName))
   END IF
   
   CALL Info(Caller,'Reading names info from file: '//TRIM(FileName),Level=10)

   DO WHILE( .TRUE. ) 
     READ(FileUnit,'(A)',IOSTAT=iostat) str
     IF( iostat /= 0 ) EXIT
     i = INDEX( str,'$')     
     j = INDEX( str,'=')
     IF( i == 0 .OR. j == 0 ) CYCLE

     i = i + 1
     DO WHILE(i<=LEN_TRIM(str) .AND. str(i:i)==' ')
       i = i + 1
     END DO     
     
     i1 = i
     i2 = j-1
     i3 = j+1

     ! Move to lowercase since the "name" in sif file is also
     ! always in lowercase. 
     DO i=i1,i2
       j = i+1-i1
       k = ICHAR(str(i:i))
       IF ( k >= A .AND. k<= Z ) THEN
         name0(j:j) = CHAR(k+U2L)
       ELSE
         name0(j:j) = str(i:i)
       END IF
     END DO

     n = str2ints( str(i3:),ivals )
     IF( n == 0 ) THEN
       CALL Fatal(Caller,'Could not find arguments for: '//str(i1:i2))
     END IF

     AlreadySet = .FALSE.

     DO i=1,Model % NumberOfBCs
       IF(.NOT. DoBCs) CYCLE
       Vlist => Model % BCs(i) % Values
       name1 = ListGetString( Vlist,'Name',Found )
       IF(.NOT. Found ) CYCLE
       IF( name0(1:i2-i1+1) == TRIM(name1) ) THEN
!        PRINT *,'Name > '//TRIM(name1)//' < matches BC '//I2S(i)
         IF( AlreadySet ) THEN
           CALL Fatal(Caller,'Mapping of name is not unique: '//TRIM(name1) )
         ELSE IF( ListCheckPresent( Vlist,'Target Boundaries') ) THEN
           CALL Info(Caller,'> Target Boundaries < already defined for BC '//I2S(i))
         ELSE
           CALL ListAddIntegerArray( Vlist,'Target Boundaries',n,ivals(1:n))
           BodyMaps = BodyMaps + 1
           AlreadySet = .TRUE.
         END IF
       END IF
     END DO

     DO i=1,Model % NumberOfBodies
       IF(.NOT. DoBodies) CYCLE
       Vlist => Model % Bodies(i) % Values
       name1 = ListGetString( Vlist,'Name',Found )
       IF(.NOT. Found ) CYCLE
       IF( name0(1:i2-i1+1) == TRIM(name1) ) THEN
!        PRINT *,'Name > '//TRIM(name1)//' < matches body '//I2S(i)
         IF( AlreadySet ) THEN
           CALL Fatal(Caller,'Mapping of name is not unique: '//TRIM(name1) )
         ELSE IF( ListCheckPresent( Vlist,'Target Bodies') ) THEN
           CALL Info(Caller,'> Target Bodies < already defined for Body '//I2S(i))
         ELSE
           CALL ListAddIntegerArray( Vlist,'Target Bodies',n,ivals(1:n))
           BCMaps = BCMaps + 1
           AlreadySet = .TRUE.
         END IF
       END IF
     END DO
     
     IF(.NOT. AlreadySet ) THEN
       IF( ParEnv % MyPe == 0 ) THEN
         CALL Info(Caller,'Could not map name to Body nor BC: '//name0(1:i2-i1+1), Level=20)
       END IF
     END IF

   END DO
   CLOSE(FileUnit)
      
   CALL Info(Caller,'Mapped '//I2S(BodyMaps)//' body names and '//I2S(BCMaps)//' bc names to elements!')
     
 END SUBROUTINE ReadTargetNames


!------------------------------------------------------------------------------
!> This subroutine reads elementwise input data from the file mesh.elements.data 
!> and inserts the data into the structured data variable 
!> Mesh % Elements(element_id) % PropertyData. The contents of the file should
!> be arranged as
!> 
!> element: element_id_1
!> data_set_name_1: a_1 a_2 ... a_n
!> data_set_name_2: b_1 b_2 ... b_m
!> data_set_name_3: ...
!> end
!> element: ...
!> ...
!> end
!------------------------------------------------------------------------------
  SUBROUTINE ReadElementPropertyFile(FileName,Mesh)
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: FileName
     TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    CHARACTER(LEN=:), ALLOCATABLE :: str
    INTEGER :: i,j,n
    INTEGER, PARAMETER :: FileUnit = 10
    REAL(KIND=dp) :: x
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementData_t), POINTER :: PD,PD1
!------------------------------------------------------------------------------
    OPEN( Unit=FileUnit, File=FileName, STATUS='old', ERR=10 )

    ALLOCATE(CHARACTER(MAX_STRING_LEN)::str)
    DO WHILE( ReadAndTrim(FileUnit,str) )
      READ( str(9:),*) i
      IF ( i < 0 .OR. i > Mesh % NumberOFBulkElements ) THEN
        CALL Fatal( 'ReadElementPropertyFile', 'Element id out of range.' )
      END IF

      IF ( SEQL( str, 'element:') ) THEN
        Element => Mesh % Elements(i)
        PD => Element % PropertyData

        DO WHILE(ReadAndTrim(FileUnit,str))
          IF ( str == 'end' ) EXIT

          i = INDEX(str, ':')
          IF ( i<=0 ) CYCLE

          IF ( .NOT.ASSOCIATED(PD)  ) THEN
            ALLOCATE( Element % PropertyData )
            PD => Element % PropertyData
            PD % Name = TRIM(str(1:i-1))
          ELSE
            DO WHILE(ASSOCIATED(PD))
              IF ( PD % Name==TRIM(str(1:i-1)) ) EXIT
              PD1 => PD
              PD => PD % Next
            END DO
            
            IF (.NOT. ASSOCIATED(PD) ) THEN
              ALLOCATE(PD1 % Next)
              PD => PD1 % Next
              PD % Name = TRIM(str(1:i-1))
            END IF
          END IF

          j = i+1
          n = 0
          DO WHILE(j<=LEN_TRIM(str))
            READ( str(j:), *, END=20,ERR=20 ) x
            n = n + 1
            DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)==' ')
              j = j + 1
            END DO
            DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)/=' ')
              j = j + 1
            END DO
          END DO
20        CONTINUE
          IF ( n>0 ) THEN
            ALLOCATE(PD % Values(n))
            j = i+1
            n = 1
            DO WHILE(j<=LEN_TRIM(str))
              READ( str(j:), *, END=30,ERR=30 ) PD % Values(n)
              n = n + 1
              DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)==' ')
                j = j + 1
              END DO
              DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)/=' ')
                j = j + 1
              END DO
            END DO
30          CONTINUE
          END IF
        END DO
      END IF
    END DO

    CLOSE(FileUnit)

10  CONTINUE

!------------------------------------------------------------------------------
  END SUBROUTINE ReadElementPropertyFile
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Writes the mesh to disk. Note that this does not include the information
!> of shared nodes needed in parallel computation. This may be used for 
!> debugging purposes and for adaptive solution, for example. 
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDisk( NewMesh, Path )
!------------------------------------------------------------------------------
    CHARACTER(LEN=*) :: Path
    TYPE(Mesh_t) :: NewMesh
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,MaxNodes,ElmCode,Parent1,Parent2
!------------------------------------------------------------------------------

    OPEN( 1,FILE=TRIM(Path) // '/mesh.header',STATUS='UNKNOWN' )
    WRITE( 1,'(i0,x,i0,x,i0)' ) NewMesh % NumberOfNodes, &
         NewMesh % NumberOfBulkElements, NewMesh % NumberOfBoundaryElements
    
    WRITE( 1,'(i0)' ) 2
    MaxNodes = 0
    ElmCode  = 0
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       IF ( NewMesh % Elements(k) % TYPE % NumberOfNodes > MaxNodes ) THEN
          ElmCode  = NewMesh % Elements(k) % TYPE % ElementCode
          MaxNodes = NewMesh % Elements(k) % TYPE % NumberOfNodes
       END IF
    END DO
    WRITE( 1,'(i0,x,i0)' ) ElmCode,NewMesh % NumberOfBoundaryElements

    MaxNodes = 0
    ElmCode  = 0
    DO i=1,NewMesh % NumberOfBulkElements
       IF ( NewMesh % Elements(i) % TYPE % NumberOfNodes > MaxNodes ) THEN
          ElmCode  = NewMesh % Elements(i) % TYPE % ElementCode
          MaxNodes = NewMesh % Elements(i) % TYPE % NumberOfNodes
       END IF
    END DO
    WRITE( 1,'(i0,x,i0)' ) ElmCode,NewMesh % NumberOfBulkElements
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.nodes', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       WRITE(1,'(i0,a,3e23.15)',ADVANCE='NO') i,' -1 ', &
            NewMesh % Nodes % x(i), &
            NewMesh % Nodes % y(i), NewMesh % Nodes % z(i)
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.elements', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBulkElements
       WRITE(1,'(3(i0,x))',ADVANCE='NO') i, &
            NewMesh % Elements(i) % BodyId, &
            NewMesh % Elements(i) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(i) % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') &
               NewMesh % Elements(i) % NodeIndexes(j)
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.boundary', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       parent1 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Left ) ) &
          parent1 = NewMesh % Elements(k) % BoundaryInfo % Left % ElementIndex
       parent2 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Right ) ) &
          parent2 = NewMesh % Elements(k) % BoundaryInfo % Right % ElementIndex
       WRITE(1,'(5(i0,x))',ADVANCE='NO') i, &
            NewMesh % Elements(k) % BoundaryInfo % Constraint, Parent1,Parent2,&
            NewMesh % Elements(k) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(k) % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') &
               NewMesh % Elements(k) % NodeIndexes(j)
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)
!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDisk
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Writes the mesh to disk, including detection of elementcodes and shared node
!> info necessary for parallel meshes.
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDisk2(Model, NewMesh, Path, Partition )
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t) :: NewMesh
    CHARACTER(LEN=*) :: Path
    INTEGER, OPTIONAL :: Partition
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,m,MaxNodes,ElmCode,NumElmCodes,ElmCodeList(100),ElmCodeCounts(100),&
        Parent1,Parent2, ElemID, nneigh, Constraint, meshBC, NumElements, NoShared, &
        iostat, BCWarns
    INTEGER, POINTER :: BList(:)
    INTEGER, ALLOCATABLE :: ElementCodes(:)
    LOGICAL :: Parallel, WarnNoTarget, Found
    CHARACTER(:), ALLOCATABLE :: headerFN, elementFN, nodeFN,&
         boundFN, sharedFN
!------------------------------------------------------------------------------

    IF(PRESENT(Partition)) THEN
       Parallel = .TRUE.
       headerFN = '/part.'//I2S(Partition+1)//'.header'
       elementFN = '/part.'//I2S(Partition+1)//'.elements'
       nodeFN =  '/part.'//I2S(Partition+1)//'.nodes'
       boundFN = '/part.'//I2S(Partition+1)//'.boundary'
       sharedFN ='/part.'//I2S(Partition+1)//'.shared'
    ELSE
       Parallel = .FALSE.
       headerFN = '/mesh.header'
       elementFN = '/mesh.elements'
       nodeFN = '/mesh.nodes'
       boundFN = '/mesh.boundary'
    END IF

    !Info for header file

    ElmCodeList = 0 !init array
    NumElmCodes = 0
    NumElements = NewMesh % NumberOfBoundaryElements + &
         NewMesh % NumberOfBulkElements
    ALLOCATE(ElementCodes(NumElements))

    !cycle to bring element code list into array-inquirable form
    DO i=1,NumElements
       ElementCodes(i) = NewMesh % Elements(i) % TYPE % ElementCode
    END DO

    DO i=NumElements,1,-1 !this should give element codes increasing value, which appears to be
                          !'standard' though I doubt it matters
       IF(ANY(ElmCodeList == ElementCodes(i))) CYCLE
       NumElmCodes = NumElmCodes + 1
       ElmCodeList(NumElmCodes) = ElementCodes(i)
    END DO

    DO j=1,NumElmCodes
       ElmCodeCounts(j) = COUNT(ElementCodes == ElmCodeList(j))
    END DO

    !Write header file
    OPEN( 1,FILE=TRIM(Path) // headerFN,STATUS='UNKNOWN', iostat=iostat)
    IF(iostat /= 0) THEN
      CALL Fatal('WriteMeshToDisk2','Could not open file: '//TRIM(Path)//headerFN)
    END IF

    WRITE( 1,'(i0,x,i0,x,i0)' ) NewMesh % NumberOfNodes, &
         NewMesh % NumberOfBulkElements, &
         NewMesh % NumberOfBoundaryElements

    WRITE( 1,'(i0)' ) NumElmCodes
    DO j=1,NumElmCodes
       WRITE( 1,'(i0,x,i0,x)' ) ElmCodeList(j),ElmCodeCounts(j)
    END DO
    IF(Parallel) THEN !need number of shared nodes
       NoShared = 0
       DO i=1,NewMesh % NumberOfNodes
          IF(SIZE(NewMesh % ParallelInfo % NeighbourList(i) % &
               Neighbours) > 1) THEN
             NoShared = NoShared + 1
          END IF
       END DO
       WRITE( 1,'(i0,x,i0)') NoShared, 0
    END IF
    CLOSE(1)

    !Write nodes file
    OPEN( 1,FILE=TRIM(Path) // nodeFN, STATUS='UNKNOWN',iostat=iostat)
    IF(iostat /= 0) THEN
      CALL Fatal('WriteMeshToDisk2','Could not open file: '//TRIM(Path)//nodeFN)
    END IF
    DO i=1,NewMesh % NumberOfNodes
       IF (Parallel) THEN
          WRITE(1,'(i0,x)', ADVANCE='NO') &
               NewMesh % ParallelInfo % GlobalDOFs(i)
       ELSE
          WRITE(1,'(i0,x)', ADVANCE='NO') i
       END IF
       WRITE(1,'(a,x,ES17.10,x,ES17.10,x,ES17.10)',ADVANCE='NO') &
            ' -1 ', NewMesh % Nodes % x(i), &
            NewMesh % Nodes % y(i), NewMesh % Nodes % z(i)
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)

    !Write elements file
    OPEN( 1,FILE=TRIM(Path) // elementFN, STATUS='UNKNOWN', iostat=iostat)
    IF(iostat /= 0) THEN
      CALL Fatal('WriteMeshToDisk2','Could not open file: '//TRIM(Path)//elementFN)
    END IF
    DO i=1,NewMesh % NumberOfBulkElements
       IF(Parallel) THEN
          ElemID = NewMesh % Elements(i) % GElementIndex
       ELSE
          ElemID = i
       END IF
       WRITE(1,'(i0,x,i0,x,i0,x)',ADVANCE='NO') ElemID, &
            NewMesh % Elements(i) % BodyId, &
            NewMesh % Elements(i) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(i) % TYPE % NumberOfNodes
          IF(Parallel) THEN
             m = NewMesh % ParallelInfo % GlobalDOFs(&
                  NewMesh % Elements(i) % NodeIndexes(j))
          ELSE
             m = NewMesh % Elements(i) % NodeIndexes(j)
          END IF
          WRITE(1,'(i0,x)', ADVANCE='NO') m
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)

    !Write boundary file
    WarnNoTarget = .FALSE.
    OPEN( 1,FILE=TRIM(Path) // boundFN, STATUS='UNKNOWN',iostat=iostat)
    IF(iostat /= 0) THEN
      CALL Fatal('WriteMeshToDisk2','Could not open file: '//TRIM(Path)//boundFN)
    END IF
    BcWarns = 0
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       parent1 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Left ) ) &
          parent1 = NewMesh % Elements(k) % BoundaryInfo % Left % ElementIndex
       parent2 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Right ) ) &
          parent2 = NewMesh % Elements(k) % BoundaryInfo % Right % ElementIndex

       IF(Parallel) THEN
          IF(parent1 /= 0) parent1 = NewMesh % Elements(parent1) % GElementIndex
          IF(parent2 /= 0) parent2 = NewMesh % Elements(parent2) % GElementIndex
       END IF

       IF(.NOT. ASSOCIATED(NewMesh % Elements(k) % BoundaryInfo ) ) THEN
         CALL Fatal('WriteMeshToDisk2','BoundaryInfo not associated for element: '//I2S(k))
       END IF
       
       Constraint = NewMesh % Elements(k) % BoundaryInfo % Constraint

       Found = .FALSE.
       IF(Constraint > 0 .AND. Constraint <= Model % NumberOfBCs ) THEN
         BList => ListGetIntegerArray( Model % BCs(Constraint) % Values, &
             'Target Boundaries', Found )
       END IF
       IF(Found) THEN
          IF(SIZE(BList) > 1) THEN
            BcWarns = BcWarns + 1
          END IF
          meshBC = BList(1)
       ELSE
          WarnNoTarget = .TRUE.
          meshBC = Constraint
       END IF

       !This meshBC stuff will *only* work if each BC has only 1 target boundary
       WRITE(1,'(i0,x,i0,x,i0,x,i0,x,i0)',ADVANCE='NO') i, & 
            meshBC, Parent1,Parent2,&
            NewMesh % Elements(k) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(k) % TYPE % NumberOfNodes
          IF(Parallel) THEN
             m = NewMesh % ParallelInfo % GlobalDOFs(&
                  NewMesh % Elements(k) % NodeIndexes(j))
          ELSE
             m = NewMesh % Elements(k) % NodeIndexes(j)
          END IF
          WRITE(1,'(x,i0)', ADVANCE='NO') m
       END DO
       WRITE(1,*) !blank write statement to create new line without extra space.
    END DO
    CLOSE(1)

    IF(BcWarns > 1 ) THEN
      CALL WARN("WriteMeshToDisk2",&
          "BC elements '//I2S(BcWarns)//' have more than one Target Boundary, SaveMesh output will not match input!")
    END IF
      
    IF(WarnNoTarget) THEN
       CALL WARN("WriteMeshToDisk2","Couldn't find a Target Boundary, assuming mapping to self")
    END IF

    IF(.NOT. Parallel) RETURN

    !Write .shared file
    !Need to create part.n.shared from Mesh % ParallelInfo %
    !NeighbourList % Neighbours.
    OPEN( 1,FILE=TRIM(Path) // sharedFN, STATUS='UNKNOWN',iostat=iostat)
    IF(iostat /= 0) THEN
      CALL Fatal('WriteMeshToDisk2','Could not open file: '//TRIM(Path)//sharedFN)
    END IF
    DO i=1,NewMesh % NumberOfNodes
       nneigh = SIZE(NewMesh % ParallelInfo % NeighbourList(i) % &
            Neighbours)
       IF(nneigh < 2) CYCLE
       WRITE(1,'(i0, x, i0, x)',ADVANCE='NO') &
            NewMesh % ParallelInfo % GlobalDOFs(i),nneigh
       DO j=1,nneigh
          WRITE(1,'(I0, x)',ADVANCE='NO') NewMesh % ParallelInfo %&
               NeighbourList(i) % Neighbours(j) + 1
       END DO
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)


!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDisk2
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Writes the mesh to disk, including detection of elementcodes and shared node
!> info necessary for parallel meshes.
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDiskPartitioned(Model, Mesh, Path, &
      ElementPart, NeighbourList )
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t) :: Mesh
    CHARACTER(LEN=*) :: Path
    INTEGER :: ElementPart(:)
    TYPE(NeighbourList_t), TARGET :: NeighbourList(:)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: NoBoundaryElements, NoBulkElements, NoNodes, NoPartitions, Partition
    INTEGER :: i,j,k,m,MaxNodes,ElmCode,NumElmCodes,ElmCodeCounts(827),&
         Parent1,Parent2, ElemID, nneigh, Constraint, meshBC, NumElements, NoShared
    LOGICAL :: Found, Hit
    CHARACTER(:), ALLOCATABLE :: DirectoryName, PrefixName
!------------------------------------------------------------------------------

    NoPartitions = MAXVAL( ElementPart ) 
    NumElmCodes = 0
    NumElements = Mesh % NumberOfBoundaryElements + Mesh % NumberOfBulkElements
        
    DirectoryName = TRIM(PATH)//'/partitioning.'//I2S(NoPartitions)
    CALL MakeDirectory( DirectoryName // CHAR(0) )
    CALL Info('WriteMeshToDiskPartitioned','Writing parallel mesh to disk: '//DirectoryName)
   

    DO Partition = 1, NoPartitions 
      
      CALL Info('WriteMeshToDiskPartitioned','Writing piece to file: '//I2S(Partition),Level=12)
      
      PrefixName = DirectoryName//'/part.'//I2S(Partition)

      CALL Info('WriteMeshToDiskPartitioned','Write nodes file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.nodes', STATUS='UNKNOWN' )
      NoNodes = 0
      DO i=1,Mesh % NumberOfNodes
        IF( ANY( NeighbourList(i) % Neighbours == Partition ) ) THEN
          WRITE(1,'(I0,x,I0,x,3ES17.10)') i,-1, &
              Mesh % Nodes % x(i), Mesh % Nodes % y(i), Mesh % Nodes % z(i)
          NoNodes = NoNodes + 1
        END IF
      END DO
      CLOSE(1)
      

      CALL Info('WriteMeshToDiskPartitioned','Write shared nodes file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.shared', STATUS='UNKNOWN' )
      NoShared = 0
      DO i=1,Mesh % NumberOfNodes
        nneigh = SIZE( NeighbourList(i) % Neighbours )
        IF( nneigh <= 1 ) CYCLE
        
        IF( ANY( NeighbourList(i) % Neighbours == Partition ) ) THEN
          NoShared = NoShared + 1
          WRITE(1,'(i0, x, i0, x)',ADVANCE='NO') i,nneigh
          DO j=1,nneigh
            WRITE(1,'(I0, x)',ADVANCE='NO') NeighbourList(i) % Neighbours(j) 
          END DO
          WRITE( 1,* ) ''
        END IF
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write elements file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.elements', STATUS='UNKNOWN' )
      NoBulkElements = 0
      ElmCodeCounts = 0      
      DO i=1,Mesh % NumberOfBulkElements
        IF( ElementPart(i) /= Partition ) CYCLE

        Element => Mesh % Elements(i)
        WRITE(1,'(i0,x,i0,x,i0,x)',ADVANCE='NO') i, &
            Element % BodyId, Element % TYPE % ElementCode
        DO j=1,Element % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') Element % NodeIndexes(j)
        END DO
        WRITE(1,*) ''
        
        ElmCode = Element % TYPE % ElementCode
        ElmCodeCounts( ElmCode ) = ElmCodeCounts( ElmCode ) + 1
        NoBulkElements = NoBulkElements + 1
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write boundary file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.boundary', STATUS='UNKNOWN' )
      NoBoundaryElements = 0
      DO i=Mesh % NumberOfBulkElements +1 ,&
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(i)
       
        parent1 = 0
        parent2 = 0
        Constraint = 0
        
        IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
          IF ( ASSOCIATED( Element % BoundaryInfo % Left ) ) &
              parent1 = Element % BoundaryInfo % Left % ElementIndex
          IF ( ASSOCIATED( Element % BoundaryInfo % Right ) ) &
              parent2 = Element % BoundaryInfo % Right % ElementIndex        
          Constraint = Element % BoundaryInfo % Constraint
        END IF

        Hit = .FALSE.
        IF( parent1 > 0 ) THEN
          IF( ElementPart( parent1 ) == Partition ) Hit = .TRUE.
        END IF
        IF( parent2 > 0 ) THEN
          IF( ElementPart( parent2 ) == Partition ) Hit = .TRUE.
        END IF

        IF( .NOT. Hit ) CYCLE

        WRITE(1,'(i0,x,i0,x,i0,x,i0,x,i0)',ADVANCE='NO') i, & 
            Constraint, Parent1, Parent2,&
            Element % TYPE % ElementCode
        DO j=1,Element % TYPE % NumberOfNodes
          WRITE(1,'(x,i0)', ADVANCE='NO') Element % NodeIndexes(j)
        END DO
        WRITE(1,*) 

        ElmCode = Element % TYPE % ElementCode
        ElmCodeCounts( ElmCode ) = ElmCodeCounts( ElmCode ) + 1
        NoBoundaryElements = NoBoundaryElements + 1
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write header file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.header',STATUS='UNKNOWN' )
      NumElmCodes = COUNT( ElmCodeCounts > 0 ) 
      WRITE( 1,'(i0,x,i0,x,i0)' ) NoNodes, &
          NoBulkElements, NoBoundaryElements      
      WRITE( 1,'(i0)' ) NumElmCodes
      DO i=SIZE(ElmCodeCounts),1,-1
        IF( ElmCodeCounts(i) == 0 ) CYCLE
        WRITE( 1,'(i0,x,i0,x)' ) i,ElmCodeCounts(i)
      END DO
      WRITE( 1,'(i0,x,i0)') NoShared, 0
      CLOSE(1)
      
      CALL Info('WriteMeshToDiskPartitioned','Done writing partition',Level=12)
    END DO

    CALL Info('WriteMeshToDiskPartitioned','Done writing parallel mesh',Level=8)

!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDiskPartitioned
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Show mesh size information
!------------------------------------------------------------------------------
  SUBROUTINE PrintMeshSize( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: na, nb, nn, ne, nf, no, ns, i
    INTEGER :: napar(0:2), nbpar(0:2), nnpar(0:2), nepar(0:2), nfpar(0:2), nopar(0:2), nspar(0:2)
    CHARACTER(*), PARAMETER :: Caller="PrintMeshSize"   
!------------------------------------------------------------------------------

    na = Mesh % NumberOfBulkElements
    nb = Mesh % NumberOfBoundaryElements
    nn = Mesh % NumberOfNodes
    ne = Mesh % NumberOfEdges
    nf = Mesh % NumberOfFaces

    IF( ParEnv % PEs > 1 .AND. .NOT. Mesh % SingleMesh ) THEN
      no = 0; ns = 0
      DO i=1,nn
        IF(Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPe) no = no+1
        IF(SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours) > 1) ns = ns+1
      END DO
      DO i=0,2
        napar(i) = ParallelReduction(na,i)
        nbpar(i) = ParallelReduction(nb,i)
        nnpar(i) = ParallelReduction(nn,i)
        nopar(i) = ParallelReduction(no,i)
        nspar(i) = ParallelReduction(ns,i)
        nepar(i) = ParallelReduction(ne,i)
        nfpar(i) = ParallelReduction(nf,i)
      END DO

      CALL Info(Caller,'Number of parallel mesh entities:   SUM       MIN       MAX')
      WRITE(Message,'(A,T30,3I10)') '  Bulk elements: ',napar
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,3I10)') '  Boundary elements: ',nbpar
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,3I10)') '  Total nodes: ',nnpar
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,3I10)') '  Owned nodes: ',nopar
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,3I10)') '  Shared nodes: ',nspar
      CALL Info(Caller,Message,Level=3)
      IF(nepar(0) > 0) THEN
        WRITE(Message,'(A,T30,3I10)') '  Element edges: ',nepar
        CALL Info(Caller,Message,Level=3)
      END IF
      IF(nfpar(0) > 0) THEN
        WRITE(Message,'(A,T30,3I10)') '  Element faces: ',nfpar
        CALL Info(Caller,Message,Level=3)
      END IF
    ELSE
      CALL Info(Caller,'Number of serial mesh entities')
      WRITE(Message,'(A,T30,1I10)') '  Bulk elements: ',na
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,1I10)') '  Boundary elements: ',nb
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,1I10)') '  Element nodes: ',nn
      CALL Info(Caller,Message,Level=3)
      IF(ne > 0) THEN
        WRITE(Message,'(A,T30,1I10)') '  Element edges: ',ne
        CALL Info(Caller,Message,Level=3)
      END IF
      IF(nf > 0) THEN
        WRITE(Message,'(A,T30,1I10)') '  Element faces: ',nf
        CALL Info(Caller,Message,Level=3)
      END IF
    END IF

  END SUBROUTINE PrintMeshSize

      
  
!------------------------------------------------------------------------------
!> Check mesh for various info. Mainly for debugging.
!------------------------------------------------------------------------------
  SUBROUTINE CheckMeshInfo( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: na, nb, nn
    INTEGER :: i,j,k,t,ii,jj,maxi,mini
    INTEGER, ALLOCATABLE :: NodeHits(:), TypeHits(:)
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: mins, maxs, s, s2
    INTEGER(KIND=8) :: Dbg(10)
    LOGICAL :: Halt
    CHARACTER(*), PARAMETER :: Caller="CheckMeshInfo"   
!------------------------------------------------------------------------------

    CALL Info(Caller,'Checking mesh information')

    na = Mesh % NumberOfBulkElements
    nb = Mesh % NumberOfBoundaryElements
    nn = Mesh % NumberOfNodes
    Halt = .FALSE.
    
    ALLOCATE(TypeHits(827))
    ALLOCATE(NodeHits(nn))
    
    CALL CheckMeshBulkHits()
    CALL CheckMeshBoundaryHits()
    CALL CheckBCTags()
    CALL CheckParentIndeces()
    CALL CheckMeshGeomSize()
    CALL CheckMeshSerendipity()
    CALL CheckMeshBodyRadius()
    CALL CheckMeshFaces()
    CALL CheckMeshEdges()
    CALL CheckParallelInfo()
    CALL CheckParallelEdgeInfo()
    CALL CheckParallelFaceInfo()

    nn = ParallelReduction(nn)
    
    CALL Info(Caller,'Finished checking mesh!')

    IF(Halt) CALL Fatal(Caller,'Some checksum was invalid, cannot continue!')


  CONTAINS

    
    SUBROUTINE CheckMeshBulkHits()
      TypeHits = 0
      NodeHits = 0
      Dbg = 0

      DO t=1,na
        Element => Mesh % Elements(t)
        IF(.NOT. ASSOCIATED( Element % TYPE ) ) THEN
          CALL Fatal(Caller,'Element type not associated for bulk elem: '//I2S(t))
        END IF
        i = Element % TYPE % ElementCode        
        TypeHits(i) = TypeHits(i)+1
        IF(ANY(Element % NodeIndexes < 1 ) ) THEN
          PRINT *,'NodeIndexes:', Element % NodeIndexes 
          CALL Fatal(Caller,'Bulk element '//I2S(t)//' has non-positive index!')
        END IF
        IF(ANY(Element % NodeIndexes > nn ) ) THEN
          PRINT *,'NodeIndexes:', Element % NodeIndexes, ' vs. ', nn 
          CALL Fatal(Caller,'Bulk element '//I2S(t)//' has too large index!')
        END IF
        IF(ANY(Element % NodeIndexes <= 0)) THEN
          PRINT *,'NodeIndexes:',Element % NodeIndexes
          CALL Fatal(Caller,'Non-positive node index encountered')
        END IF
        IF(ANY(Element % NodeIndexes < 1) ) THEN
          PRINT *,'Too small bulk element index: ',t, Element % NodeIndexes
        ELSE IF(ANY(Element % NodeIndexes > nn) ) THEN
          PRINT *,'Too large bulk element index: ',t, Element % NodeIndexes
        ELSE
          NodeHits(Element % NodeIndexes) = NodeHits(Element % NodeIndexes) + 1
        END IF
      END DO

      Dbg(1) = na
      Dbg(2) = SUM(NodeHits) 
      DO i=1,SIZE(NodeHits)
        Dbg(3) = dbg(3) + i*NodeHits(i)
      END DO

      DO i=1,SIZE(TypeHits)
        j = TypeHits(i)
        IF(j>0) CALL Info(Caller,'Bulk element type '//I2S(i)//' count: '//I2S(j))        
      END DO
      
      t=MAXVAL(NodeHits)

      IF( InfoActive(25)) THEN
        DO i=0,t
          j = COUNT( NodeHits == i)
          IF(j>0) PRINT *,'Bulk node hits '//I2S(i)//' count: ',j
        END DO
      END IF
      dbg(4) = t      
      Dbg(5) = COUNT(TypeHits>0)


      WRITE(Message,*) 'Bulk Checksum: ',Dbg(1:5)
      CALL Info(Caller,Message)      
      
      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckMeshBulkHits


    SUBROUTINE CheckMeshBoundaryHits()
      TypeHits = 0
      NodeHits = 0
      dbg = 0

      DO t=na+1,na+nb
        Element => Mesh % Elements(t)
        IF(.NOT. ASSOCIATED( Element % TYPE ) ) THEN
          CALL Fatal(Caller,'Element type not associated for bc elem: '//I2S(t-na))
        END IF
        i = Element % TYPE % ElementCode
        TypeHits(i) = TypeHits(i)+1
        IF(ANY(Element % NodeIndexes < 1 ) ) THEN
          PRINT *,'NodeIndexes:', Element % NodeIndexes 
          CALL Fatal(Caller,'Boundary element '//I2S(t)//' has non-positive index!')
        END IF
        IF(ANY(Element % NodeIndexes > nn ) ) THEN
          PRINT *,'NodeIndexes:', Element % NodeIndexes, ' vs. ', nn 
          CALL Fatal(Caller,'Boundary element '//I2S(t)//' has too large index!')
        END IF
        NodeHits(Element % NodeIndexes) = NodeHits(Element % NodeIndexes) + 1

        IF(ASSOCIATED(Element % BoundaryInfo) ) THEN
          IF(.NOT. ( ASSOCIATED(Element % BoundaryInfo % Left) .OR. &
              ASSOCIATED(Element % BoundaryInfo % Right) ) ) THEN
            PRINT *,'Boundary Info present but no left/right parent: ',t
          END IF
        END IF
      END DO
      DO i=1,SIZE(TypeHits)
        j = TypeHits(i)
        IF(j>0) CALL Info(Caller,'Boundary element type '//I2S(i)//' count: '//I2S(j))        
      END DO

      t=MAXVAL(NodeHits)
      IF(InfoActive(25)) THEN
        DO i=0,t
          j = COUNT( NodeHits == i)
          IF(j>0) PRINT *,'Boundary node hits '//I2S(i)//' count: ',j
        END DO
      END IF
        
      Dbg(1) = nb
      Dbg(2) = SUM(NodeHits) 
      DO i=1,SIZE(NodeHits)
        Dbg(3) = dbg(3) + i*NodeHits(i)
      END DO
      Dbg(4) = COUNT(TypeHits>0)
      Dbg(5) = t
      WRITE(Message,*) 'Boundary Checksum: ',Dbg(1:5)
      CALL Info(Caller,Message)      
     
      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckMeshBoundaryHits


    SUBROUTINE CheckParentIndeces()
      INTEGER :: Misses
      TYPE(Element_t), POINTER :: Parent

      Misses = 0
      dbg = 0
      
      DO t=na+1,na+nb
        Element => Mesh % Elements(t)
        i = Element % TYPE % NumberOfNodes
        IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE
        DO j=1,2
          IF(j==1) THEN
            Parent => Element % BoundaryInfo % Left
          ELSE
            Parent => Element % BoundaryInfo % Right
          END IF
          IF(.NOT. ASSOCIATED(Parent)) CYCLE

          dbg(3) = dbg(3) + 1
          dbg(4) = dbg(4) + Element % ElementIndex
          dbg(5) = dbg(5) + Element % BoundaryInfo % Constraint

          k = 0
          DO i=1,Element % TYPE % NumberOfNodes
            IF( .NOT. ANY(Parent % NodeIndexes == Element % NodeIndexes(i) ) ) k=k+1
          END DO
          IF( k > 0 ) THEN
            Misses = Misses + 1
            IF( Misses <= 10 ) THEN
              PRINT *,'Indeces missing in parent: ',ParEnv % Mype, Element % ElementIndex,Parent % ElementIndex, &
                  Element % TYPE % NumberOfNodes, k
              PRINT *,'Element codes:',Element % TYPE % ElementCode, &
                  Parent % TYPE % elementCode
              PRINT *,'bc elem inds:',Element % NodeIndexes
              PRINT *,'bulk elem inds:',Parent % NodeIndexes 
            END IF
          END IF
        END DO
      END DO

      IF(Misses>0) PRINT *,'Parent elements missing nodes:',ParEnv % Mype, Misses      
      dbg(1) = nb
      dbg(2) = Misses
      
      WRITE(Message,*) 'Parent Checksum: ',Dbg(1:5)
      CALL Info(Caller,Message)

      IF(Misses > 0) CALL Fatal(Caller,'We need all parent indeces!')

      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckParentIndeces


    SUBROUTINE CheckBCTags()
      INTEGER :: Misses, Tag, MinTag, MaxTag
      INTEGER, ALLOCATABLE :: TagCount(:), BCNodeCount(:)

      ! Not BCs to go through
      IF(nb==0) RETURN

      Misses = 0
      MinTag = HUGE(MinTag)
      MaxTag = -HUGE(MaxTag)

      DO k=1,2
        DO t=na+1,na+nb
          Element => Mesh % Elements(t)
          IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) THEN
            IF(k==1) Misses = Misses + 1
            CYCLE
          END IF
          Tag = Element % BoundaryInfo % Constraint
          IF(k==1) THEN
            MinTag = MIN(MinTag,Tag)
            MaxTag = MAX(MaxTag,Tag)
          ELSE
            TagCount(Tag) = TagCount(Tag) + 1
          END IF
        END DO
        IF(k==1) THEN
          ! Not tags defined in this partition
          IF(MinTag > MaxTag) THEN
            EXIT
          ELSE
            ALLOCATE(TagCount(MinTag:MaxTag))
            TagCount = 0
          END IF
        END IF
      END DO

      ALLOCATE(BCNodeCount(MinTag:MaxTag))
      BCNodeCount = 0
      
      DO k=1, MaxTag
        IF(TagCount(k)==0) CYCLE
        NodeHits = 0
        DO t=na+1,na+nb
          Element => Mesh % Elements(t)
          IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE
          Tag = Element % BoundaryInfo % Constraint
          IF(Tag==k) NodeHits(Element % NodeIndexes) = 1
        END DO
        BCNodeCount(k) = COUNT(NodeHits == 1)
      END DO

      DO k=MinTag,MaxTag
        IF(TagCount(k) > 0) THEN
          PRINT *,'BC'//I2S(k)//': elems '//I2S(TagCount(k))//' nodes '//I2S(BCNodeCount(k))
        END IF
      END DO           
      
    END SUBROUTINE CheckBCTags

    
    
    SUBROUTINE CheckMeshGeomSize()

      IF(.NOT. InfoActive(25)) RETURN
      
      PRINT *,'Coordinate x: ',MINVAL(Mesh % Nodes % x), MAXVAL(Mesh % Nodes % x)
      PRINT *,'Coordinate y: ',MINVAL(Mesh % Nodes % y), MAXVAL(Mesh % Nodes % y)
      PRINT *,'Coordinate z: ',MINVAL(Mesh % Nodes % z), MAXVAL(Mesh % Nodes % z)

      mins = HUGE(mins); maxs = 0.0_dp
      DO t=1,na+nb
        Element => Mesh % Elements(t)
        DO i=1,Element % TYPE % NumberOfNodes
          ii = Element % NodeIndexes(i)
          DO j=i+1, Element % TYPE % NumberOfNodes
            jj = Element % NodeIndexes(j)
            s2 = (Mesh % Nodes % x(ii)-Mesh % Nodes % x(jj))**2 + &
                (Mesh % Nodes % y(ii)-Mesh % Nodes % y(jj))**2 + &
                (Mesh % Nodes % z(ii)-Mesh % Nodes % z(jj))**2
            IF( s2 < mins ) THEN
              mins = s2
              mini = t 
            END IF
            IF( s2 > maxs ) THEN
              maxs = s2
              maxi = t
            END IF
          END DO
        END DO

        IF( t==na .OR. t==na+nb) THEN
          mins = SQRT(mins)
          maxs = SQRT(maxs)            
          IF(t==na) THEN
            PRINT *,'Bulk element h range:',mins,maxs          
            mins = HUGE(mins); maxs = 0.0_dp
          ELSE
            PRINT *,'Boundary element h range:',mins,maxs
          END IF

          Element => Mesh % Elements(maxi)
          PRINT *,'Maximum element:',maxi
          PRINT *,'x:',Mesh % Nodes % x(Element % NodeIndexes)
          PRINT *,'y:',Mesh % Nodes % y(Element % NodeIndexes)
          PRINT *,'z:',Mesh % Nodes % z(Element % NodeIndexes)

        END IF
      END DO
      
    END SUBROUTINE CheckMeshGeomSize

    SUBROUTINE CheckMeshSerendipity()
      INTEGER :: ElemCode 
      INTEGER :: Indexes0(27),EdgeInds(2),n,ne
      INTEGER, POINTER :: Indexes(:)
      REAL(KIND=dp) :: Coord(3),Coord0(3)
      
      DO t=1,na
        Element => Mesh % Elements(t)

        n = Element % Type % NumberOfNodes
        ne = Element % Type % NumberOfEdges
        
        ElemCode = Element % TYPE % ElementCode
        Indexes => Element % NodeIndexes
        Indexes0(1:n) = Indexes(1:n)
        
        SELECT CASE( ElemCode )
        CASE( 306, 408 )

          DO i=1,ne
            EdgeInds(1) = Indexes(i)
            IF(i==ne) THEN
              EdgeInds(2) = Indexes(1)
            ELSE
              EdgeInds(2) = Indexes(i+1)
            END IF

            ! Center of edge 
            Coord0(1) = SUM( Mesh % Nodes % x(EdgeInds)) / 2
            Coord0(2) = SUM( Mesh % Nodes % y(EdgeInds)) / 2
            Coord0(3) = SUM( Mesh % Nodes % z(EdgeInds)) / 2

            ! Is there some node closer to center of edge?
            maxs = HUGE(maxs)
            DO j=ne+1,n
              Coord(1) = Mesh % Nodes % x(Indexes(j)) 
              Coord(2) = Mesh % Nodes % y(Indexes(j)) 
              Coord(3) = Mesh % Nodes % z(Indexes(j)) 
              s2 = SUM((Coord-Coord0)**2)
              IF(s2 < maxs ) THEN
                Indexes0(ne+i) = Indexes(j)
                maxs = s2
              END IF              
            END DO
          END DO

        END SELECT
          
        j = COUNT( Indexes(1:n) /= Indexes0(1:n) )
        IF( j > 0 ) THEN
          !PRINT *,'Discrepancy: ',Indexes(ne+1:n), Indexes0(ne+1:n)
          Element % NodeIndexes(1:n) = Indexes0(1:n)
          CALL Warn('CheckMeshInfo','Node order wrong for '//I2S(j)//' nodes in element '//I2S(t))
        END IF
          
      END DO
    END SUBROUTINE CheckMeshSerendipity


    SUBROUTINE CheckMeshBodyRadius()
      REAL(KIND=dp) :: r
      REAL(KIND=dp), ALLOCATABLE :: RadRange(:,:)
      INTEGER, POINTER :: Indexes(:)
      INTEGER, ALLOCATABLE :: BodyHits(:)
      INTEGER :: nob

      IF(.NOT. InfoActive(25)) RETURN
      
      nob = CurrentModel % NumberOfBodies
      ALLOCATE(RadRange(0:nob,2),BodyHits(0:nob))
      RadRange(:,1) = HUGE(r)
      RadRange(:,2) = 0.0_dp
      BodyHits = 0
      
      DO t=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        Indexes => Element % NodeIndexes
        r = 0.0_dp
        k = Element % BodyId
        IF(k<0 .OR. k>nob) CYCLE
        DO i=1, Element % TYPE % NumberOfNodes
          j = Indexes(i)
          r = Mesh % Nodes % x(j)**2
          r = r + Mesh % Nodes % y(j)**2
          r = r + Mesh % Nodes % z(j)**2
        END DO
        RadRange(k,1) = MIN(RadRange(k,1),r)
        RadRange(k,2) = MAX(RadRange(k,2),r)                   
        BodyHits(k) = BodyHits(k)+1
      END DO
      RadRange = SQRT( RadRange )
      DO i=0,nob
        IF(BodyHits(i)==0) CYCLE
        PRINT *,'Radius range: ',i,RadRange(i,:)
      END DO
    END SUBROUTINE CheckMeshBodyRadius


    SUBROUTINE CheckMeshEdges()
      INTEGER, POINTER :: Indexes(:)
      INTEGER :: m

      IF(Mesh % NumberOfEdges == 0 ) RETURN      
      dbg = 0
      dbg(1) = Mesh % NumberOfEdges

      DO t=1,Mesh % NumberOfEdges
        Element => Mesh % Edges(t)
        IF(.NOT. ASSOCIATED(Element)) THEN
          CALL Fatal(Caller,'Edge not associated on edge list: '//I2S(t))
        END IF
        Indexes => Element % NodeIndexes          
        IF(.NOT. ASSOCIATED(Indexes)) THEN
          CALL Fatal(Caller,'NodeIndexes not associated on edge: '//I2S(t))
        END IF
        IF(.NOT. ASSOCIATED(Element % TYPE)) THEN
          CALL Fatal(Caller,'Edge type '//I2S(t)//' not associated!')
        END IF
        m = Element % Type % NumberOfNodes
        IF(SIZE(Indexes) /= m) THEN
          CALL Fatal(Caller,'Invalid size of edge '//I2S(t)//&
              ' NodeIndexes: '//I2S(SIZE(Indexes))//' vs. '//I2S(m))
        END IF
        IF(SIZE(Indexes)>0) dbg(2) = dbg(2) + SUM(Indexes)
        dbg(3) = dbg(3) + Element % ElementIndex
        dbg(4) = dbg(4) + Element % GElementIndex        
      END DO

      WRITE(Message,*) 'Edges Checksum: ',Dbg(1:5)
      CALL Info(Caller,Message)                  

      !IF(ANY(Dbg < 0) ) Halt = .TRUE.

    END SUBROUTINE CheckMeshEdges

    SUBROUTINE CheckMeshFaces()
      INTEGER, POINTER :: Indexes(:)
      INTEGER :: m

      IF(Mesh % NumberOfFaces == 0 ) RETURN      
      dbg = 0
      dbg(1) = Mesh % NumberOfFaces

      DO t=1,Mesh % NumberOfFaces
        Element => Mesh % Faces(t)
        IF(.NOT. ASSOCIATED(Element)) THEN
          CALL Fatal(Caller,'Face not associated on face list: '//I2S(t))
        END IF
        Indexes => Element % NodeIndexes          
        IF(.NOT. ASSOCIATED(Indexes)) THEN
          CALL Fatal(Caller,'NodeIndexes not associated on face: '//I2S(t))
        END IF
        IF(.NOT. ASSOCIATED(Element % TYPE)) THEN
          CALL Fatal(Caller,'Face type '//I2S(t)//' not associated!')
        END IF
        m = Element % TYPE % NumberOfNodes
        IF(SIZE(Indexes) /= m) THEN
          CALL Fatal(Caller,'Invalid size of face '//I2S(t)//&
              ' NodeIndexes: '//I2S(SIZE(Indexes))//' vs. '//I2S(m))
        END IF
        IF(SIZE(Indexes)>0) dbg(2) = dbg(2) + SUM(Indexes)
        dbg(3) = dbg(3) + Element % ElementIndex
        dbg(4) = dbg(4) + Element % GElementIndex        

      END DO

      WRITE(Message,*) 'Faces Checksum: ',Dbg(1:5)
      CALL Info(Caller,Message)                  

      !IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckMeshFaces


    SUBROUTINE CheckParallelInfo()

      IF( ParEnv % PEs == 1) RETURN
      IF(.NOT. ASSOCIATED( Mesh % ParallelInfo % NeighbourList) ) RETURN
      
      dbg = 0
      dbg(1) = SIZE(Mesh % ParallelInfo % NeighbourList)

      dbg(2) = COUNT(Mesh % ParallelInfo % Ginterface)
      DO i=1, SIZE(Mesh % ParallelInfo % Ginterface)        
        IF( Mesh % ParallelInfo % Ginterface(i) ) dbg(3) = dbg(3) + i 
      END DO
      
      DO i=1, SIZE(Mesh % ParallelInfo % NeighbourList)
        IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)) THEN
          dbg(7) = dbg(7) + 1
          CYCLE
        END IF
        j = SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
        dbg(4) = dbg(4) + j
        dbg(5) = dbg(5) + SUM(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
        dbg(6) = dbg(6) + i*j
      END DO

      WRITE(Message,*) 'ParallelInfo Checksum: ',Dbg(1:7)
      CALL Info(Caller,Message)                         

      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
     END SUBROUTINE CheckParallelInfo
       
    SUBROUTINE CheckParallelEdgeInfo()

      IF( ParEnv % PEs == 1) RETURN
      IF( Mesh % NumberOfEdges == 0) RETURN
      IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % EdgeNeighbourList)) RETURN
      
      dbg = 0      
      dbg(1) = SIZE(Mesh % ParallelInfo % EdgeNeighbourList)

      IF(ASSOCIATED(Mesh % ParallelInfo % EdgeInterface ) ) THEN
        j = SIZE(Mesh % ParallelInfo % EdgeInterface )        
        IF(j>1) THEN
          dbg(2) = j
          DO i=1, SIZE(Mesh % ParallelInfo % Edgeinterface)        
            IF( Mesh % ParallelInfo % Edgeinterface(i) ) dbg(3) = dbg(3) + i 
          END DO
        END IF
      END IF
        
      DO i=1, SIZE(Mesh % ParallelInfo % EdgeNeighbourList)
        IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours)) THEN
          dbg(7) = dbg(7) + 1
          CYCLE
        END IF
        j = SIZE(Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours)
        dbg(4) = dbg(4) + j
        dbg(5) = dbg(5) + SUM(Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours)
        dbg(6) = dbg(6) + i*j
      END DO

      WRITE(Message,*) 'ParallelEdges Checksum: ',Dbg(1:7)
      CALL Info(Caller,Message)                         
      
      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckParallelEdgeInfo

    SUBROUTINE CheckParallelFaceInfo()

      IF( ParEnv % PEs == 1) RETURN
      IF( Mesh % NumberOfFaces == 0) RETURN
      IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % FaceNeighbourList)) RETURN
      
      dbg = 0
      dbg(1) = SIZE(Mesh % ParallelInfo % FaceNeighbourList)

      IF( ASSOCIATED( Mesh % ParallelInfo % FaceInterface ) ) THEN
        j = SIZE(Mesh % ParallelInfo % FaceInterface )
        IF(j>1) THEN
          dbg(2) = j
          DO i=1, j
            IF( Mesh % ParallelInfo % Faceinterface(i) ) dbg(3) = dbg(3) + i 
          END DO
        END IF
      END IF

      DO i=1, SIZE(Mesh % ParallelInfo % FaceNeighbourList)
        IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours)) THEN
          dbg(7) = dbg(7) + 1
          CYCLE
        END IF
        j = SIZE(Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours)
        dbg(4) = dbg(4) + j
        dbg(5) = dbg(5) + SUM(Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours)
        dbg(6) = dbg(6) + i*j
      END DO

      WRITE(Message,*) 'ParallelFaces Checksum: ',Dbg(1:7)
      CALL Info(Caller,Message)                         
      
      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckParallelFaceInfo
           
!------------------------------------------------------------------------------
  END SUBROUTINE CheckMeshInfo
!------------------------------------------------------------------------------


END MODULE MeshIO
