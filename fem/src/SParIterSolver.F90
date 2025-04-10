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
! *  Authors: Jouni Malinen, Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2000
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!> Module for the parallel iterative solution in Elmer.
!-------------------------------------------------------------------------------


#include "huti_fdefs.h"
#include "../config.h"

MODULE SParIterSolve
  USE Types
  USE Lists
  USE ListMatrix
  USE SParIterGlobals
  USE SParIterComm
  USE SParIterPrecond
  USE IterSolve, ONLY : NumericalError

  USE CRSMatrix

  IMPLICIT NONE

CONTAINS

  
  SUBROUTINE Dummy
  END SUBROUTINE Dummy


  SUBROUTINE DPCond( u,v,ipar )
    REAL(KIND=dp) :: u(*),v(*)
    INTEGER :: ipar(*)

    u(1:HUTI_NDIM) = v(1:HUTI_NDIM)
  END SUBROUTINE DPCond


  SUBROUTINE CPCond( u,v,ipar )
    COMPLEX(KIND=dp) :: u(*),v(*)
    INTEGER :: ipar(*)

    u(1:HUTI_NDIM) = v(1:HUTI_NDIM)
  END SUBROUTINE CPCond

  !********************************************************************
  !
  ! Initialize the Matrix structures for parallel environment
  !
  FUNCTION ParInitMatrix( SourceMatrix, ParallelInfo ) RESULT ( SParMatrixDesc )
    TYPE (Matrix_t),TARGET :: SourceMatrix
    TYPE (ParallelInfo_t), TARGET :: ParallelInfo
    TYPE (SParIterSolverGlobalD_t), POINTER :: SParMatrixDesc

    TYPE (ParEnv_t), POINTER :: ParallelEnv
    !******************************************************************

    ALLOCATE( SParMatrixDesc )

    SParMatrixDesc % ParEnv = ParEnv
    ALLOCATE(SParMatrixDesc % ParEnv % Active(ParEnv % PEs))
    SParMatrixDesc % ParEnv % Active = ParEnv % Active
    SParMatrixDesc % ParEnv % IsNeighbour => Null()
    ParEnv => SParMatrixDesc % ParEnv

    CALL ParEnvInit(SParMatrixDesc, ParallelInfo, SourceMatrix)

    SParMatrixDesc % Matrix => SourceMatrix
    SParMatrixDesc % DOFs = 1
    SParMatrixDesc % ParallelInfo => ParallelInfo

    ParEnv % ActiveComm = SourceMatrix % Comm

    SParMatrixDesc % SplittedMatrix => &
                        SplitMatrix( SourceMatrix, ParallelInfo )

  END FUNCTION ParInitMatrix


!********************************************************************
!********************************************************************
!
! Split the given matrix (SourceMatrix) to InsideMatrix part (elements
! inside a partition) and InterfaceMatrix table (interface blocks).
!
  FUNCTION SplitMatrix( SourceMatrix, ParallelInfo ) &
       RESULT ( SplittedMatrix )

    USE Types
    IMPLICIT NONE

    TYPE (Matrix_t) :: SourceMatrix       ! Original matrix in this partition
    TYPE (ParallelInfo_t) :: ParallelInfo
    TYPE (SplittedMatrixT), POINTER :: SplittedMatrix

  ! External routines

! EXTERNAL SearchNode
! INTEGER :: SearchNode

  ! Local variables

    INTEGER :: i, j, k, l, m, n, sz, dd, RowInd, ColInd, RowStart, RowEnd, Gcol
    INTEGER :: InsideMRows, InsideMCols, Row, Col, currifi, RowOwner
    INTEGER, DIMENSION(:), ALLOCATABLE :: OwnIfMRows, OwnIfMCols, OwnOCOR
    INTEGER, DIMENSION(:), ALLOCATABLE :: NbsIfMRows, NbsIfMCols, NbsOCOR
    INTEGER, DIMENSION(:), ALLOCATABLE :: OwnOldCols, NbsOldCols, Perm(:)
    INTEGER :: FoundIndex, AtLeastOneCol, OldCols, NewCol, iii
    INTEGER, POINTER :: ColBuffer(:), RowBuffer(:), IndBuffer(:), Ncol(:)

    TYPE (NeighbourList_t), POINTER :: CurrNbsL

    TYPE(Matrix_t), POINTER :: A
    TYPE (BasicMatrix_t), DIMENSION(:), ALLOCATABLE :: OwnIfMatrix
    TYPE (BasicMatrix_t), POINTER :: NbsIfMatrix(:), CurrIf
    TYPE (BasicMatrix_t), DIMENSION(:), ALLOCATABLE :: RecvdIfMatrix

    LOGICAL, ALLOCATABLE :: isNeighbour(:)
    LOGICAL :: NeedMass, NeedDamp, NeedPrec, NeedILU, GotNewCol, Found
    CHARACTER(:), ALLOCATABLE :: Prec
    REAL(kind=dp) :: st
  !******************************************************************
    st = realtime()
    ALLOCATE( SplittedMatrix )
    SplittedMatrix % InsideMatrix => AllocateMatrix()

    ALLOCATE( SplittedMatrix % GlueTable )

    NULLIFY( SplittedMatrix % IfMatrix )
    NULLIFY( SplittedMatrix % NbsIfMatrix )
    NULLIFY( SplittedMatrix % Vecindices )
    NULLIFY( SplittedMatrix % IfVecs )
    NULLIFY( SplittedMatrix % RHS )
    NULLIFY( SplittedMatrix % Work )
    NULLIFY( SplittedMatrix % ResBuf )
    NULLIFY( SplittedMatrix % TmpXVec )
    NULLIFY( SplittedMatrix % TmpRVec )

  !----------------------------------------------------------------------
  !
  ! Copy the glueing table row and column indices
  !
  !----------------------------------------------------------------------

    ALLOCATE(SplittedMatrix % GlueTable % Inds(SIZE(SourceMatrix % Cols)))
    ALLOCATE( SplittedMatrix % VecIndices(ParEnv % PEs) )
    CALL CountNeighbourConns( SourceMatrix, SplittedMatrix, ParallelInfo )

#ifdef HAVE_HYPRE
! the blocksolver needs the parallel structures for matrix-vector prod.,
! even if using HYPRE
!   IF (ListGetLogical(CurrentModel % Solver % Values, &
!      'Linear System Use HYPRE', Found )) RETURN
#endif

  !----------------------------------------------------------------------
  !
  ! Allocate some temporary work space
  !
  !----------------------------------------------------------------------

    ALLOCATE( OwnIfMRows(ParEnv % PEs) )
    ALLOCATE( OwnIfMCols(ParEnv % PEs) )
    ALLOCATE( NbsIfMRows(ParEnv % PEs) )
    ALLOCATE( NbsIfMCols(ParEnv % PEs) )
    ALLOCATE( OwnOCOR(ParEnv % PEs) )
    ALLOCATE( NbsOCOR(ParEnv % PEs) )
    ALLOCATE( OwnOldCols(ParEnv % PEs) )
    ALLOCATE( NbsOldCols(ParEnv % PEs) )
    OwnIfMRows(:) = 0; OwnIfMCols(:) = 0; NbsIfMRows(:) = 0; NbsIfMCols(:) = 0

    j = 0
    k = SIZE(ParallelInfo % NeighbourList) 
    DO i=1,k
      IF(.NOT. ASSOCIATED( ParallelInfo % NeighbourList(i) % Neighbours ) ) THEN
        ALLOCATE(ParallelInfo % NeighbourList(i) % Neighbours(1))
        ParallelInfo % NeighbourList(i) % Neighbours(1) = ParEnv % MyPe        
        j = j+1
      END IF
    END DO
    IF( j > 0 ) THEN
      CALL Info('SplitMatrix','Added mention to self in '//I2S(j)//' neighbours ouf of '//TRIM(I2S(k)))
    END IF
    
  !----------------------------------------------------------------------
  !
  ! Compute the memory allocations for split matrix blocks
  !
  !----------------------------------------------------------------------
    InsideMRows = 0; InsideMCols = 0
    DO i = 1, SourceMatrix % NumberOfRows

       AtLeastOneCol = 0; OwnOCOR(:) = 0; NbsOCOR(:) = 0
       RowInd = i
       RowOwner = ParallelInfo % NeighbourList(RowInd) % Neighbours(1)

       IF ( RowOwner == ParEnv % MyPE ) THEN
!
!         Connection is local-local or local-if:
!         --------------------------------------
          DO j = SourceMatrix % Rows(i), SourceMatrix % Rows(i+1) - 1

             ColInd = SourceMatrix % Cols(j)                          
             IF ( ParallelInfo % NeighbourList(ColInd) % Neighbours(1) == &
                                ParEnv % MyPE ) THEN

              !----------------------------------------------------------
              !
              ! Connection is local-local
              !
              !----------------------------------------------------------

                InsideMCols = InsideMCols + 1
                AtLeastOneCol = 1
             ELSE

              !----------------------------------------------------------
              !
              ! Connection is local-if
              !
              !----------------------------------------------------------

                CurrNbsL => ParallelInfo % NeighbourList( ColInd )
                NbsIfMCols( CurrNbsL % Neighbours(1) + 1 ) = &
                     NbsIfMCols( CurrNbsL % Neighbours(1) + 1 ) + 1
                NbsOCOR( CurrNbsL % Neighbours(1) + 1 ) = 1
          END IF
        END DO

     ELSE

        !----------------------------------------------------------------
        !
        ! Connection is if-local or if-if
        !
        !----------------------------------------------------------------

        DO j = SourceMatrix % Rows(i), SourceMatrix % Rows(i+1) - 1

           ColInd = SourceMatrix % Cols(j)

           IF ( ParallelInfo % NeighbourList(ColInd) % Neighbours(1) == &
                            ParEnv % MyPE ) THEN

              !----------------------------------------------------------
              !
              ! Connection is if-local
              !
              !----------------------------------------------------------

              CurrNbsL => ParallelInfo % NeighbourList(RowInd)
              OwnIfMCols( CurrNbsL % Neighbours(1) + 1) = &
                   OwnIfMCols( CurrNbsL % Neighbours(1) + 1) + 1
              OwnOCOR( CurrNbsL % Neighbours(1) + 1) = 1
              
           ELSE

              !----------------------------------------------------------
              !
              ! Connection is if-if
              !
              !----------------------------------------------------------
              CurrNbsL => ParallelInfo % NeighbourList(ColInd)
              NbsIfMCols( CurrNbsL % Neighbours(1) + 1) = &
                   NbsIfMCols( CurrNbsL % Neighbours(1) + 1) + 1
              NbsOCOR( CurrNbsL % Neighbours(1) + 1) = 1
           END IF
        END DO

     END IF

     InsideMRows = InsideMRows + AtLeastOneCol
     NbsIfMRows(:) = NbsIfMRows(:) + NbsOCOR(:)
     OwnIfMRows(:) = OwnIfMRows(:) + OwnOCOR(:)
  END DO

  !----------------------------------------------------------------------
  !
  ! Allocate the block inside this partition
  !
  !----------------------------------------------------------------------
  ALLOCATE( SplittedMatrix % InsideMatrix % Rows( InsideMRows + 1 ) )
  ALLOCATE( SplittedMatrix % InsideMatrix % Cols( InsideMCols ) )
  ALLOCATE( SplittedMatrix % InsideMatrix % GRows( InsideMRows ) )
  ALLOCATE( SplittedMatrix % InsideMatrix % Gorder( InsideMRows ) )
  ALLOCATE( SplittedMatrix % InsideMatrix % Diag( InsideMRows ) )
  SplittedMatrix % InsideMatrix % Diag=0

  NULLIFY( SplittedMatrix % InsideMatrix % Values )
  NeedPrec = .FALSE.
  IF (  ASSOCIATED( SourceMatrix % PrecValues ) ) THEN
     IF ( SIZE(SourceMatrix % Values) == &
           SIZE(SourceMatrix % PrecValues) ) NeedPrec = .TRUE.
  END IF
  NULLIFY( SplittedMatrix % InsideMatrix % PrecValues )

  NeedMass = .FALSE.
  IF (  ASSOCIATED( SourceMatrix % MassValues ) ) THEN
     IF ( SIZE(SourceMatrix % Values) == &
           SIZE(SourceMatrix % MassValues) ) NeedMass = .TRUE.
  END IF
  NULLIFY( SplittedMatrix % InsideMatrix % MassValues )

  NeedDamp = .FALSE.
  IF (  ASSOCIATED( SourceMatrix % DampValues ) ) THEN
     IF ( SIZE(SourceMatrix % Values) == &
          SIZE( SourceMatrix % DampValues) ) NeedDamp = .TRUE.
  END IF
  NULLIFY( SplittedMatrix % InsideMatrix % DampValues )

  NeedIlu = .FALSE.
  IF( ListGetString(CurrentModel % Solver % Values, &
      'Linear System Solver', Found ) == 'iterative' ) THEN
    Prec = ListGetString(CurrentModel % Solver % Values, &
        'Linear System Preconditioning', Found ) 
    IF(Found) NeedILU = SEQL(Prec,'vanka')
  END IF
  
  SplittedMatrix % InsideMatrix % Ordered = .FALSE.
  NULLIFY( SplittedMatrix % InsideMatrix % ILUValues )

  !----------------------------------------------------------------------
  !
  ! Allocate the interface blocks (both the own parts and the ones to be
  ! sent to the neighbours)
  !
  !----------------------------------------------------------------------

  ALLOCATE( SplittedMatrix % IfMatrix(ParEnv % PEs) )
  ALLOCATE( SplittedMatrix % IfVecs(ParEnv % PEs) )
  ALLOCATE( SplittedMatrix % NbsIfMatrix(ParEnv % PEs) )
  ALLOCATE( OwnIfMatrix(ParEnv % PEs) )
  ALLOCATE( RecvdIfMatrix(ParEnv % PEs) )

  NbsIfMatrix => SplittedMatrix  % NbsIfMatrix
 
  SplittedMatrix % InsideMatrix % NumberOfRows = InsideMRows

  DO i = 1, ParEnv % PEs

     RecvdIfMatrix(i) % NumberOfRows = 0

     NbsIfMatrix(i) % NumberOfRows = NbsIfMRows(i)
     IF ( NbsIfMRows(i) /= 0 ) THEN
        ALLOCATE( NbsIfMatrix(i) % Rows(NbsIfMRows(i)+1) )
        ALLOCATE( NbsIfMatrix(i) % Cols(NbsIfMCols(i)) )
        ALLOCATE( NbsIfMatrix(i) % Diag(NbsIfMRows(i)) )
        ALLOCATE( NbsIfMatrix(i) % GRows(NbsIfMRows(i)) )
        ALLOCATE( NbsIfMatrix(i) % Values(NbsIfMCols(i)) )
        ALLOCATE( NbsIfMatrix(i) % RowOwner(NbsIfMRows(i)) )

        IF ( NeedPrec ) ALLOCATE( NbsIfMatrix(i) % PrecValues( NbsIfMCols(i) ) )

!       NULLIFY( NbsIfMatrix(i) % MassValues )
        IF ( NeedMass ) ALLOCATE( NbsIfMatrix(i) % MassValues( NbsIfMCols(i) ) )

!       NULLIFY( NbsIfMatrix(i) % DampValues )
        IF ( NeedDamp ) ALLOCATE( NbsIfMatrix(i) % DampValues( NbsIfMCols(i) ) )

!       NULLIFY( NbsIfMatrix(i) % ILUValues )
        IF ( NeedILU ) ALLOCATE( NbsIfMatrix(i) % ILUValues( NbsIfMCols(i) ) )
     END IF

     OwnIfMatrix(i) % NumberOfRows = OwnIfMRows(i)
     IF ( OwnIfMRows(i) /= 0 ) THEN
        ALLOCATE( OwnIfMatrix(i) % Rows(OwnIfMRows(i)+1) )
        ALLOCATE( OwnIfMatrix(i) % Cols(OwnIfMCols(i)) )
        ALLOCATE( OwnIfMatrix(i) % Diag(OwnIfMRows(i)) )
        ALLOCATE( OwnIfMatrix(i) % GRows(OwnIfMRows(i)) )
        ALLOCATE( OwnIfMatrix(i) % Values(OwnIfMCols(i)) )
        ALLOCATE( OwnIfMatrix(i) % RowOwner(OwnIfMRows(i)) )

        IF ( NeedPrec ) ALLOCATE( OwnIfMatrix(i) % PrecValues(OwnIfMCols(i)) )
        IF ( NeedMass ) ALLOCATE( OwnIfMatrix(i) % MassValues(OwnIfMCols(i)) )
        IF ( NeedDamp ) ALLOCATE( OwnIfMatrix(i) % DampValues(OwnIfMCols(i)) )
        IF ( NeedILU  ) ALLOCATE( OwnIfMatrix(i) % ILUValues(OwnIfMCols(i)) )
     END IF
  END DO

  !----------------------------------------------------------------------
  !
  ! Copy the actual indices into correct blocks
  !
  !----------------------------------------------------------------------

  InsideMRows = 1; InsideMCols = 1; OldCols = 1
  NbsIfMRows(:) = 1; NbsIfMCols(:) = 1; NbsOldCols(:) = 1
  OwnIfMRows(:) = 1; OwnIfMCols(:) = 1; OwnOldCols(:) = 1
  DO i = 1, SourceMatrix % NumberOfRows
     AtLeastOneCol = 0; NbsOCOR(:) = 0; OwnOCOR(:) = 0

     RowInd   = i
     RowOwner = ParallelInfo % NeighbourList(RowInd) % Neighbours(1)

     IF ( RowOwner == ParEnv % MyPE ) THEN
        !----------------------------------------------------------------
        !
        ! Connection is local-local or local-if
        !
        !----------------------------------------------------------------
        
        DO j = SourceMatrix % Rows(i), SourceMatrix % Rows(i+1) - 1

           ColInd = SourceMatrix % Cols(j)

           IF ( ParallelInfo % NeighbourList(ColInd) % Neighbours(1) == &
                ParEnv % MyPE ) THEN

              !----------------------------------------------------------
              !
              ! Connection is local-local
              !
              !----------------------------------------------------------

              SplittedMatrix % InsideMatrix % Cols( InsideMCols ) = &
                             SourceMatrix % Cols(j)
              SplittedMatrix % GlueTable % Inds(j) = InsideMCols
              AtLeastOneCol = 1
              InsideMCols = InsideMCols + 1
           ELSE

              !----------------------------------------------------------
              !
              ! Connection is local-if
              !
              !----------------------------------------------------------

              currifi = ParallelInfo % NeighbourList(Colind) % Neighbours(1) + 1

              NbsIfMatrix(currifi) % Cols(NbsIfMcols(currifi)) =  &
                  ParallelInfo  % GlobalDOFs(ColInd)
              NbsIfMcols(currifi) = NbsIfMcols(currifi) + 1
              NbsOCOR(currifi) = 1
              SplittedMatrix % GlueTable % Inds(j) = -(ParEnv % PEs + currifi)

           END IF
        END DO

     ELSE

        !----------------------------------------------------------------
        !
        ! Connection is if-local or if-if
        !
        !----------------------------------------------------------------

        DO j = SourceMatrix % Rows(i), SourceMatrix % Rows(i+1) - 1

           ColInd = SourceMatrix % Cols(j)

           IF ( ParallelInfo % NeighbourList(ColInd) % Neighbours(1) == &
              ParEnv % MyPE ) THEN

              !----------------------------------------------------------
              !
              ! Connection is if-local
              !
              !----------------------------------------------------------

              currifi = ParallelInfo % NeighbourList(RowInd) % Neighbours(1) + 1

              OwnIfMatrix(currifi) % Cols(OwnIfMcols(currifi)) =  &
                  ParallelInfo  % GlobalDOFs(ColInd)
              SplittedMatrix % GlueTable % Inds(j) = -currifi
              OwnOCOR(currifi) = 1
              OwnIfMcols(currifi) = OwnIfMcols(currifi) + 1

           ELSE

              !----------------------------------------------------------
              !
              ! Connection is if-if
              !
              !----------------------------------------------------------
              currifi = ParallelInfo % NeighbourList(ColInd) % Neighbours(1) + 1

              NbsIfMatrix(currifi) % Cols(NbsIfMcols(currifi)) =  &
                ParallelInfo  % GlobalDOFs(ColInd)
              SplittedMatrix % GlueTable % Inds(j) = -(ParEnv % PEs + currifi)
              NbsOCOR(currifi) = 1
              NbsIfMcols(currifi) = NbsIfMcols(currifi) + 1
           END IF
        END DO

     END IF

     !-------------------------------------------------------------------
     !
     ! Update the row indices to keep the CRS structures valid
     !
     !-------------------------------------------------------------------

     RowInd = ParallelInfo % GlobalDOFs(RowInd)

     IF ( AtLeastOneCol /= 0 ) THEN
        SplittedMatrix % InsideMatrix % GRows( InsideMRows ) = RowInd
        SplittedMatrix % InsideMatrix % Rows(  InsideMRows ) = OldCols
     END IF

     DO j = 1, ParEnv % PEs
        IF ( OwnOCOR(j) /= 0 ) THEN
           OwnIfMatrix(j) % Rows(OwnIfMRows(j))  = OwnOldCols(j)
           OwnIfMatrix(j) % GRows(OwnIfMRows(j)) = RowInd
           OwnIfMatrix(j) % RowOwner(OwnIfMRows(j)) = RowOwner
        END IF

        IF ( NbsOCOR(j) /= 0 ) THEN
           NbsIfMatrix(j) % Rows(NbsIfMRows(j))  = NbsOldCols(j)
           NbsIfMatrix(j) % GRows(NbsIfMRows(j)) = RowInd
           NbsIfMatrix(j) % RowOwner(NbsIfMRows(j)) = RowOwner
        END IF
     END DO

     InsideMRows = InsideMRows + AtLeastOneCol
     NbsIfMRows(:) = NbsIfMRows(:) + NbsOCOR(:)
     OwnIfMRows(:) = OwnIfMRows(:) + OwnOCOR(:)

     OldCols = InsideMCols
     OwnOldCols(:) = OwnIfMCols(:)
     NbsOldCols(:) = NbsIfMCols(:)
  END DO

  SplittedMatrix % InsideMatrix % Rows( InsideMRows ) = InsideMCols
  DO j = 1, ParEnv % PEs
     IF ( OwnIfMatrix(j) % NumberOfRows /= 0 ) &
          OwnIfMatrix(j) % Rows(OwnIfMRows(j)) = OwnIfMCols(j)

     IF ( NbsIfMatrix(j) % NumberOfRows /= 0 ) &
          NbsIfMatrix(j) % Rows(NbsIfMRows(j)) = NbsIfMCols(j)
  END DO

  DO j = 1, ParEnv % PEs
     IF ( OwnIfMatrix(j) % NumberOfRows /= (OwnIfMRows(j) - 1) .OR. &
          NbsIfMatrix(j) % NumberOfRows /= (NbsIfMRows(j) - 1) ) THEN
        WRITE( Message, * ) OwnIfMRows, NbsIfMRows, OwnIfMCols, NbsIfMCols
        CALL Error( 'SplitMatrix', Message )
     END IF
  END DO

  !----------------------------------------------------------------------
  !
  ! Exchange the interface blocks and glue received interface parts into
  ! the interface blocks already at this processor.
  !
  !----------------------------------------------------------------------
!if ( parenv % mype==0 ) print*,'SPLIT INIT TIME: ', realtime()-st; st = realtime()
  CALL ExchangeInterfaces( NbsIfMatrix, RecvdIfMatrix )
!if ( parenv % mype==0 ) print*, 'EXCHANGE INTERF TIME: ', realtime()-st; st=realtime()

  !----------------------------------------------------------------------
  !
  ! sort inside matrix global tags to speed lookup
  !
  !----------------------------------------------------------------------
  A=>SplittedMatrix % InsideMatrix
  n = A % NumberOFRows
  DO i=1,n
    A % Gorder(i) = i
  END DO

  CALL SortI( n, A % GRows, A % GOrder )

  !----------------------------------------------------------------------
  !
  ! allocate stuff
  !
  !----------------------------------------------------------------------
  sz = A % Rows(n+1)-1
  ALLOCATE(A % Values(sz))

!-----------------------------------------------------------------------------

  SplittedMatrix % IfMatrix(:) % NumberOfRows = 0
  DO i = 1, ParEnv % PEs
     IF ( ( OwnIfMatrix(i) % NumberOfRows   /= 0 ) .OR. &
          ( RecvdIfMatrix(i) % NumberOfRows /= 0 ) ) THEN
        CALL CombineCRSMatIndices( OwnIfMatrix(i), RecvdIfMatrix(i), &
                 SplittedMatrix % IfMatrix(i) )
     END IF
  END DO

  DO i=1,ParEnv % PEs
    IF ( i-1==Parenv % MyPE ) CYCLE
    Currif => SplittedMatrix % IfMatrix(i)
    DO j=1,Currif % NumberOfRows
      IF ( CurrIf % RowOwner(j)==Parenv % MyPE ) CYCLE
      IF( .NOT. ParEnv % IsNeighbour(Currif % RowOwner(j)+1) ) THEN
         ParEnv % NumOfNeighbours = ParEnv % NumOfNeighbours+1
         ParEnv % isNeighbour(CurrIf % RowOwner(j)+1) = .TRUE.
      END IF
    END DO
  END DO

  ! Sync neighbour information, if changed by the above ^:
  ! ------------------------------------------------------
! ALLOCATE(isNeighbour(ParEnv % PEs))
! isNeighbour = Parenv % IsNeighbour
! CALL MPI_ALLREDUCE( isNeighbour, Parenv % IsNeighbour, Parenv % Pes, &
!     MPI_LOGICAL, MPI_LOR, SourceMatrix % Comm, i )
! Parenv % IsNeighbour(Parenv % myPE+1) = .FALSE.
! Parenv % NumOfNeighbours = COUNT(Parenv % IsNeighbour)
  CALL SyncNeighbours(ParEnv)

  ! -
  ALLOCATE(SplittedMatrix % IfLCols(ParEnv % PEs))
  ALLOCATE(SplittedMatrix % IfORows(ParEnv % PEs))
  DO i = 1, ParEnv % PEs
     IF ( SplittedMatrix % IfMatrix(i) % NumberOfRows /= 0 ) THEN
        ALLOCATE( SplittedMatrix % IfVecs(i) % IfVec( &
             SplittedMatrix % IfMatrix(i) % NumberOfRows ), &
             SplittedMatrix % IfMatrix(i) % Values( SIZE(  &
             SplittedMatrix % IfMatrix(i) % Cols)), &
             SplittedMatrix % IfLCols(i) % IfVec( SIZE( &
             SplittedMatrix % IfMatrix(i) % Cols ) ), &
             SplittedMatrix % IfORows(i) % IfVec( &
             SplittedMatrix % IfMatrix(i) % NumberOfRows))
        SplittedMatrix % IfORows(i) % IfVec   = 0
        SplittedMatrix % IfLCols(i) % IfVec   = 0
        SplittedMatrix % IfVecs(i) % IfVec    = 0
        SplittedMatrix % IfMatrix(i) % Values = 0._dp

        IF ( NeedPrec ) THEN
           ALLOCATE( SplittedMatrix % IfMatrix(i) % PrecValues(  &
              SIZE(SplittedMatrix % IfMatrix(i) % Cols) ) )
           SplittedMatrix % IfMatrix(i) % PrecValues = 0._dp
        END IF

        IF ( NeedMass ) THEN
           ALLOCATE( SplittedMatrix % IfMatrix(i) % MassValues(  &
              SIZE(SplittedMatrix % IfMatrix(i) % Cols) ) )
           SplittedMatrix % IfMatrix(i) % MassValues = 0._dp
        END IF

        IF ( NeedDamp ) THEN
           ALLOCATE( SplittedMatrix % IfMatrix(i) % DampValues(  &
              SIZE(SplittedMatrix % IfMatrix(i) % Cols) ) )
           SplittedMatrix % IfMatrix(i) % DampValues = 0._dp
        END IF

        IF ( NeedILU ) THEN
           ALLOCATE( SplittedMatrix % IfMatrix(i) % ILUValues(  &
              SIZE(SplittedMatrix % IfMatrix(i) % Cols) ) )
           SplittedMatrix % IfMatrix(i) % ILUValues = 0._dp
        END IF
     END IF
  END DO
!if ( parenv % mype==0 ) print*, 'COMBINE AND STUFF TIME: ', realtime()-st; st=realtime()

  !----------------------------------------------------------------------
  !
  ! Renumber the Degrees of Freedom in SplittedMatrix % InsideMatrix
  !
  !----------------------------------------------------------------------
  CALL RenumberDOFs( SourceMatrix, SplittedMatrix, ParallelInfo )
!if ( parenv % mype==0 ) print*, 'RENUMBER TIME: ', realtime()-st; st=realtime()

  !----------------------------------------------------------------------
  !
  ! Build indirect indexing for vectors to speed up mat-vec multiply
  !
  !----------------------------------------------------------------------
  CALL BuildRevVecIndices( SplittedMatrix )
!if ( parenv % mype==0 ) print*, 'BUILD REV TIME: ', realtime()-st; st=realtime()

  !----------------------------------------------------------------------
  !
  ! Add incoming parts to insidematrix. Use ListMatrix as not all
  ! connections might be known in advance...
  !
  !----------------------------------------------------------------------
  sz = SIZE(A % Values)

  ! Check whether we need to create List matrix and add new column entries. 
  GotNewCol = .FALSE.
  DO i=1,Parenv % PEs
    CurrIF => SplittedMatrix % IfMatrix(i)
    DO j = 1, CurrIf % NumberOfRows
      IF ( Currif % RowOwner(j) /= ParEnv % MyPE ) CYCLE
      RowInd = SplittedMatrix % IfORows(i) % IfVec(j)
      IF ( rowind<=0 ) CYCLE
      DO k = CurrIf % Rows(j), CurrIf % Rows(j+1) - 1
        ColInd = SplittedMatrix % IfLCols(i) % IfVec(k)
        IF ( colind<=0 ) CYCLE
        IF( .NOT. CRS_CheckMatrixElement(A,RowInd,ColInd) ) THEN
          GotNewCol = .TRUE.
          GOTO 1
        END IF        
      END DO      
    END DO
  END DO

1 IF(GotNewCol) THEN
    CALL List_toListMatrix(A)
    DO i=1,Parenv % PEs
      CurrIF => SplittedMatrix % IfMatrix(i)
      DO j = 1, CurrIf % NumberOfRows
        IF ( Currif % RowOwner(j) /= ParEnv % MyPE ) CYCLE
        RowInd = SplittedMatrix % IfORows(i) % IfVec(j)
        IF ( rowind<=0 ) CYCLE
        DO k = CurrIf % Rows(j), CurrIf % Rows(j+1) - 1
          ColInd = SplittedMatrix % IfLCols(i) % IfVec(k)
          IF ( colind<=0 ) CYCLE
          l = l+1
          CALL List_AddMatrixIndex(A % ListMatrix,RowInd,ColInd)
        END DO
      END DO
    END DO
    CALL List_toCRSMatrix(A)
  END IF
    
  !----------------------------------------------------------------------
  !
  ! If need be, rebuild the inside part of GlueTable (place in the parallel
  ! 'insidematrix' where the 'partition' matrix entry is to be added).
  !
  !----------------------------------------------------------------------
  IF(sz /= SIZE(A % Values)) THEN
    ALLOCATE(Perm(A % NumberOfRows))
    j = 0;
    DO i=1,SourceMatrix % NumberOfRows
      IF (ParallelInfo % NeighbourList(i) % Neighbours(1)==ParEnv % MyPE) THEN 
        j=j+1
        Perm(j)=i
      END IF
    END DO

    DO i=1,A % NumberOfRows
      RowInd = Perm(i)
      K = SourceMatrix % Rows(RowInd)
      DO j=A % Rows(i),A % Rows(i+1)-1
        ColInd = Perm(A % Cols(j)) 
        DO WHILE(K<SourceMatrix % Rows(RowInd+1))
          IF(SourceMatrix % Cols(K)>=ColInd) THEN
            IF(SourceMatrix % Cols(K)==ColInd) &
              SplittedMatrix % GlueTable % Inds(K)=J
            EXIT
          END IF
          K=K+1
        END DO
      END DO
    END DO
    DEALLOCATE(Perm)
  END IF
!if ( parenv % mype==0 ) print*, 'REBUILD GLUETABLE: ', realtime()-st; st=realtime()

  ! Allocate more insidematrix stuff
  ! ---------------------------------
  sz = SIZE(A % Values)
  NULLIFY( A % PrecValues )
  IF ( NeedPrec ) ALLOCATE(A % PrecValues(sz))

  NULLIFY( A % MassValues )
  IF ( NeedMass ) ALLOCATE(A % MassValues(sz))

  NULLIFY( A % DampValues )
  IF ( NeedDamp ) ALLOCATE(A % DampValues(sz))

  A % ILUValues => NULL()
  IF ( NeedILU ) ALLOCATE(A % ILUValues(sz))


  !----------------------------------------------------------------------
  !
  ! Clean up temporary work space
  !
  !----------------------------------------------------------------------

  DEALLOCATE( OwnIfMRows, OwnIfMCols, OwnOCOR, OwnOldCols, &
       NbsIfMRows, NbsIfMCols, NbsOCOR, NbsOldCols )

  DO i= 1,ParEnv % PEs
     IF ( OwnIfMatrix(i) % NumberOfRows > 0 ) THEN
        DEALLOCATE( OwnIfMatrix(i) % Rows )
        DEALLOCATE( OwnIfMatrix(i) % Diag )
        DEALLOCATE( OwnIfMatrix(i) % Cols )
        DEALLOCATE( OwnIfMatrix(i) % GRows )
        DEALLOCATE( OwnIfMatrix(i) % RowOwner )
     END IF

     IF ( RecvdIfMatrix(i) % NumberOfRows > 0 ) THEN
        DEALLOCATE( RecvdIfMatrix(i) % Rows )
        DEALLOCATE( RecvdIfMatrix(i) % Diag )
        DEALLOCATE( RecvdIfMatrix(i) % Cols )
        DEALLOCATE( RecvdIfMatrix(i) % GRows )
        DEALLOCATE( RecvdIfMatrix(i) % RowOwner )
     END IF
  END DO

  DEALLOCATE( OwnIfMatrix, RecvdIfMatrix )
!if ( parenv % mype==0 ) print*, 'REST TIME: ', realtime()-st; st=realtime()
!********************************************************************
END FUNCTION SplitMatrix
!********************************************************************



!----------------------------------------------------------------------
!> Zero the splitted matrix (for new non-linear iteration)
!----------------------------------------------------------------------
SUBROUTINE ZeroSplittedMatrix( SplittedMatrix )
!----------------------------------------------------------------------
  USE Types
  IMPLICIT NONE
!----------------------------------------------------------------------
  TYPE (SplittedMatrixT), POINTER :: SplittedMatrix
!----------------------------------------------------------------------

  ! Local variables

  INTEGER :: i

  LOGICAL :: NeedMass, NeedDamp, NeedPrec, NeedILU

  !*******************************************************************

  NeedMass = ASSOCIATED(SplittedMatrix % InsideMatrix % MassValues)
  NeedDamp = ASSOCIATED(SplittedMatrix % InsideMatrix % DampValues)
  NeedPrec = ASSOCIATED(SplittedMatrix % InsideMatrix % PrecValues)

  SplittedMatrix % InsideMatrix % Values = 0._dp
  IF ( NeedMass ) SplittedMatrix % InsideMatrix % MassValues = 0._dp
  IF ( NeedDamp ) SplittedMatrix % InsideMatrix % DampValues = 0._dp
  IF ( NeedPrec ) SplittedMatrix % InsideMatrix % PrecValues = 0._dp

  NeedILU = .FALSE.
  DO i = 1, ParEnv % PEs
     IF ( SplittedMatrix % IfMatrix(i) % NumberOfRows /= 0 ) THEN
       NeedILU = ALLOCATED(SplittedMatrix % IfMatrix(i) % ILUValues)
       SplittedMatrix % IfMatrix(i) % Values = 0._dp
       IF ( NeedILU ) SplittedMatrix % IfMatrix(i) % ILUValues  = 0._dp
       IF ( NeedMass.AND.ALLOCATED(SplittedMatrix % IfMatrix(i) % MassValues) ) &
         SplittedMatrix % IfMatrix(i) % MassValues = 0._dp
       IF ( NeedPrec.AND.ALLOCATED(SplittedMatrix % IfMatrix(i) % PrecValues) ) &
         SplittedMatrix % IfMatrix(i) % PrecValues = 0._dp
       IF ( NeedDamp.AND.ALLOCATED(SplittedMatrix % IfMatrix(i) % DampValues) ) &
         SplittedMatrix % IfMatrix(i) % DampValues = 0._dp
     END IF

     IF ( SplittedMatrix % NbsIfMatrix(i) % NumberOfRows /= 0 ) THEN
       IF(ALLOCATED(SplittedMatrix % NbsIfMatrix(i) % Values)) &
         SplittedMatrix % NbsIfMatrix(i) % Values = 0._dp

       IF ( NeedILU.AND.ALLOCATED(SplittedMatrix % NbsIfMatrix(i) % ILUvalues) ) &
         SplittedMatrix % NbsIfMatrix(i) % ILUValues  = 0._dp
       IF ( NeedPrec.AND.ALLOCATED(SplittedMatrix % NbsIfMatrix(i) % Precvalues) ) &
         SplittedMatrix % NbsIfMatrix(i) % PrecValues  = 0._dp
       IF ( NeedMass.AND.ALLOCATED(SplittedMatrix % NbsIfMatrix(i) % MassValues) ) &
         SplittedMatrix % NbsIfMatrix(i) % MassValues = 0._dp
       IF ( NeedDamp.AND.ALLOCATED(SplittedMatrix % NbsIfMatrix(i) % DampValues) ) &
         SplittedMatrix % NbsIfMatrix(i) % DampValues = 0._dp
     END IF
  END DO
  IF(NeedILU) SplittedMatrix % InsideMatrix % ILUValues = 0._dp
!----------------------------------------------------------------------
END SUBROUTINE ZeroSplittedMatrix
!----------------------------------------------------------------------



!----------------------------------------------------------------------
! Get Hypre parameters in one central place. 
!----------------------------------------------------------------------
  SUBROUTINE SetHypreParameters(Params, hypremethod, Ilun, hypre_dppara, hypre_intpara )
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: hypremethod
    INTEGER :: ilun
    REAL(KIND=dp) :: hypre_dppara(10), r
    INTEGER :: hypre_intpara(20)

    CHARACTER(:), ALLOCATABLE, TARGET :: PrecMethod, IterativeMethod    
    CHARACTER(:), POINTER :: Method
    INTEGER :: hypre_sol, hypre_pre, i, j, k
    LOGICAL :: BPC, Found
    CHARACTER(*), PARAMETER :: Caller = 'HypreParameters' 
       
    CALL Info(Caller,'Setting parameters for Hypre solvers',Level=12)

    hypremethod = 0
    ilun = 0
    hypre_dppara = 0.0_dp
    hypre_intpara = 0

    hypre_sol = 0
    hypre_pre = 0

    !---------------------------------------
    !              No  Prec
    ! none         0    x   -
    ! BoomerAMG    1    x   x 
    ! AMS          2    x   x
    ! ILU          3    x   x
    ! Parasails    4    x   -
    ! FSAI         5    x   -
    ! PCG          6    -   x
    ! BiCGStab     7    -   x
    ! GMRes        8    -   x
    ! FlexGMRes    9    -   x
    ! LGMRes       10   -   x
    ! COGMRes      11   -   x
    !---------------------------------------    

    ! 1st round set the method
    ! 2nd round set the preconditioner!
    ! Many hypre methods can be both preconditioners and solvers. 
    DO i=0,1
      j = 0
      
      IF(i==0) THEN        
        hypre_sol = ListGetInteger( Params,'Linear System Method Hypre Index',Found )
        IF(Found ) CYCLE                    
        IterativeMethod = ListGetString( Params,'Linear System Iterative Method' )      
        Method => IterativeMethod
      ELSE
        hypre_pre = ListGetInteger( Params,'Linear System Preconditioning Hypre Index',Found )

        IF(Found ) CYCLE                    
        PrecMethod = ListGetString(Params,'Linear System Preconditioning',Found)

        IF(.NOT. Found) THEN
          CALL Info(Caller,'"Linear System Preconditioning" not given, assuming "none"')
          PrecMethod = "none"
        END IF
        Method => PrecMethod
      END IF
      
      IF( Method == 'none' ) THEN
        j = 0
      ELSE IF( Method == 'boomeramg' ) THEN
        j = 1
      ELSE IF( Method == 'ams' ) THEN
        j = 2
      ELSE IF ( SEQL(Method,'ilu') ) THEN
        Ilun = 0
        IF(LEN(Method)>=4) READ( Method(4:), *, END=10 ) ILUn
10      CONTINUE
        j = 3
      ELSE IF( Method == 'parasails' ) THEN
        j = 4
      ELSE IF( Method == 'fsai' ) THEN
        j = 5
      ELSE IF( Method == 'cg' .OR. Method == 'pcg' ) THEN
        j = 6
      ELSE IF( Method == 'bicgstab' ) THEN
        j = 7
      ELSE IF( Method == 'gmres' ) THEN
        j = 8
      ELSE IF( Method == 'flexgmres' ) THEN
        j = 9
      ELSE IF( Method == 'lgmres' ) THEN
        j = 10
      ELSE IF( Method == 'cogmres' ) THEN
        j = 11
      ELSE
        CALL Fatal(Caller,'Invalid method for Hypre: '//TRIM(Method))
      END IF
      
      IF(i==0) THEN
        CALL Info(Caller,'Using HYPRE iterative method: '//TRIM(IterativeMethod),Level=7)
        hypre_sol = j
      ELSE
        CALL Info(Caller,'Using HYPRE preconditioning method: '//TRIM(PrecMethod),Level=7)
        hypre_pre = j
      END IF
    END DO
    
    ! Some methods (Krylov methods) can not act as preconditioners!
    IF( ANY( hypre_pre == [6,7,8,9,10,11] ) ) THEN
      CALL Fatal(Caller,'Invalid preconditioner for Hypre: '//TRIM(PrecMethod))
    END IF

    ! Some methods cannot act as solvers!
    IF( ANY( hypre_sol == [0,4,5] ) ) THEN
      CALL Fatal(Caller,'Invalid solver for Hypre: '//TRIM(IterativeMethod))
    END IF

    ! We map the precondtioner + solver to one figure. 
    hypremethod = 100 * hypre_sol + hypre_pre
    CALL Info(Caller,'Hypre method index: '//I2S(hypremethod),Level=6)
    
    DO i=0,1
      IF(i==0) THEN
        j = hypre_pre
      ELSE
        j = hypre_sol
      END IF
      
      SELECT CASE(j)
      CASE(1) ! BoomerAMG
        hypre_intpara(1) = ListGetInteger( Params,&
            'BoomerAMG Relax Type', Found, DefValue = 3 )

        hypre_intpara(2) = ListGetInteger( Params,&
            'BoomerAMG Coarsen Type', Found, DefValue = 0)

        hypre_intpara(3) = ListGetInteger( Params,&
            'BoomerAMG Num Sweeps', Found, DefValue = 1 )

        hypre_intpara(4) = ListGetInteger( Params,&
            'Boomeramg Max Levels', Found, DefValue = 25 )

        hypre_intpara(5) = ListGetInteger( Params,&
            'BoomerAMG Interpolation Type', Found, DefValue = 0 )

        hypre_intpara(6) = ListGetInteger( Params,&
            'BoomerAMG Smooth Type', Found, DefValue = 6 )

        hypre_intpara(7) = ListGetInteger( Params,&
            'BoomerAMG Cycle Type', Found, DefValue = 1)

        BPC = ListGetLogical( Params, 'Block Preconditioner', Found )

        hypre_intpara(8) = ListGetInteger( Params, 'BoomerAMG Num Functions', Found )
        k = CurrentModel % Solver % Variable % DOFs            
        IF (.NOT.Found)  THEN
          IF (BPC) THEN
            hypre_intpara(8) = 1
          ELSE
            hypre_intpara(8) = k
          END IF
        ELSE IF ( .NOT.BPC .AND. (hypre_intpara(8) /= k) ) THEN
          WRITE(Message,'(A,I0,A,I0)') 'Read > BoomerAMG Num Functions < value ',&
              hypre_intpara(8), ', not equal to DOFs of solver variable ',k
          CALL Warn(Caller,Message)
        END IF

        hypre_dppara(1) = ListGetCReal( Params,&
            'BoomerAMG Strong Threshold', Found, DefValue = 0.25_dp)

      CASE(2) ! AMS
        ! The numbering follows the old convention. 
        hypre_intpara(1) = ListGetInteger( Params,&
            'AMS Max Iterations', Found, DefValue = 1 )

        hypre_dppara(1) = ListGetCReal( Params,&
            'AMS Tolerance', Found, DefValue = 1.0e-6_dp )
        
        hypre_intpara(2) = ListGetInteger( Params,&
            'AMS Cycle Type', Found, DefValue = 1)  ! 1-14
       
        hypre_intpara(3) = ListGetInteger( Params,&
            'AMS Relax Type', Found, DefValue = 2 )
        hypre_intpara(4) = ListGetInteger( Params,&
            'AMS Relax Times', Found, DefValue = 1 )

        hypre_dppara(2) = ListGetCReal( Params,&
            'AMS Relax Weight', Found, DefValue = 1.0_dp )
        hypre_dppara(3) = ListGetCReal( Params,&
            'AMS Relax Omega', Found, DefValue = 1.0_dp )
        hypre_dppara(4) = ListGetCReal( Params,&
            'AMS Alpha Threshold', Found, DefValue = 0.25_dp )
        hypre_dppara(5) = ListGetCReal( Params,&
            'AMS Beta Threshold', Found, DefValue = 0.25_dp )
        
        IF( ListGetLogical( Params,&
            'AMS Singular Matrix', Found ) ) hypre_intpara(5) = 1
                
      CASE(3) ! ILU
        hypre_intpara(1) = ListGetInteger( Params,&
            'ILU Max Iterations', Found, DefValue = 1 )

        hypre_intpara(2) = ListGetInteger( Params,&
            'ILU Reordering', Found, DefValue = 0 )

        hypre_intpara(3) = ListGetInteger( Params,&
            'ILU Max Nonzeros Per Row', Found, DefValue = 0 )

        hypre_intpara(4) = ListGetInteger( Params,&
            'ILU Schur Max Iterations', Found, DefValue = 0 )

        hypre_intpara(5) = ListGetInteger( Params,&
            'ILU Trisolve Direct', Found, DefValue = 0 )
        
        hypre_dppara(1) = ListGetCReal( Params,&
            'ILU Tolerance', Found, DefValue = 0.0_dp)

        hypre_dppara(2) = ListGetCReal( Params,&
            'ILU Drop Threshold', Found, DefValue = 0.0_dp)

        hypre_dppara(3) = ListGetCReal( Params,&
            'ILU NSH Drop Threshold', Found, DefValue = 0.0_dp)     
        
      CASE(4) ! Parasails
        hypre_intpara(1) = ListGetInteger( Params,&
            'ParaSails Symmetry', Found, DefValue = 0 )

        hypre_intpara(2) = ListGetInteger( Params,&
            'ParaSails Maxlevel', Found, DefValue = 1)

        hypre_dppara(1) = ListGetCReal( Params,&
            'ParaSails Threshold', Found, DefValue = -0.95_dp)

        hypre_dppara(2) = ListGetCReal( Params,&
            'ParaSails Filter', Found, DefValue = -0.95_dp)

      CASE(5) ! FSAI
        hypre_intpara(1) = ListGetInteger( Params,&
            'FSAI Max Steps', Found, DefValue = 5 )
        hypre_intpara(2) = ListGetInteger( Params,&
            'FSAI Max Step Size', Found, DefValue = 3 )
        
        hypre_dppara(1) = ListGetCReal( Params,&
            'FSAI Kap Tolerance', Found, DefValue = 1.0e-3_dp )
        

      CASE(6) ! PCG 
        IF( ListGetLogical(Params,'PCG Two Norm',Found ) ) THEN
          hypre_intpara(11) = 1
        END IF
        !IF( ListGetLogical(Params,'PCG Set Flex',Found ) ) THEN
        !  hypre_intpara(12) = 1
        !END IF

      CASE(7) ! BiCGStab

      CASE(8) ! GMRes        

      CASE(9) ! FlexGMRes
        
      CASE(10) ! LGMRes
        hypre_intpara(11) = ListGetInteger( Params,&
            'HYPRE LGmRes Aug Dim', Found, DefValue = 2 )
        
      CASE(11) ! COGGMRes
        hypre_intpara(11) = ListGetInteger( Params,&
            'HYPRE COGMRes Unroll', Found, DefValue = 0 )
        hypre_intpara(12) = ListGetInteger( Params,&
            'HYPRE COGMRes CGS', Found, DefValue = 0 )
        
      END SELECT

      ! These are used by many
      IF(j>5) THEN
        hypre_dppara(6) = ListGetCReal( Params,&
            'HYPRE Absolute Tolerance', Found, DefValue = 0.0_dp )
      END IF

      IF(j>7) THEN
        hypre_intpara(10) = ListGetInteger( Params,&
            'HYPRE GmRes Dimension', Found, DefValue = 100 )        
      END IF
      
    END DO
      
  
  END SUBROUTINE SetHypreParameters



!----------------------------------------------------------------------
  SUBROUTINE SParInitSolve( SourceMatrix,XVec,RHSVec,RVec, &
          ParallelInfo, UpdateMatrix )
!----------------------------------------------------------------------
! Initialize Parallel Solve
!----------------------------------------------------------------------

    TYPE (Matrix_t) :: SourceMatrix
    LOGICAL :: UpdateMatrix
    TYPE (ParallelInfo_t), TARGET :: ParallelInfo
    REAL(KIND=dp), DIMENSION(:) :: XVec, RHSVec,RVec

!----------------------------------------------------------------------

    ! Local variables

    REAL(KIND=dp), POINTER :: TmpRHSVec(:)
    INTEGER :: i, j, k, l, grow, gcol
    INTEGER :: nodeind, ifind, dd, rowind, ierr
    TYPE (BasicMatrix_t), POINTER :: CurrIf
    TYPE (GlueTableT), POINTER :: GT
    TYPE (SplittedMatrixT), POINTER :: SplittedMatrix

    LOGICAL :: NeedMass, NeedDamp, NeedPrec, NeedILU, found
!----------------------------------------------------------------------

    GlobalData     => SourceMatrix % ParMatrix
    ParEnv         => GlobalData % ParEnv
    ParEnv % ActiveComm = SourceMatrix % Comm
    SplittedMatrix => Globaldata % SplittedMatrix

    GlobalData % DOFs  = 1
    GlobalData % ParallelInfo => ParallelInfo

    IF ( UpdateMatrix ) THEN
      NeedMass = ASSOCIATED(SplittedMatrix % InsideMatrix % MassValues)
      NeedDamp = ASSOCIATED(SplittedMatrix % InsideMatrix % DampValues)
      NeedPrec = ASSOCIATED(SplittedMatrix % InsideMatrix % PrecValues)

      NeedILU = .FALSE.
      DO i=1,ParEnv % PEs
        IF (SplittedMatrix % IfMatrix(i) % NumberOfRows>0) THEN
          NeedILU = ALLOCATED(SplittedMatrix % IfMatrix(i) % ILUValues)
          EXIT
        END IF
      END DO

      CALL ZeroSplittedMatrix( SplittedMatrix )
      !------------------------------------------------------------------
      !
      ! Copy the Matrix % Values into SplittedMatrix
      !
      !------------------------------------------------------------------

      GT => SplittedMatrix % GlueTable
      DO i = 1, SourceMatrix % NumberOfRows
       
         GRow = ParallelInfo % GlobalDOFs(i)
         DO j = SourceMatrix % Rows(i),SourceMatrix % Rows(i+1) - 1

            GCol = ParallelInfo % GlobalDOFs(SourceMatrix % Cols(j))

            IF ( GT % Inds(j)>0 ) THEN
               SplittedMatrix % InsideMatrix % Values( GT % Inds(j) ) = &
                    SplittedMatrix % InsideMatrix % Values( &
                       GT % Inds(j) ) + SourceMatrix % Values(j)

               IF ( NeedPrec ) &
                  SplittedMatrix % InsideMatrix % PrecValues( GT % Inds(j) ) = &
                       SplittedMatrix % InsideMatrix % PrecValues( &
                        GT % Inds(j) ) + SourceMatrix % PrecValues(j)
               IF ( NeedMass ) &
                  SplittedMatrix % InsideMatrix % MassValues( GT % Inds(j) ) = &
                       SplittedMatrix % InsideMatrix % MassValues( &
                        GT % Inds(j) ) + SourceMatrix % MassValues(j)
               IF ( NeedDamp ) &
                  SplittedMatrix % InsideMatrix % DampValues( GT % Inds(j) ) = &
                       SplittedMatrix % InsideMatrix % DampValues( &
                          GT % Inds(j) ) + SourceMatrix % DampValues(j)
               IF ( NeedILU ) &
                  SplittedMatrix % InsideMatrix % ILUValues( GT % Inds(j) ) = &
                       SplittedMatrix % InsideMatrix % ILUValues( &
                          GT % Inds(j) ) + SourceMatrix % ILUValues(j)

            ELSE IF ( (GT % Inds(j) + ParEnv % PEs) >= 0 ) THEN

               ifind = ABS(GT % Inds(j))
               CurrIf => SplittedMatrix % IfMatrix(ifind)

               RowInd = -1
               IF ( CurrIf % NumberOfRows > 0 ) THEN
                  RowInd = SearchIAItem( CurrIf % NumberOfRows, CurrIf % GRows, GRow )
               END IF

               IF ( RowInd /= -1 ) THEN
                  DO l = CurrIf % Rows(rowind), CurrIf % Rows(rowind+1)-1

                     IF ( GCol == CurrIf % Cols(l) ) THEN
                        CurrIf % Values(l) = CurrIf % Values(l) + SourceMatrix % Values(j)
                        IF ( NeedPrec ) &
                          CurrIf % PrecValues(l) = CurrIf % PrecValues(l) + &
                               SourceMatrix % PrecValues(j)
                        IF ( NeedMass ) &
                          CurrIf % MassValues(l) = CurrIf % MassValues(l) + &
                               SourceMatrix % MassValues(j)
                        IF ( NeedDamp ) &
                          CurrIf % DampValues(l) = CurrIf % DampValues(l) + &
                               SourceMatrix % DampValues(j)
                        IF ( NeedILU ) &
                          CurrIf % ILUValues(l) = CurrIf % ILUValues(l) + &
                               SourceMatrix % ILUValues(j)
                        EXIT
                     END IF

                  END DO
               END IF

            ELSE IF ((GT % Inds(j) + (2*ParEnv % PEs)) >= 0) THEN

               ifind = -ParEnv % PEs + ABS(GT % Inds(j))
               CurrIf => SplittedMatrix % NbsIfMatrix(ifind)
                  
               RowInd = -1
               IF ( CurrIf % NumberOfRows > 0 ) THEN
                  RowInd = SearchIAItem( CurrIf % NumberOfRows, CurrIf % GRows, GRow )
               END IF

               IF ( RowInd /= -1 ) THEN
                  DO l = CurrIf % Rows(RowInd), CurrIf % Rows(RowInd+1)-1
                     IF ( GCol == CurrIf % Cols(l) ) THEN
                        CurrIf % Values(l) = CurrIf % Values(l) + SourceMatrix % Values(j)
                        IF ( NeedPrec ) &
                           CurrIf % PrecValues(l) = CurrIf % PrecValues(l) + &
                                SourceMatrix % PrecValues(j)
                        IF ( NeedMass ) &
                           CurrIf % MassValues(l) = CurrIf % MassValues(l) + &
                                SourceMatrix % MassValues(j)
                        IF ( NeedDamp ) &
                           CurrIf % DampValues(l) = CurrIf % DampValues(l) + &
                                SourceMatrix % DampValues(j)
                        IF ( NeedILU ) &
                           CurrIf % ILUValues(l) = CurrIf % ILUValues(l) + &
                                SourceMatrix % ILUValues(j)
                        EXIT
                     END IF
                  END DO
               END IF
            END IF
         END DO
      END DO

      CALL GlueFinalize( SourceMatrix, SplittedMatrix, ParallelInfo )

      !
      ! Initialize Right-Hand-Side:
      ! ---------------------------
      IF ( .NOT. ASSOCIATED( SplittedMatrix % TmpXVec ) ) THEN
         ALLOCATE( SplittedMatrix % InsideMatrix % RHS(  &
           SplittedMatrix % InsideMatrix % NumberOfRows ) )
         SPlittedMatrix % InsideMatrix % RHS = 0
  
         ALLOCATE( SplittedMatrix % TmpXVec( SplittedMatrix %  &
                    InsideMatrix % NumberOfRows ) )
         SplittedMatrix % TmpXVec = 0

         ALLOCATE( SplittedMatrix % TmpRVec( SplittedMatrix %  &
                    InsideMatrix % NumberOfRows ) )
         SplittedMAtrix % TmpRVec = 0
      END IF
    END IF

    CALL SParUpdateRHS( SourceMatrix, RHSVec, ParallelInfo )
 
    !
    ! Initialize temporary XVec and RVec for iterator. The
    ! originals contain also the items on interfaces.
    ! -------------------------------------------------------------
    CALL SParUpdateSolve( SourceMatrix, XVec, RVec )
    !
    ! Set up the preconditioner:
    ! --------------------------
    IF ( .NOT.SPlittedMatrix % InsideMatrix % Ordered ) &
       CALL CRS_SortMatrix( SplittedMatrix % InsideMatrix )
!----------------------------------------------------------------------
  END SUBROUTINE SParInitSolve
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  SUBROUTINE SParUpdateRHS( SourceMatrix, RHSVec, ParallelInfo )
!----------------------------------------------------------------------
    TYPE (Matrix_t) :: SourceMatrix
    TYPE(ParallelInfo_t) :: ParallelInfo
    REAL(KIND=dp) :: RHSVec(:)
!----------------------------------------------------------------------
    ! Local variables

    REAL(KIND=dp), POINTER :: TmpRHSVec(:)
    TYPE (SplittedMatrixT), POINTER :: SplittedMatrix

!----------------------------------------------------------------------
    SplittedMatrix => SourceMatrix % ParMatrix % SplittedMatrix

    IF (.NOT. ASSOCIATED(SplittedMatrix % InsideMatrix % RHS)) &
         ALLOCATE(SplittedMatrix % InsideMatrix % RHS(SplittedMatrix % InsideMatrix % NumberOfRows))
    TmpRHSVec => SplittedMatrix % InsideMatrix % RHS

    CALL ExchangeRHSIf( SourceMatrix, SplittedMatrix, &
            ParallelInfo, RHSVec, TmpRHSVec )
!----------------------------------------------------------------------
  END SUBROUTINE SParUpdateRHS
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  SUBROUTINE SParUpdateSolve( SourceMatrix, x, r )
!----------------------------------------------------------------------
    REAL(KIND=dp) :: x(:),r(:)
    TYPE(Matrix_t) :: SourceMatrix
!----------------------------------------------------------------------
    TYPE(ParallelInfo_t), POINTER :: ParallelInfo
    INTEGER :: i,j,k
    REAL(KIND=dp), POINTER :: TmpXVec(:),TmpRVec(:)
!----------------------------------------------------------------------
    ParallelInfo => SourceMatrix % ParMatrix % ParallelInfo
    IF (.NOT. ASSOCIATED(SourceMatrix % ParMatrix % SplittedMatrix % TmpXVec)) &
         ALLOCATE(SourceMatrix % ParMatrix % SplittedMatrix % TmpXVec( &
	 SourceMatrix % ParMatrix % SplittedMatrix % InsideMatrix % NumberOfRows))
    TmpXVec => SourceMatrix % ParMatrix % SplittedMatrix % TmpXVec
    IF (.NOT. ASSOCIATED(SourceMatrix % ParMatrix % SplittedMatrix % TmpRVec)) &
         ALLOCATE(SourceMatrix % ParMatrix % SplittedMatrix % TmpRVec( &
	 SourceMatrix % ParMatrix % SplittedMatrix % InsideMatrix % NumberOfRows))
    TmpRVec => SourceMatrix % ParMatrix % SplittedMatrix % TmpRVec

    j = 0
    DO i = 1, SourceMatrix % NumberOfRows
       IF ( ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
          j = j + 1
          TmpXVec(j) = x(i)
          TmpRVec(j) = r(i)
       END IF
    END DO
!----------------------------------------------------------------------
  END SUBROUTINE SParUpdateSolve
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  SUBROUTINE SParUpdateResult( SourceMatrix, XVec, RVec, GlobalUpdate )
!----------------------------------------------------------------------
    TYPE (Matrix_t) :: SourceMatrix
    LOGICAL :: GlobalUpdate
    REAL(KIND=dp) :: XVec(:), RVec(:)
!----------------------------------------------------------------------
    ! Local variables

    REAL(KIND=dp), POINTER :: TmpXVec(:), TmpRVec(:)
    INTEGER, ALLOCATABLE :: VecEPerNB(:)
    TYPE(ParallelInfo_t), POINTER :: ParallelInfo
    INTEGER :: i, j, k, nbind
    TYPE (SplittedMatrixT), POINTER :: SplittedMatrix

!----------------------------------------------------------------------

    ! Collect the result:
    ! -------------------
    SplittedMatrix => SourceMatrix % ParMatrix % SplittedMatrix
    ParallelInfo => SourceMatrix % ParMatrix % ParallelInfo

    TmpXVec => SplittedMatrix % TmpXVec
    TmpRVec => SplittedMatrix % TmpRVec

    IF(.NOT. ASSOCIATED(SourceMatrix % Rhs) ) THEN
      CALL Info('SParUpdateResult','Rhs is not yet associated!?',Level=20)
      ALLOCATE( SourceMatrix % Rhs(SourceMatrix % NumberOfRows) )
      SourceMatrix % Rhs = 0.0_dp
    END IF
           
    j = 0
    DO i = 1, SourceMatrix % NumberOfRows
       IF ( ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
         j = j + 1
         XVec(i) = TmpXVec(j)
         RVec(i) = TmpRVec(j)
       ELSE
         Rvec(i) = SourceMatrix % RHS(i)
       END IF
    END DO

    IF ( .NOT. GlobalUpdate ) RETURN

    ALLOCATE( VecEPerNB( ParEnv % PEs ) )
    VecEPerNB = 0
    DO i = 1, SourceMatrix % NumberOfRows
       IF ( SIZE(ParallelInfo % NeighbourList(i) % Neighbours) > 1 ) THEN
          IF ( ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
             DO j = 1, SIZE(ParallelInfo % NeighbourList(i) % Neighbours)
               IF (ParallelInfo % NeighbourList(i) % Neighbours(j)/=ParEnv % MyPE) THEN
                 nbind = ParallelInfo % NeighbourList(i) % Neighbours(j) + 1
                 VecEPerNB(nbind) = VecEPerNB(nbind) + 1

                 SplittedMatrix % ResBuf(nbind) % ResVal(VecEPerNB(nbind)) = &
                      XVec(i)
                 SplittedMatrix % ResBuf(nbind) % ResInd(VecEPerNB(nbind)) = &
                   ParallelInfo % GlobalDOFs(i)
               END IF
             END DO
          END IF
       END IF
    END DO

    CALL ExchangeResult( SourceMatrix, SplittedMatrix, ParallelInfo, XVec )

#if 0
    VecEPerNB = 0
    DO i = 1, SourceMatrix % NumberOfRows
       k = (SourceMatrix % INVPerm(i) + DOFs-1) / DOFs
       IF ( SIZE(ParallelInfo % NeighbourList(k) % Neighbours) > 1 ) THEN
          IF ( ParallelInfo % NeighbourList(k) % Neighbours(1) == ParEnv % MyPE ) THEN
             DO j = 1, SIZE(ParallelInfo % NeighbourList(k) % Neighbours)
               IF (ParallelInfo % NeighbourList(k) % Neighbours(j)/=ParEnv % MyPE) THEN
                 nbind = ParallelInfo % NeighbourList(k) % Neighbours(j) + 1
                 VecEPerNB(nbind) = VecEPerNB(nbind) + 1

                 SplittedMatrix % ResBuf(nbind) % ResVal(VecEPerNB(nbind)) = &
                      RVec(i)
                 SplittedMatrix % ResBuf(nbind) % ResInd(VecEPerNB(nbind)) = &
                   ParallelInfo % GlobalDOFs(k)*DOFs - (DOFs-1-MOD( i-1, DOFs))
               END IF
             END DO
          END IF
       END IF
     END DO

     CALL ExchangeResult( SourceMatrix, SplittedMatrix, ParallelInfo, RVec, DOFs ) 
#endif
     !
     ! Clean the work space:
     !----------------------
     DEALLOCATE( VecEPerNB )
!----------------------------------------------------------------------
  END SUBROUTINE SParUpdateResult
!----------------------------------------------------------------------



!----------------------------------------------------------------------
!> Create continuous numbering for the dofs expected by some linear solvers.
!----------------------------------------------------------------------
  SUBROUTINE ContinuousNumbering(ParallelInfo, Mperm, Aperm, Owner, &
              nin, Mesh, nOwn, gStart )
!--------------------------------------------------------------------
     INTEGER :: Mperm(:), Aperm(:), Owner(:)
     TYPE(Mesh_t), OPTIONAL :: Mesh
     INTEGER, OPTIONAL :: nin
     TYPE(ParallelInfo_t) :: ParallelInfo
     INTEGER, OPTIONAL :: nOwn, gStart

     INTEGER, ALLOCATABLE :: neigh(:), sz(:),buf_a(:,:),buf_g(:,:), &
                  buf_aa(:), buf_gg(:)
     INTEGER, ALLOCATABLE :: n_nownbuf(:),g_nownbuf(:),i_nownbuf(:)
     INTEGER, POINTER :: nb(:)
     LOGICAL, POINTER :: isNeighbour(:)
     INTEGER :: min_id,max_id,my_id,next_id,prev_id, active_neighbours
     INTEGER :: i,j,k,src,n,nob,gind,gindp,ssz,status(MPI_STATUS_SIZE),ierr,nneigh

     Owner = 0; Aperm = 0;
     IF(PRESENT(nin)) THEN
       n=nin
     ELSE
       n = SIZE(ParallelInfo % GlobalDOFs)
     END IF

     IF(PRESENT(Mesh)) THEN
       ALLOCATE(isNeighbour(ParEnv % PEs))
       nneigh = MeshNeighbours(Mesh,isNeighbour)
     ELSE
       nneigh = ParEnv % NumOfNeighbours
       isNeighbour => ParEnv % IsNeighbour
     ENDIF
 
     DO min_id=0,Parenv % PEs-1
       IF ( ParEnv % Active(min_id+1) ) EXIT
     END DO

     DO max_id=Parenv % PEs-1,0,-1
       IF ( ParEnv % Active(max_id+1) ) EXIT
     END DO

     my_id = ParEnv % MyPE

     ! receive our lowest index:
     ! -------------------------
     gind = 0
     IF ( my_id /= min_id ) THEN
       DO prev_id=my_id-1,0,-1
         IF ( ParEnv % Active(prev_id+1) ) EXIT
       END DO
       CALL MPI_RECV( gind, 1, MPI_INTEGER, &
           prev_id, 801, ELMER_COMM_WORLD, status, ierr )
     END IF

     IF(PRESENT(gStart)) gStart = gind

     ! give a number to dofs owned by us:
     ! -----------------------------------
     gindp = gind
     DO i=1,n
       nb => ParallelInfo % NeighbourList(i) % Neighbours
       IF ( nb(1)==my_id .OR. ALL(nb/=my_id) ) THEN
         Owner(i) = 1
         gind = gind + 1
         Aperm(i) = gind
       END IF
     END DO


     ! Compute the number of dofs owned                                         
     IF (PRESENT(nOwn)) nOwn = gind - gindp

     ! next pe in line needs it's base:
     ! --------------------------------
     IF ( ParEnv % MyPE < max_id ) THEN
       DO next_id=my_id+1,ParEnv % PEs-1
         IF ( ParEnv % Active(next_id+1) ) EXIT
       END DO
       CALL MPI_BSEND( gind, 1, MPI_INTEGER, &
          next_id, 801, ELMER_COMM_WORLD, ierr )
     END IF

     ! the rest is to communicate the numbering of shared dofs
     ! from owners to others that share the dofs:
     ! --------------------------------------------------------
     ALLOCATE(neigh(Parenv % PEs),sz(nneigh))
     neigh=-1
     k = 0
     DO i=1,ParEnv % PEs
       IF ( .NOT. ParEnv % Active(i) ) CYCLE

       IF( IsNeighbour(i) ) THEN
         k = k + 1
         neigh(i) = k
       END IF
     END DO

     sz = 0
     DO i=1,n
       IF ( Owner(i)==1 ) THEN
         nb => ParallelInfo % NeighbourList(i) % Neighbours
         IF(nb(1) == my_id) THEN
           DO j=2,SIZE(nb)
             k = neigh(nb(j)+1)
             IF(k<=0) CYCLE
             sz(k) = sz(k) + 1
           END DO
         END IF
       END IF
     END DO


     ALLOCATE( buf_a(MAXVAL(sz),nneigh), &
               buf_g(MAXVAL(sz),nneigh))

     nob = COUNT(owner==0)
     DO i=1,n
       nb => ParallelInfo % NeighbourList(i) % Neighbours
       IF(Owner(i)==1.AND. nb(1)/=my_id) nob=nob+1
     END DO
     ALLOCATE( n_nownbuf(nob), g_nownbuf(nob), i_nownbuf(nob) )

     sz  = 0
     nob = 0
     DO i=1,n
       IF ( Owner(i)==1 ) THEN
         nb => ParallelInfo % NeighbourList(i) % Neighbours
         IF(nb(1)==my_id) THEN
           DO j=2,SIZE(nb)
             k = neigh(nb(j)+1)
             IF(k<=0) CYCLE
             sz(k) = sz(k) + 1
             buf_a(sz(k),k) = Aperm(i)
             buf_g(sz(k),k) = ParallelInfo % GlobalDOFs(i)
           END DO
         END IF
       ELSE
         nob = nob + 1
         i_nownbuf(nob) = nob
         n_nownbuf(nob) = i
         g_nownbuf(nob) = ParallelInfo % GlobalDOFs(i)
       END IF
     END DO
     IF ( nob>0 ) CALL SortI( nob,g_nownbuf,n_nownbuf )

     active_neighbours = 0
     DO i=1,ParEnv % PEs
       IF ( .NOT. ParEnv % Active(i) ) CYCLE

       IF ( isNeighbour(i) ) THEN
         active_neighbours = active_neighbours+1
         k = neigh(i)
         ssz = sz(k)
         CALL MPI_BSEND( ssz,1,MPI_INTEGER,i-1,802,ELMER_COMM_WORLD,ierr )
         IF ( ssz>0 ) THEN
           CALL MPI_BSEND( buf_a(1:ssz,k),ssz,MPI_INTEGER,i-1,803,ELMER_COMM_WORLD,ierr )
           CALL MPI_BSEND( buf_g(1:ssz,k),ssz,MPI_INTEGER,i-1,804,ELMER_COMM_WORLD,ierr )
         END IF
       END IF
     END DO 

     DEALLOCATE( buf_a, buf_g, neigh, sz )

     DO i=1,active_neighbours
       CALL MPI_RECV( ssz,1,MPI_INTEGER,MPI_ANY_SOURCE,802,ELMER_COMM_WORLD,status,ierr )
       IF ( ssz>0 ) THEN
         src = status(MPI_SOURCE)
         ALLOCATE( buf_aa(ssz), buf_gg(ssz) )

         CALL MPI_RECV( buf_aa,ssz,MPI_INTEGER,src,803,ELMER_COMM_WORLD,status,ierr )
         CALL MPI_RECV( buf_gg,ssz,MPI_INTEGER,src,804,ELMER_COMM_WORLD,status,ierr )

         DO j=1,ssz
           k = SearchIAItem( nob, g_nownbuf, buf_gg(j), i_nownbuf )
           IF ( k>0 ) Aperm(n_nownbuf(k)) = buf_aa(j)
         END DO

         DEALLOCATE( buf_aa, buf_gg )
       END IF
     END DO

     IF(PRESENT(Mesh)) DEALLOCATE(isNeighbour)
     DEALLOCATE( n_nownbuf, g_nownbuf, i_nownbuf )

     CALL SParIterActiveBarrier()
!--------------------------------------------------------------------
   END SUBROUTINE ContinuousNumbering
!--------------------------------------------------------------------




!--------------------------------------------------------------------
!> Call the parallel iterative solver
!--------------------------------------------------------------------
SUBROUTINE SParIterSolver( SourceMatrix, ParallelInfo, XVec, &
            RHSVec, Solver, SParMatrixDesc )

  USE, INTRINSIC :: iso_c_binding                

  TYPE (ParallelInfo_t) :: ParallelInfo
  TYPE (Matrix_t) :: SourceMatrix
  REAL(KIND=dp), DIMENSION(:) :: XVec, RHSVec
  TYPE (Solver_t) :: Solver
  TYPE (SParIterSolverGlobalD_t), POINTER :: SParMatrixDesc,SaveGlobalData

  ! Local variables

  TYPE (ErrInfoT) :: ErrInfo
  INTEGER :: i, j, n,grow,gcol, nbind
  TYPE (SplittedMatrixT), POINTER :: SplittedMatrix
  INTEGER :: k, l, nodeind, ifind, dd, rowind, ierr
  TYPE (BasicMatrix_t), POINTER :: CurrIf
  TYPE (GlueTableT), POINTER :: GT

  CHARACTER(:), ALLOCATABLE :: XmlFile
  REAL(KIND=dp) :: hypre_dppara(10) = 0, TOL
  INTEGER :: ILUn, BILU, Rounds, buf(2), src, status(MPI_STATUS_SIZE), ssz,nob, &
      hypremethod,  gind, hypre_intpara(20) = 0

  INTEGER, DIMENSION(:), ALLOCATABLE :: VecEPerNB
  INTEGER, ALLOCATABLE :: Owner(:), Aperm(:)

  REAL(KIND=dp), TARGET :: DummyPrecVals(1)

  REAL(KIND=dp), POINTER :: Vals(:)
  INTEGER, POINTER :: nb(:), Rows(:), Cols(:)

  LOGICAL :: NeedMass, NeedDamp, NeedPrec, NeedILU, Found
  LOGICAL :: NewSetup, UpdateTolerance
  INTEGER :: verbosity
  
  INTEGER :: nrows, ncols, nnz
  TYPE(ValueList_t), POINTER :: Params
  INTEGER,ALLOCATABLE::revdoflist(:)
  INTEGER::inside
   
#ifdef HAVE_TRILINOS
    INTERFACE
      !! create Trilinos matrix and setup solver/preconditioner
      SUBROUTINE SolveTrilinos1( n, nnz, Rows, Cols, Vals, &
           GDOFs, Owner, &
           xmlfilename, verbosity, triliContainer, n_nodes, &
           xcoords,ycoords,zcoords,ierr) BIND(C,name='SolveTrilinos1')
        
        USE, INTRINSIC :: iso_c_binding
        
        INTEGER(KIND=c_int) :: n ! number of nodes belonging to local elements
        INTEGER(KIND=c_int) :: nnz   ! number of local nonzeros
        INTEGER(KIND=c_int) :: Rows(n+1), Cols(nnz), GDOFs(n), &
                   PE, Owner(n)
        REAL(KIND=c_double) :: Vals(nnz)
        INTEGER(KIND=c_int) :: verbosity
        CHARACTER(c_char) :: xmlfilename
        INTEGER(KIND=C_INTPTR_T) :: triliContainer
        INTEGER(KIND=C_INT) :: n_nodes
        REAL(KIND=c_double) :: xcoords(n_nodes), ycoords(n_nodes), zcoords(n_nodes)
        INTEGER(KIND=C_INT) :: ierr
      END SUBROUTINE SolveTrilinos1

      !! solve linear system with same matrix as in SolveTrilinos1.
      SUBROUTINE SolveTrilinos2( n, Xvec, RHSVec, Rounds, TOL, &
      verbosity, triliContainer, ierr) &
                BIND(C,name='SolveTrilinos2')
        USE, INTRINSIC :: iso_c_binding
        INTEGER(KIND=c_int) :: n, Rounds, verbosity
        REAL(KIND=c_double) :: Xvec(n),RHSvec(n),TOL
        INTEGER(KIND=C_INTPTR_T) :: triliContainer
        INTEGER(KIND=C_INT) :: ierr
      END SUBROUTINE SolveTrilinos2

      !! note: SolveTrilinos3 does not yet exist, it would update the matrix
      !!       but not the preconditioner.

      !! destroy the data structures (should be called when the matrix has
      !! to be updated and SolveTrilinos1 has to be called again).
      SUBROUTINE SolveTrilinos4(triliContainer) BIND(C,name='SolveTrilinos4')
        USE, INTRINSIC :: iso_c_binding
        INTEGER(KIND=C_INTPTR_T) :: triliContainer
      END SUBROUTINE SolveTrilinos4

    END INTERFACE
#endif

  CHARACTER(*), PARAMETER :: Caller = 'SParIterSolver' 
    
  !******************************************************************
  SaveGlobalData => GlobalData
  GlobalData     => SParMatrixDesc
  ParEnv         => GlobalData % ParEnv
  ParEnv % ActiveComm = SourceMatrix % Comm
  SplittedMatrix => SParMatrixDesc % SplittedMatrix

  Rows => SourceMatrix % Rows
  Cols => SourceMatrix % Cols
  Vals => SourceMatrix % Values

  !******************************************************************

  CALL Info(Caller,'Solving linear in parallel with iterative methods',Level=8)

  Params => Solver % Values

  IF (ListGetLogical(Params,'Linear System Use HYPRE', Found )) THEN
#ifdef HAVE_HYPRE
    CALL SolveHypre(SourceMatrix,XVec,RHSVec,Solver,&
        ParallelInfo,SplittedMatrix)    
    RETURN
#else
    CALL Fatal(Caller,'This version has been compiled without HYPRE!')
#endif
  END IF

  IF (ListGetLogical( Params,'Linear System Use Trilinos', Found )) THEN
#ifdef HAVE_TRILINOS
    CALL SolveTrilinos()      
    RETURN
#else
    CALL Fatal(Caller,'This version has been compiled without Trilinos!')
#endif
  END IF


  NeedILU = .FALSE.
  DO i=1,ParEnv % PEs
    IF ( SplittedMatrix % IfMatrix(i) % NumberOFRows>0 ) THEN
      NeedILU = ALLOCATED( SplittedMatrix % IfMatrix(i) % ILUValues )
      EXIT
    END IF
  END DO
  NeedMass = ASSOCIATED( SplittedMatrix % InsideMatrix % MassValues )
  NeedDamp = ASSOCIATED( SplittedMatrix % InsideMatrix % DampValues )
  NeedPrec = ASSOCIATED( SplittedMatrix % InsideMatrix % PrecValues )

  CALL ZeroSplittedMatrix( SplittedMatrix )

  !------------------------------------------------------------------
  !
  ! Copy the Matrix % Values into SplittedMatrix
  !
  !------------------------------------------------------------------

  CALL Info(Caller,'Copying Matrix values into SplittedMatrix',Level=20)
  CALL ResetTimer('SplittedMatrix')

  
  GT => SplittedMatrix % GlueTable
  DO i = 1, SourceMatrix % NumberOfRows
     GRow = ParallelInfo % GlobalDOFs(i)

     DO j = SourceMatrix % Rows(i),SourceMatrix % Rows(i+1) - 1

        GCol = ParallelInfo % GlobalDOFs(SourceMatrix % Cols(j))
        
        IF ( GT % Inds(j) > 0 ) THEN
           SplittedMatrix % InsideMatrix % Values( GT % Inds(j) ) = &
                SplittedMatrix % InsideMatrix % Values( &
                   GT % Inds(j) ) + SourceMatrix % Values(j)

           IF ( NeedPrec ) &
              SplittedMatrix % InsideMatrix % PrecValues( GT % Inds(j) ) = &
                   SplittedMatrix % InsideMatrix % PrecValues( &
                      GT % Inds(j) ) + SourceMatrix % PrecValues(j)
           IF ( NeedMass ) &
              SplittedMatrix % InsideMatrix % MassValues( GT % Inds(j) ) = &
                   SplittedMatrix % InsideMatrix % MassValues( &
                      GT % Inds(j) ) + SourceMatrix % MassValues(j)
           IF ( NeedDamp ) &
              SplittedMatrix % InsideMatrix % DampValues( GT % Inds(j) ) = &
                   SplittedMatrix % InsideMatrix % DampValues( &
                      GT % Inds(j) ) + SourceMatrix % DampValues(j)
           IF ( NeedILU ) &
              SplittedMatrix % InsideMatrix % ILUValues( GT % Inds(j) ) = &
                   SplittedMatrix % InsideMatrix % ILUValues( &
                      GT % Inds(j) ) + SourceMatrix % ILUValues(j)

        ELSE IF ( (GT % Inds(j) + ParEnv % PEs) >= 0 ) THEN

           ifind = ABS(GT % Inds(j))
           CurrIf => SplittedMatrix % IfMatrix(ifind)

           RowInd = -1
           IF ( CurrIf % NumberOfRows > 0 ) THEN
              RowInd = SearchIAItem( CurrIf % NumberOfRows, CurrIf % GRows, GRow )
           END IF

           IF ( RowInd /= -1 ) THEN
              DO l = CurrIf % Rows(rowind), CurrIf % Rows(rowind+1)-1

                 IF ( GCol == CurrIf % Cols(l) ) THEN
                    CurrIf % Values(l) = CurrIf % Values(l) + &
                         SourceMatrix % Values(j)
                    IF ( NeedPrec ) &
                       CurrIf % PrecValues(l) = CurrIf % PrecValues(l) + &
                            SourceMatrix % PrecValues(j)
                    IF ( NeedMass ) &
                       CurrIf % MassValues(l) = CurrIf % MassValues(l) + &
                            SourceMatrix % MassValues(j)
                    IF ( NeedDamp ) &
                       CurrIf % DampValues(l) = CurrIf % DampValues(l) + &
                            SourceMatrix % DampValues(j)
                    IF ( NeedILU ) &
                       CurrIf % ILUValues(l) = CurrIf % ILUValues(l) + &
                            SourceMatrix % ILUValues(j)
                    EXIT
                 END IF
              END DO
           END IF

        ELSE IF ((GT % Inds(j) + (2*ParEnv % PEs)) >= 0) THEN

           ifind = -ParEnv % PEs + ABS(GT % Inds(j))
           CurrIf => SplittedMatrix % NbsIfMatrix(ifind)
                
           RowInd = -1
           IF ( CurrIf % NumberOfRows > 0 ) THEN
              RowInd = SearchIAItem( CurrIf % NumberOfRows, CurrIf % GRows, GRow )
           END IF

           IF ( RowInd /= -1 ) THEN
              DO l = CurrIf % Rows(RowInd), CurrIf % Rows(RowInd+1)-1
                 IF ( GCol == CurrIf % Cols(l) ) THEN
                    CurrIf % Values(l) = CurrIf % Values(l) + &
                         SourceMatrix % Values(j)
                    IF ( NeedPrec ) &
                       CurrIf % PrecValues(l) = CurrIf % PrecValues(l) + &
                            SourceMatrix % PrecValues(j)
                    IF ( NeedMass ) &
                       CurrIf % MassValues(l) = CurrIf % MassValues(l) + &
                            SourceMatrix % MassValues(j)
                    IF ( NeedDamp ) &
                       CurrIf % DampValues(l) = CurrIf % DampValues(l) + &
                            SourceMatrix % DampValues(j)
                    IF ( NeedILU ) &
                       CurrIf % ILUValues(l) = CurrIf % ILUValues(l) + &
                            SourceMatrix % ILUValues(j)
                    EXIT
                 END IF
              END DO
           END IF
        END IF
     END DO
  END DO

  CALL GlueFinalize( SourceMatrix, SplittedMatrix, ParallelInfo )
  CALL CheckTimer('SplittedMatrix',Level=7,Delete=.TRUE.)


  !------------------------------------------------------------------
  ! Call the actual solver routine (based on older design)
  !------------------------------------------------------------------
  CALL Info(Caller,'Going into actual parallel solution',Level=20)

  CALL SolveHutiter( SourceMatrix, SParMatrixDesc % SplittedMatrix, &
      ParallelInfo, RHSVec, XVec, Solver, Errinfo )

  GlobalData => SaveGlobalData
  

CONTAINS

  
#ifdef HAVE_TRILINOS
  SUBROUTINE SolveTrilinos()

    ! we attempt to read Trilinos settings from an XML file,
    ! which is their usual way of getting parameters. If no 
    ! file is given, we use default settings and issue a    
    ! warning.
    xmlfile = ListGetString( Params, 'Trilinos Parameter File', Found )
    IF (.NOT. Found) THEN
      xmlfile = 'none'
    END IF
    xmlfile = TRIM(xmlfile)//C_NULL_CHAR

    ! tolerance and max iter are taken from the Elmer
    ! internal list to overrule the settings in the XML file
    TOL = ListGetConstReal( Params, &
        'Linear System Convergence Tolerance', Found )
    IF ( .NOT. Found ) TOL=-1.0

    Rounds = ListGetInteger( Params, &
        'Linear System Max Iterations', Found )
    IF ( .NOT. Found ) Rounds=-1


    ! I think nrows == ncols here?
    n = SourceMatrix%NumberOfRows ! number of nodes of local elements
    nnz = SourceMatrix%Rows(n+1) ! number of local nonzeros

    verbosity=0 !TODO: get this from the simulation verbosity or something
    NewSetup=ListGetLogical( Params, 'Linear System Refactorize',Found ) 
    IF (NewSetup) THEN
      IF (SourceMatrix % Trilinos/=0) THEN
        CALL SolveTrilinos4(SourceMatrix % Trilinos)
      END IF
    END IF
    ! setup solver/preconditioner
    IF (SourceMatrix % Trilinos==0) THEN

      ALLOCATE( Owner(n))
      Owner = 0
      DO i=1,n
        IF (ParallelInfo % NeighbourList(i) % Neighbours(1)== ParEnv % MyPE) THEN
          Owner(i) = 1
        END IF
      END DO

      CALL SolveTrilinos1(n, nnz, & 
          SourceMatrix % Rows, &
          SourceMatrix % Cols, &
          SourceMatrix % Values, &
          ParallelInfo%GlobalDofs, Owner, &
          xmlfile, verbosity, SourceMatrix % Trilinos, &
          CurrentModel%Nodes%NumberOfNodes, &
          CurrentModel%Nodes%x, CurrentModel%Nodes%y, CurrentModel%Nodes%z, ierr)     

      DEALLOCATE( Owner )
      IF (ierr<0) THEN
        CALL Fatal(Caller,'Failed to construct Trilinos solver')
      ELSE IF (ierr>0) THEN
        CALL Warn(Caller,'Warning issued when trying to construct Trilinos solver')          
      END IF
    END IF
    ! solve using previously computed Trilinos data structures.
    ! NOTE: this is only correct if the matrix has not changed,
    ! otherwise we should use the SolveTrilinos3 function, which is
    ! not implemented, yet. This function will not update the matrix
    ! in the solver and thus solve an old system if A has changed. 
    CALL SolveTrilinos2( n, Xvec, RHSvec, &
        Rounds, TOL, verbosity, SourceMatrix % Trilinos, ierr)

    IF (ierr<0) THEN
      CALL Fatal(Caller,&
          'Linear system solve using Trilinos caused an error');
    ELSE IF (ierr>0) THEN
      CALL NumericalError(Caller,&
          'Linear system solve using Trilinos issued a warning')
    END IF

    ALLOCATE( VecEPerNB( ParEnv % PEs ) )
    VecEPerNB = 0
    DO i = 1, SourceMatrix % NumberOfRows
      IF ( SIZE(ParallelInfo % NeighbourList(i) % Neighbours) > 1 ) THEN
        IF ( ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
          DO j = 1, SIZE(ParallelInfo % NeighbourList(i) % Neighbours)
            IF (ParallelInfo % NeighbourList(i) % Neighbours(j)/=ParEnv % MyPE) THEN
              nbind = ParallelInfo % NeighbourList(i) % Neighbours(j) + 1
              VecEPerNB(nbind) = VecEPerNB(nbind) + 1

              SplittedMatrix % ResBuf(nbind) % ResVal(VecEPerNB(nbind)) = XVec(i)
              SplittedMatrix % ResBuf(nbind) % ResInd(VecEPerNB(nbind)) = &
                  ParallelInfo % GlobalDOFs(i)
            END IF
          END DO
        END IF
      END IF
    END DO

    CALL ExchangeResult( SourceMatrix, SplittedMatrix, ParallelInfo, XVec )
    DEALLOCATE( VecEPerNB )

    !     CALL ExchangeSourceVec( SourceMatrix, SplittedMatrix, ParallelInfo, RHSVec )

  END SUBROUTINE SolveTrilinos
#endif
   
  
!*********************************************************************
END SUBROUTINE SParIterSolver
!*********************************************************************



SUBROUTINE SolveHypre(Matrix, XVec, RHSVec, Solver, ParallelInfo, SplittedMatrix )

  USE, INTRINSIC :: iso_c_binding                

  TYPE (Matrix_t) :: Matrix
  REAL(KIND=dp), DIMENSION(:) :: XVec, RHSVec
  TYPE (Solver_t) :: Solver
  TYPE (ParallelInfo_t), OPTIONAL :: ParallelInfo
  TYPE (SplittedMatrixT), POINTER, OPTIONAL :: SplittedMatrix


  CHARACTER(*), PARAMETER :: Caller = 'HypreSolver' 

#ifndef HAVE_HYPRE
  CALL Fatal(Caller,'Serial Hypre requested but the library is not linked in!')
#else
    
  ! Local variables
  LOGICAL :: Parallel
  INTEGER :: i, j, k, l, n
  REAL(KIND=dp) :: TOL, hypre_dppara(10) = 0
  INTEGER :: ILUn, BILU, Rounds, buf(2), src, status(MPI_STATUS_SIZE), ssz,nob, &
      hypremethod,  hypre_intpara(20) = 0
  INTEGER, DIMENSION(:), ALLOCATABLE :: VecEPerNB

  INTEGER, ALLOCATABLE :: Owner(:), Aperm(:)
  REAL(KIND=dp), TARGET :: DummyPrecVals(1)
  REAL(KIND=dp), POINTER :: Vals(:)
  INTEGER, POINTER :: Rows(:), Cols(:)

  LOGICAL :: Found, NewSetup, UpdateTolerance, DoAMS
  INTEGER :: verbosity, myverb
  INTEGER :: nrows, ncols, nnz
  TYPE(ValueList_t), POINTER :: Params

  TYPE(Matrix_t), POINTER :: GM
  INTEGER:: nnd,ind(2), precond
  REAL(KIND=dp), POINTER :: PrecVals(:)
  REAL(KIND=dp), ALLOCATABLE :: xx_d(:),yy_d(:),zz_d(:)  
  INTEGER, ALLOCATABLE :: nodeowner(:),nodeperm(:),bperm(:), bowner(:)

  
  INTERFACE
    !! create HYPRE matrix and setup solver/preconditioner
    SUBROUTINE SolveHYPRE1( n, Rows, Cols, Vals, Precond, PrecVals, GDOFs, &
        Owner, ILUn, BILU, hypremethod, hypre_intpara, hypre_dppara, rounds, &
        TOL, verbosity, hypreContainer, fcomm ) BIND(C,name="solvehypre1")
      USE, INTRINSIC :: iso_c_binding
      INTEGER(KIND=c_int) :: n, Rows(n), Cols(n), GDOFs(n), &
          PE, Owner(n), ILUn, BILU, hypremethod, precond, &
          symmetry, maxlevel, hypre_intpara(20), verbosity,fcomm,rounds
      REAL(KIND=c_double) :: Vals(n),Xvec(n),RHSvec(n),TOL,threshold,filter, &
          hypre_dppara(10), PrecVals(*)
      INTEGER(KIND=C_INTPTR_T) :: hypreContainer
    END SUBROUTINE SolveHYPRE1

    !! solve linear system with same matrix as in SolveHYPRE1
    SUBROUTINE SolveHYPRE2( n, GDOFs, &
        Owner, Xvec, RHSVec, Rounds, TOL, verbosity, hypreContainer, fcomm) BIND(C,name="solvehypre2")
      USE, INTRINSIC :: iso_c_binding
      INTEGER(KIND=c_int) :: n, GDOFs(n), Owner(n), Rounds, verbosity, fcomm
      REAL(KIND=c_double) :: Xvec(n),RHSvec(n),TOL
      INTEGER(KIND=C_INTPTR_T) :: hypreContainer
    END SUBROUTINE SolveHYPRE2

    !! note: SolveHYPRE3 does not yet exist, it would update the matrix
    !!       but not the preconditioner.

    !! destroy the data structures (should be called when the matrix has
    !! to be updated and SolveHYPRE1 has to be called again).
    SUBROUTINE SolveHYPRE4(hypreContainer, verbosity ) BIND(C,name="solvehypre4")
      USE, INTRINSIC :: iso_c_binding
      INTEGER(KIND=C_INTPTR_T) :: hypreContainer
      INTEGER(KIND=c_int) :: verbosity
    END SUBROUTINE SolveHYPRE4
    
    SUBROUTINE CreateHypreAMS(nrows,rows,cols,vals,n,grows,gcols,gvals, &
        perm, invperm, globaldofs, owner, Bperm,nodeowner,xvec, rhsvec, pe, ILUn, rounds, &
        TOL, xx_d, yy_d, zz_d, hypremethod, hypre_intpara, hypre_dppara,verbosity,hyprecontainer,fcomm ) & 
        BIND(C,name="createhypreams")
      
      USE, INTRINSIC :: iso_c_binding
      INTEGER(KIND=c_int) :: nrows, n, Rows(*), Cols(*), Perm(*), INVPerm(*), &
          Grows(*), gcols(*), PE, Owner(*), Rounds, ILUn, hypremethod, fcomm, &
          symmetry, maxlevel, hypre_intpara(20),bperm(*),nodeOwner(*),globaldofs(*),verbosity
      REAL(KIND=c_double) :: Vals(*),Xvec(*),RHSvec(*),TOL,threshold,filter, &
          hypre_dppara(10), Gvals(*), xx_d(*), yy_d(*), zz_d(*)
      INTEGER(KIND=C_INTPTR_T) :: hypreContainer
    END SUBROUTINE CreateHypreAMS

    SUBROUTINE UpdateHypre(TOL, hypremethod, hypreContainer, verbosity, fcomm ) BIND(C,name="updatehypre")
      USE, INTRINSIC :: iso_c_binding
      REAL(KIND=c_double) :: TOL
      INTEGER(KIND=c_int) :: hypremethod
      INTEGER(KIND=C_INTPTR_T) :: hypreContainer
      INTEGER(KIND=c_int) :: verbosity, fcomm
    END SUBROUTINE UpdateHypre
    
  END INTERFACE

  
  CALL Info(Caller,'Solving linear system using HYPRE library',Level=6)

  Rows => Matrix % Rows
  Cols => Matrix % Cols
  Vals => Matrix % Values
  Params => Solver % Values
  
  CALL SetHypreParameters(Params, hypremethod, ilun, hypre_dppara, hypre_intpara )

  TOL = ListGetCReal( Params,'Linear System Convergence Tolerance', Found )
  IF ( .NOT. Found ) TOL = 1.0d-6
  
  Rounds = ListGetInteger( Params,'Linear System Max Iterations', Found )
  IF ( .NOT. Found ) Rounds = 1000  
  
  ! Hypre wants to have a continuous ascending numbering across
  ! partitions, try creating such a beast:
  ! ------------------------------------------------------------
  Parallel = PRESENT(ParallelInfo) 
  IF( Parallel ) THEN
    n = SIZE(ParallelInfo % GlobalDOFs)
    ALLOCATE( Owner(n), Aperm(n) )
    CALL ContinuousNumbering(ParallelInfo,Matrix % Perm,APerm,Owner)
    ! Newer hypre libraries require zero based indexing    
    Aperm = Aperm-1
  ELSE
    n = Matrix % NumberOfRows
    ALLOCATE( Owner(n), Aperm(n) )
    DO i=1,n
      Aperm(i) = i-1
    END DO
    Owner = 1
  END IF

  CALL SParIterActiveBarrier()

  
  !------------------------------------------------------------
  verbosity = ListGetInteger( CurrentModel % Simulation,'Max Output Level',Found )
  IF( .NOT. Found ) verbosity = 10

  NewSetup = ListGetLogical( Params, 'Linear System Refactorize',Found ) 
  IF(.NOT.Found) NewSetup = .TRUE.
  
  IF (ListGetLogical(Params, 'HYPRE Block Diagonal', Found)) THEN
    bilu = Solver % Variable % Dofs
  ELSE
    bilu = 1
  END IF

  ! AMS requites additional information and is therefore treated separately.
  DoAMS = ( ( MODULO(hypremethod,100) == 2 ) .OR. (hypremethod/100 == 2) )
  
  IF (NewSetup) THEN
    IF (Matrix % Hypre /= 0) THEN
      CALL Info(Caller,'Destroy old Hypre solver structures',Level=10)
      myverb = verbosity
      IF(ParEnv % MyPe /= 0) myverb = 0
      CALL SolveHYPRE4(Matrix % Hypre, myverb )
      Matrix % Hypre = 0
    END IF
  END IF

  IF (Matrix % Hypre == 0) THEN
    CALL Info(Caller,'Setting up Hypre solver structures',Level=10)
    PrecVals => Matrix % PrecValues
    IF(ASSOCIATED(PrecVals)) THEN
      precond=1
    ELSE
      precond=0
      PrecVals => DummyPrecVals
    END IF

    IF( DoAMS ) THEN
      IF(Matrix % NumberOfRows > Solver % Mesh % NumberOfEdges ) THEN
!       CALL Fatal(Caller,'More unknowns that edges, current Hypre AMS can not be used!')
      END IF

      CALL PrepareHypreAMS() 
      nnd = Solver % Mesh % NumberOfNodes      
      CALL CreateHYPREAMS( Matrix % NumberOfRows, Rows, Cols, Vals, &
          nnd,GM % Rows,GM % Cols,GM % Values,Aperm,Aperm,Aperm,Owner, &
          Bperm,NodeOwner,Xvec,RHSvec,ParEnv % myPE, ILUn, Rounds,TOL,  &
          xx_d,yy_d,zz_d,hypremethod,hypre_intpara, hypre_dppara,verbosity, &
          Matrix % Hypre, Matrix % Comm)
      CALL CleanHypreAMS() 
    END IF
    CALL SolveHYPRE1( Matrix % NumberOfRows, Rows, Cols, Vals, Precond, &
        PrecVals, Aperm, Owner,  ILUn, BILU, hypremethod, hypre_intpara, hypre_dppara,&
        rounds, TOL, verbosity, Matrix % Hypre, Matrix % Comm)
  END IF


  ! In some cases the stopping tolerance is adapted during the solution procedure.
  ! Make an update if needed:
  UpdateTolerance = ListGetLogical( Params, 'Linear System Adaptive Tolerance', Found )
  IF (UpdateTolerance) THEN
    WRITE(Message,'(A,ES12.3)') 'Updating Hypre tolerances to:',TOL
    CALL Info(Caller,Message,Level=10)
    CALL UpdateHypre( TOL, hypremethod, Matrix % Hypre, verbosity, Matrix % Comm )
  END IF

  ! solve using previously computed HYPRE data structures.
  ! NOTE: this is only correct if the matrix has not changed,
  ! otherwise we should use the SolveHYPRE3 function, which is
  ! not implemented, yet. This function will not update the matrix
  ! in the solver and thus solve an old system if A has changed. 
  !-----------------------------------------------------------------------
  CALL Info(Caller,'Solving previously created Hypre setup',Level=10)
  CALL SolveHYPRE2( Matrix % NumberOfRows, Aperm, Owner, Xvec, RHSvec, &
      Rounds, TOL, verbosity, Matrix % Hypre, Matrix % Comm )

  IF(Parallel) CALL ExchangeHypreResults()    
  
  CALL SParIterActiveBarrier()
  DEALLOCATE( Owner, Aperm )

  CALL Info(Caller,'Finished solving linear system with HYPRE',Level=12)


CONTAINS

  SUBROUTINE PrepareHypreAMS()

    TYPE(Mesh_t), POINTER :: Mesh

    Mesh => Solver % Mesh
    
    nnd = Mesh % NumberOfNodes
    ALLOCATE( NodeOwner(nnd), NodePerm(nnd), Bperm(nnd))

    DO i=1,nnd
      NodePerm(i)=i
    END DO

    IF( Parallel ) THEN
      NodeOwner = 0
      CALL ContinuousNumbering( Mesh % ParallelInfo, &
          NodePerm, BPerm, NodeOwner, nnd, Mesh)
      bPerm = bPerm-1 
    ELSE
      NodeOwner = 1
      DO i=1,nnd
        bperm(i) = i-1
      END DO
    END IF

    GM => AllocateMatrix()
    GM % FORMAT = MATRIX_LIST

    DO i=Mesh % NumberofEdges,1,-1
      ind = Mesh % Edges(i) % NodeIndexes
      IF( Parallel ) THEN
        IF (Mesh % ParallelInfo % GlobalDOFs(ind(1))> &
            Mesh % ParallelInfo % GlobalDOFs(ind(2))) THEN
          k=ind(1); ind(1)=ind(2);ind(2)=k
        END IF
      ELSE
        IF (ind(1)> ind(2)) THEN
          k=ind(1); ind(1)=ind(2);ind(2)=k
        END IF
      END IF        
      IF(Matrix % Complex) THEN
        CALL List_AddToMatrixElement(gm % listmatrix,2*i-1,ind(1),-1._dp)
        CALL List_AddToMatrixElement(gm % listmatrix,2*i-1,ind(2), 1._dp)
        CALL List_AddToMatrixElement(gm % listmatrix,2*i,ind(1),-1._dp)
        CALL List_AddToMatrixElement(gm % listmatrix,2*i,ind(2), 1._dp)
      ELSE
        CALL List_AddToMatrixElement(gm % listmatrix,i,ind(1),-1._dp)
        CALL List_AddToMatrixElement(gm % listmatrix,i,ind(2), 1._dp)
      END IF
    END DO
    CALL List_tocrsMatrix(gm)
    
    nnd = Mesh % NumberOfEdges
    IF(Matrix % Complex) nnd=nnd*2
    ALLOCATE(xx_d(nnd),yy_d(nnd),zz_d(nnd) )
    CALL CRS_MatrixVectorMultiply(gm,Mesh % Nodes % x,xx_d)
    CALL CRS_MatrixVectorMultiply(gm,Mesh % Nodes % y,yy_d)
    CALL CRS_MatrixVectorMultiply(gm,Mesh % Nodes % z,zz_d)
    
    nnd = Mesh % NumberOfNodes
  END SUBROUTINE PrepareHypreAMS


  SUBROUTINE CleanHypreAMS() 
    
    IF(.NOT. ASSOCIATED(GM)) THEN
      CALL Fatal(Caller,'Matrix "GM" should be allocated!')
    END IF
    
    DEALLOCATE(GM % Rows) 
    DEALLOCATE(GM % Cols) 
    DEALLOCATE(GM % Diag) 
    DEALLOCATE(GM % Values) 
    DEALLOCATE(GM)

  END SUBROUTINE CleanHypreAMS


  SUBROUTINE ExchangeHypreResults()
    INTEGER :: nbind
    
    ALLOCATE( VecEPerNB( ParEnv % PEs ) )
    VecEPerNB = 0
    DO i = 1, Matrix % NumberOfRows
      IF ( SIZE(ParallelInfo % NeighbourList(i) % Neighbours) > 1 ) THEN
        IF ( ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
          DO j = 1, SIZE(ParallelInfo % NeighbourList(i) % Neighbours)
            IF (ParallelInfo % NeighbourList(i) % Neighbours(j)/=ParEnv % MyPE) THEN
              nbind = ParallelInfo % NeighbourList(i) % Neighbours(j) + 1
              VecEPerNB(nbind) = VecEPerNB(nbind) + 1

              SplittedMatrix % ResBuf(nbind) % ResVal(VecEPerNB(nbind)) = XVec(i)
              SplittedMatrix % ResBuf(nbind) % ResInd(VecEPerNB(nbind)) = &
                  ParallelInfo % GlobalDOFs(i)
            END IF
          END DO
        END IF
      END IF
    END DO

    CALL ExchangeResult( Matrix, SplittedMatrix, ParallelInfo, XVec )
    DEALLOCATE( VecEPerNB )

  END SUBROUTINE ExchangeHypreResults


#endif
  
END SUBROUTINE SolveHypre




!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
SUBROUTINE SolveHutiter( SourceMatrix, SplittedMatrix, ParallelInfo, &
    RHSVec, XVec, Solver, ErrInfo )

  USE IterSolve

  TYPE (SplittedMatrixT), POINTER :: SplittedMatrix
  TYPE(Matrix_t), TARGET :: SourceMatrix
  TYPE (ParallelInfo_t) :: ParallelInfo
  TYPE (Solver_t) :: Solver
  TYPE (ErrInfoT) :: ErrInfo
  REAL(KIND=dp), DIMENSION(:) :: RHSVec, XVec

  ! Local variables

  LOGICAL :: stat, EdgeBasis = .FALSE.
  INTEGER :: i, j, k, l, nbind, dof
  INTEGER, DIMENSION(:), ALLOCATABLE :: VecEPerNB
  REAL(KIND=dp) :: dpar(HUTI_DPAR_DFLTSIZE)
  INTEGER :: ipar(HUTI_IPAR_DFLTSIZE)
  REAL(KIND=dp), DIMENSION(:,:), POINTER :: Work
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: TmpXVec, TmpRHSVec, r
  INTEGER(KIND=AddrInt) :: AddrFunc
  EXTERNAL :: AddrFunc
  REAL(KIND=dp) :: ILUT_TOL
  INTEGER :: ILUn

  TYPE(Matrix_t), POINTER :: CM,SaveMatrix
  INTEGER, POINTER :: SPerm(:)
  INTEGER, POINTER CONTIG :: SCols(:)
  REAL(kind=dp)::tt,rt

  !*******************************************************************

  PIGpntr => GlobalData

  !----------------------------------------------------------------------
  ! Initialize Right-Hand-Side
  !----------------------------------------------------------------------
  ALLOCATE(TmpRHSVec(SplittedMatrix % InsideMatrix % NumberOfRows))
  TmpRHSVec = 0

  ALLOCATE(r(SIZE(rhsvec)))
  r = RHSvec
  CALL ExchangeRHSIf( SourceMatrix, SplittedMatrix, &
         ParallelInfo, r, TmpRHSVec )
  DEALLOCATE(r)

  !----------------------------------------------------------------------
  !
  ! Initialize temporary xvec for iterator. The original XVec contains
  ! also the items on interfaces. Initialize also global pointer.
  !
  !----------------------------------------------------------------------

!  EdgeBasis = .NOT.ListGetLogicalAnyBC( CurrentModel, 'Stride Projector' )
!  EdgeBasis = EdgeBasis.AND.ListGetLogical( Solver % Values, 'Edge Basis', stat )
!  IF ( EdgeBasis ) THEN
!    IF( ASSOCIATED( SplittedMatrix % InsideMatrix % Eperm) ) &
!      DEALLOCATE( SplittedMatrix % InsideMatrix % EPerm )
!
!    ALLOCATE( SplittedMatrix % InsideMatrix % EPerm(SplittedMatrix % InsideMatrix % NumberOfRows) )
!    SplittedMatrix % InsideMatrix % EPerm = 0
!  END IF
 
 
  ALLOCATE( TmpXVec(SplittedMatrix % InsideMatrix % NumberOfRows) )
  j = 0
  SplittedMatrix % InsideMatrix % ExtraDOFs=0
  DO i = 1, SourceMatrix % NumberOfRows
    IF ( ParallelInfo % NeighbourList(i) % Neighbours(1)==ParEnv % MyPE ) THEN
      j = j + 1
      TmpXVec(j) = XVec(i)

      IF(i>SourceMatrix % NumberOfRows-SourceMatrix % ExtraDOFs) THEN
        SplittedMatrix % InsideMatrix % ExtraDOFs = &
           SplittedMatrix % InsideMatrix % ExtraDOFs+1
      END IF

      IF( EdgeBasis) THEN
        k = SourceMatrix % InvPerm(i)
        IF( k>Solver % Variable % DOFs*Solver % Mesh % NumberOfNodes ) &
             SplittedMatrix % InsideMatrix % EPerm(j)=1
      END IF
    END IF
  END DO

  SaveMatrix => GlobalMatrix
  GlobalMatrix => SplittedMatrix % InsideMatrix
  GlobalData % SplittedMatrix => SplittedMatrix
  GlobalMatrix % MatVecSubr = SourceMatrix % MatVecSubr

  GlobalMatrix % Ematrix => SourceMatrix
  GlobalMatrix % COMPLEX = SourceMatrix % COMPLEX

 !----------------------------------------------------------------------
 ! Set up the preconditioner
 !----------------------------------------------------------------------
 IF (SplittedMatrix % InsideMatrix % NumberOFRows>0) THEN
!  IF (SplittedMatrix % InsideMatrix % Diag(1)==0) THEN
     DO i = 1, SplittedMatrix % InsideMatrix % NumberOfRows
       DO j = SplittedMatrix % InsideMatrix % Rows(i), &
            SplittedMatrix % InsideMatrix % Rows(i+1) - 1
          IF ( SplittedMatrix % InsideMatrix % Cols(j) == i ) THEN
             SplittedMatrix % InsideMatrix % Diag(i) = j
             EXIT
          END IF
        END DO
      END DO
!   END IF
  END IF

#if 0
  DO i = 1, SplittedMatrix % InsideMatrix % NumberOfRows
    j = SplittedMatrix % InsideMatrix % Diag(i)
    IF ( SplittedMatrix % InsideMatrix % Values(j)==0 ) THEN
      l=0
      DO k=SplittedMatrix % InsideMatrix % Rows(i),SplittedMatrix % InsideMatrix % Rows(i+1)-1
        IF ( SplittedMatrix % InsideMatrix % Values(k) /= 0 ) l = 1
      END DO
      IF ( l == 0 ) THEN
         SplittedMatrix % InsideMatrix % Values(j) = 1
         TmpRHSVec(i) = 0; TMPXVec(i) = 0
      END IF
    END IF
  END DO
#endif


 !----------------------------------------------------------------------
 ! Call the main iterator routine
 !----------------------------------------------------------------------
  CM => SourceMatrix % ConstraintMatrix
  IF (ASSOCIATED(CM)) THEN
    ALLOCATE(SPerm(SourceMatrix % NumberOfRows)); SPerm=0
    j = 0
    DO i = 1, SourceMatrix % NumberOfRows
      IF (ParallelInfo % NeighbourList(i) % Neighbours(1)==ParEnv % MyPE) THEN
         j = j + 1
         SPerm(i)=j
      END IF
    END DO
    SCols => CM % Cols
    ALLOCATE(CM % Cols(SIZE(SCols)))
    CM % Cols = SPerm(SCols)
    SplittedMatrix % InsideMatrix % Perm => SPerm
    SplittedMatrix % InsideMatrix % ConstraintMatrix => CM
  END IF

  IF ( .NOT. SourceMatrix % COMPLEX ) THEN
     CALL IterSolver( SplittedMatrix % InsideMatrix, TmpXVec, &
        TmpRHSVec, Solver, DotF=AddrFunc(SParDotProd), NormF=AddrFunc(SParNorm), &
          matVecF=AddrFunc(SParMatrixVector) )
  ELSE
     CALL IterSolver( SplittedMatrix % InsideMatrix, TmpXVec, &
       TmpRHSVec, Solver, DotF=AddrFunc(SParCDotProd), NormF=AddrFunc(SParCNorm), &
         MatVecF=AddrFunc(SParCMatrixVector) )
  END IF

  IF (ASSOCIATED(CM)) THEN
    DEALLOCATE(CM % Cols)
    CM % Cols => SCols
    DEALLOCATE(SPerm)
    SplittedMatrix % InsideMatrix % Perm => NULL()
    SplittedMatrix % InsideMatrix % ConstraintMatrix => NULL()
  END IF

  GlobalMatrix => SaveMatrix

  !----------------------------------------------------------------------
  ! Collect the result
  !----------------------------------------------------------------------
  ALLOCATE( VecEPerNB( ParEnv % PEs ) )
  VecEPerNB = 0

  j = 0
  DO i = 1, SourceMatrix % NumberOfRows
     IF ( ParallelInfo % NeighbourList(i) % Neighbours(1)==ParEnv % MyPE ) THEN
        j = j + 1
        XVec(i) = TmpXVec(j)
     END IF
  END DO

  DO i = 1, SourceMatrix % NumberOfRows
     IF ( SIZE(ParallelInfo % NeighbourList(i) % Neighbours) > 1 ) THEN
        IF ( ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
           DO j = 1, SIZE(ParallelInfo % NeighbourList(i) % Neighbours)
              IF (ParallelInfo % NeighbourList(i) % Neighbours(j)/=ParEnv % MyPE) THEN
                 nbind = ParallelInfo % NeighbourList(i) % Neighbours(j) + 1
                 VecEPerNB(nbind) = VecEPerNB(nbind) + 1

                 SplittedMatrix % ResBuf(nbind) % ResVal(VecEPerNB(nbind)) = XVec(i)
                 SplittedMatrix % ResBuf(nbind) % ResInd(VecEPerNB(nbind)) = &
                   ParallelInfo % GlobalDOFs(i)
              END IF
           END DO
        END IF
     END IF
  END DO

  CALL ExchangeResult( SourceMatrix,SplittedMatrix,ParallelInfo,XVec )

  !----------------------------------------------------------------------
  ! Clean the work space
  !----------------------------------------------------------------------
  DEALLOCATE( TmpXVec, TmpRHSVec, VecEPerNB )
!----------------------------------------------------------------------
END SUBROUTINE SolveHutiter
!----------------------------------------------------------------------



!---------------------------------------------------------------------
!> External Matrix - Vector operations (Parallel real, version)
!> Multiply vector u with the global matrix,  return the result in v
!> Called from HUTIter library
!----------------------------------------------------------------------
  SUBROUTINE SParMatrixVector( u,v,ipar )
!----------------------------------------------------------------------

  IMPLICIT NONE

  ! Input parameters

  INTEGER, DIMENSION(*) :: ipar
  REAL(KIND=dp), DIMENSION(*) :: u, v

  ! Local parameters

  INTEGER :: i, j, k, l, n, ii, colind, rowind
  TYPE( IfVecT), POINTER :: IfV
  TYPE( IfLColsT), POINTER :: IfL, IfO
  TYPE (BasicMatrix_t), POINTER :: CurrIf
  TYPE (Matrix_t), POINTER :: InsideMatrix

  INTEGER :: nneigh
  TYPE(Buff_t), POINTER ::  buffer(:)
  REAL(KIND=dp) :: rsum
  REAL(KIND=dp), POINTER CONTIG :: Vals(:)
  INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
  INTEGER, ALLOCATABLE :: neigh(:), recv_size(:), requests(:)
  REAL(kind=dp) :: s
  !*******************************************************************

  InsideMatrix => GlobalData % SplittedMatrix % InsideMatrix
  n = InsideMatrix % NumberOfRows

  !----------------------------------------------------------------------
  !
  ! Compute the interface contribution for each neighbour
  !
  !----------------------------------------------------------------------

  nneigh = ParEnv % NumOfNeighbours
  ALLOCATE( neigh(nneigh) )

  j = 0
  DO i=1,Parenv % Pes
    IF ( ParEnv % IsNeighbour(i) ) THEN
      j = j + 1
      neigh(j) = i-1
      IF ( j==nneigh ) EXIT
    END IF
  END DO

  CALL Send_LocIf_size( GlobalData % SplittedMatrix, nneigh, neigh )
  ALLOCATE( recv_size(nneigh) )
  CALL Recv_LocIf_size( nneigh, neigh, recv_size )

  ALLOCATE( buffer(nneigh) )
  DO i=1,nneigh
    ALLOCATE( buffer(i) % rbuf(recv_size(i)) ) 
  END DO

  ALLOCATE( requests(nneigh) )
  CALL Recv_LocIf( GlobalData % SplittedMatrix, &
      nneigh, neigh, recv_size, requests, buffer )

  !$OMP PARALLEL DO
  DO i=1,n
     v(i) = 0.0
  END DO
  !$OMP END PARALLEL DO
  DO i = 1, ParEnv % PEs
     CurrIf => GlobalData % SplittedMatrix % IfMatrix(i)

     IF ( CurrIf % NumberOfRows /= 0 ) THEN
        IfV => GlobalData % SplittedMatrix % IfVecs(i)
        IfL => GlobalData % SplittedMatrix % IfLCols(i)
        IfO => GlobalData % SplittedMatrix % IfORows(i)

        !$OMP PARALLEL PRIVATE(ColInd,j,k)
        !$OMP DO
        DO j=1,CurrIf % NumberOfRows
           IfV % IfVec(j) = 0.0
        END DO
        !$OMP END DO
        !$OMP DO
        DO j = 1, CurrIf % NumberOfRows
           IF ( Currif % RowOwner(j) /= ParEnv % MyPE ) THEN
             DO k = CurrIf % Rows(j), CurrIf % Rows(j+1) - 1
               Colind = IfL % IfVec(k)
               IF ( ColInd > 0 ) &
                 IfV % IfVec(j)=IfV % IfVec(j)+CurrIf % Values(k)*u(colind)
             END DO
           END IF
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
     END IF
  END DO

  CALL Send_LocIf( GlobalData % SplittedMatrix, nneigh, neigh )
  !----------------------------------------------------------------------
  !
  ! Compute the local part
  !
  !----------------------------------------------------------------------
  CALL CRS_MatrixVectorMultiply( InsideMatrix, u, v )

  CALL Recv_LocIf_Wait( GlobalData % SplittedMatrix, n, v,  &
       nneigh, neigh, recv_size, requests, buffer )

  DEALLOCATE( neigh, recv_size, requests )

  DO i=1,SIZE(buffer)
    DEALLOCATE( buffer(i) % rbuf )
  END DO
  DEALLOCATE( buffer )

  CALL SParIterActiveBarrier 
!----------------------------------------------------------------------
END SUBROUTINE SParMatrixVector
!----------------------------------------------------------------------


!---------------------------------------------------------------------
!> External Matrix - Vector operations (Parallel real, version)
!> Multiply vector u with the matrix which has as entries the absolute
!> values of the global matrix entries, return the result in v.
!> (HUTIter library convention)
!----------------------------------------------------------------------
  SUBROUTINE SParABSMatrixVector( u,v,ipar )
!----------------------------------------------------------------------

  IMPLICIT NONE

  ! Input parameters

  INTEGER, DIMENSION(*) :: ipar
  REAL(KIND=dp), DIMENSION(*) :: u, v

  ! Local parameters

  INTEGER :: i, j, k, l, n, ii, colind, rowind
  TYPE( IfVecT), POINTER :: IfV
  TYPE( IfLColsT), POINTER :: IfL, IfO
  TYPE (BasicMatrix_t), POINTER :: CurrIf
  TYPE (Matrix_t), POINTER :: InsideMatrix

  INTEGER :: nneigh

  TYPE(Buff_t), POINTER ::  buffer(:)

  REAL(KIND=dp) :: rsum

  REAL(KIND=dp), POINTER CONTIG :: Vals(:), Abs_Vals(:)
! REAL(KIND=dp), POINTER :: Vals(:)

  INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
! INTEGER, POINTER :: Cols(:),Rows(:)

  INTEGER, ALLOCATABLE :: neigh(:), recv_size(:), requests(:)

REAL(kind=dp) :: s, RealTime
  !*******************************************************************

  InsideMatrix => GlobalData % SplittedMatrix % InsideMatrix
  n = InsideMatrix % NumberOfRows

  !----------------------------------------------------------------------
  !
  ! Compute the interface contribution for each neighbour
  !
  !----------------------------------------------------------------------

  nneigh = ParEnv % NumOfNeighbours
  ALLOCATE( neigh(nneigh) )

  j = 0
  DO i=1,Parenv % Pes
    IF ( ParEnv % IsNeighbour(i) ) THEN
      j = j + 1
      neigh(j) = i-1
      IF ( j==nneigh ) EXIT
    END IF
  END DO

  CALL Send_LocIf_size( GlobalData % SplittedMatrix, nneigh, neigh )
  ALLOCATE( recv_size(nneigh) )
  CALL Recv_LocIf_size( nneigh, neigh, recv_size )

  ALLOCATE( buffer(nneigh) )
  DO i=1,nneigh
    ALLOCATE( buffer(i) % rbuf(recv_size(i)) ) 
  END DO

  ALLOCATE( requests(nneigh) )
  CALL Recv_LocIf( GlobalData % SplittedMatrix, &
      nneigh, neigh, recv_size, requests, buffer )

  v(1:n) = 0.0
  DO i = 1, ParEnv % PEs
     CurrIf => GlobalData % SplittedMatrix % IfMatrix(i)

     IF ( CurrIf % NumberOfRows /= 0 ) THEN
       IfV => GlobalData % SplittedMatrix % IfVecs(i)
       IfL => GlobalData % SplittedMatrix % IfLCols(i)
       IfO => GlobalData % SplittedMatrix % IfORows(i)

        IfV % IfVec(1:CurrIf % NumberOfRows) = 0.0
        DO j = 1, CurrIf % NumberOfRows
           IF ( Currif % RowOwner(j) /= ParEnv % MyPE ) THEN
             DO k = CurrIf % Rows(j), CurrIf % Rows(j+1) - 1
               Colind = IfL % IfVec(k)
               IF ( ColInd > 0 ) &
                 IfV % IfVec(j)=IfV % IfVec(j)+ABS(CurrIf % Values(k))*u(colind)
             END DO
           END IF
        END DO
     END IF
  END DO

  CALL Send_LocIf( GlobalData % SplittedMatrix, nneigh, neigh )
  !----------------------------------------------------------------------
  !
  ! Compute the local part
  !
  !----------------------------------------------------------------------


  Rows => InsideMatrix % Rows
  Cols => InsideMatrix % Cols
  Vals => InsideMatrix % Values

  IF  ( GlobalMatrix % MatvecSubr /= 0 ) THEN
    ALLOCATE(Abs_Vals(SIZE(InsideMatrix % Values)))
    Abs_Vals = ABS(Vals)
    CALL MatVecSubrExt(GlobalMatrix % MatVecSubr, &
        GlobalMatrix % SpMV, n,Rows,Cols,Abs_Vals,u,v,0)
    DEALLOCATE(Abs_Vals)
  ELSE
!$omp parallel do private(j,rsum)
    DO i = 1, n
      rsum = 0._dp
      DO j = Rows(i), Rows(i+1) - 1
        rsum = rsum + ABS(Vals(j)) * u(Cols(j))
      END DO
      v(i)=v(i)+rsum
    END DO
!$omp end parallel do
  END IF


  CALL Recv_LocIf_Wait( GlobalData % SplittedMatrix, n, v,  &
       nneigh, neigh, recv_size, requests, buffer )

  DEALLOCATE( neigh, recv_size, requests )

  DO i=1,SIZE(buffer)
    DEALLOCATE( buffer(i) % rbuf )
  END DO
  DEALLOCATE( buffer )

  CALL SParIterActiveBarrier 
!----------------------------------------------------------------------
END SUBROUTINE SParABSMatrixVector
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!> External Matrix - Vector operations (Parallel, complex version).
!> Multiply vector u with the matrix in A_val return the result in v
!> Called from HUTIter library
!----------------------------------------------------------------------
SUBROUTINE SParCMatrixVector( u, v, ipar )

  IMPLICIT NONE

  ! Input parameters

  INTEGER, DIMENSION(*) :: ipar
  COMPLEX(KIND=dp), DIMENSION(*) :: u, v


  ! Local parameters

  INTEGER :: i, j, k, l, n, colind, rowind
  TYPE (BasicMatrix_t), POINTER :: CurrIf
  TYPE (Matrix_t), POINTER :: InsideMatrix
  TYPE( IfVecT), POINTER :: IfV
  TYPE( IfLColsT), POINTER :: IfL, IfO

  COMPLEX(KIND=dp) :: A, rsum

  REAL(KIND=dp), POINTER :: Vals(:)
  INTEGER, POINTER :: Cols(:),Rows(:)
  REAL(KIND=dp), ALLOCATABLE :: buf(:)

  !*******************************************************************

  InsideMatrix => GlobalData % SplittedMatrix % InsideMatrix
  n = InsideMatrix % NumberOfRows/2

  !----------------------------------------------------------------------
  !
  ! Compute the interface contribution for each neighbour
  !
  !----------------------------------------------------------------------

  v(1:n) = 0.0
  DO i = 1, ParEnv % PEs
     CurrIf => GlobalData % SplittedMatrix % IfMatrix(i)

     IF ( CurrIf % NumberOfRows /= 0 ) THEN
        IfV => GlobalData % SplittedMatrix % IfVecs(i)
        IfL => GlobalData % SplittedMatrix % IfLCols(i)
        IfO => GlobalData % SplittedMatrix % IfORows(i)

        IfV % IfVec(1:CurrIf % NumberOfRows) = 0._dp
        DO j = 1, CurrIf % NumberOfRows / 2
          IF ( Currif % RowOwner(2*j-1) == ParEnv % MyPE ) THEN
!
! This is now dealt with in the matrix initialization, by adding the entries
! to 'InsideMatrix'...
!            rowind = (IfO % IfVec(2*j-1)+1)/2
!            DO k = CurrIf % Rows(2*j-1), CurrIf % Rows(2*j)-1,2
!              Colind = (IfL % IfVec(k)+1)/2
!              IF ( ColInd > 0 ) THEN
!                A = CMPLX( CurrIf % Values(k), -CurrIf % Values(k+1), KIND=dp )
!                v(rowind) = v(rowind) + A * u(ColInd)
!              END IF
!            END DO
          ELSE
            DO k = CurrIf % Rows(2*j-1), CurrIf % Rows(2*j)-1, 2
               IF ( IfL % IfVec(k) > 0 ) THEN
                 ColInd = (IfL % IfVec(k)+1)/2
                 A = CMPLX( CurrIf % Values(k), -CurrIf % Values(k+1), KIND=dp )
                 A = A * u(ColInd)
                 IfV % IfVec(2*j-1) = IfV % IfVec(2*j-1) + REAL(A)
                 IfV % IfVec(2*j-0) = IfV % IfVec(2*j-0) + AIMAG(A)
               END IF
            END DO
          END IF
        END DO
     END IF
  END DO

  CALL Send_LocIf_Old( GlobalData % SplittedMatrix )

  !----------------------------------------------------------------------
  !
  ! Compute the local part
  !
  !----------------------------------------------------------------------

  Rows => InsideMatrix % Rows
  Cols => InsideMatrix % Cols
  Vals => InsideMatrix % Values

!$omp parallel do private(j,A,rsum)
  DO i = 1, n
     rsum=0._dp
     DO j = Rows(2*i-1), Rows(2*i)-1, 2
        A = CMPLX( Vals(j), -Vals(j+1), KIND=dp )
        rsum = rsum + A * u(Cols(j+1)/2)
     END DO
     v(i)=v(i)+rsum
  END DO
!$omp end parallel do
  !----------------------------------------------------------------------
  !
  ! Receive interface parts of the vector and sum them to vector v
  !
  !----------------------------------------------------------------------
  ALLOCATE( buf(2*n) )
  buf = 0._dp
  CALL Recv_LocIf_Old( GlobalData % SplittedMatrix, 2*n, buf )

  DO i=1,n
    v(i) = v(i) + CMPLX( buf(2*i-1), buf(2*i), KIND=dp )
  END DO
  DEALLOCATE( buf )
!*********************************************************************
END SUBROUTINE SParCMatrixVector
!*********************************************************************




!--------------------------------------------------------------------------
!> This routine is used to count the nodes which have connections to
!> neighbours the count of connections for neighbour. Then this information
!> is used to allocate some communication buffers.
!--------------------------------------------------------------------------
SUBROUTINE CountNeighbourConns( SourceMatrix, SplittedMatrix, ParallelInfo )
!--------------------------------------------------------------------------
  USE Types
  IMPLICIT NONE

  TYPE (SplittedMatrixT) :: SplittedMatrix
  TYPE (ParallelInfo_t) :: ParallelInfo
  TYPE (Matrix_t) :: SourceMatrix

  ! Local variables

  INTEGER :: i, j, k
  INTEGER, DIMENSION(:), ALLOCATABLE :: ResEPerNB, RHSEPerNB

  !*******************************************************************

  IF ( .NOT. ASSOCIATED( SplittedMatrix % ResBuf ) ) THEN
     ALLOCATE( SplittedMatrix % ResBuf(ParEnv % PEs) )
  END IF
  IF ( .NOT. ASSOCIATED( SplittedMatrix % RHS ) ) THEN
     ALLOCATE( SplittedMatrix % RHS(ParEnv % PEs) )
  END IF

  !----------------------------------------------------------------------
  !
  ! Count the nodes per neighbour which are shared between neighbours
  !
  !----------------------------------------------------------------------

  ALLOCATE( ResEPerNB( ParEnv % PEs ) )
  ALLOCATE( RHSEPerNB( ParEnv % PEs ) )
  ResEPerNB = 0; RHSEPerNB = 0

  DO i = 1, SourceMatrix % NumberOfRows
     IF ( ParallelInfo % GInterface(i) ) THEN
        IF ( ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
           DO j = 1, SIZE( ParallelInfo % NeighbourList(i) % Neighbours )
               IF ( ParallelInfo % NeighbourList(i) % Neighbours(j)/=ParEnv % MyPE ) THEN
                ResEPerNB(ParallelInfo % NeighbourList(i) % Neighbours(j)+1) = &
                ResEPerNB(ParallelInfo % NeighbourList(i) % Neighbours(j)+1) + 1
             END IF
            END DO
        ELSE
          RHSEPerNB(ParallelInfo % NeighbourList(i) % Neighbours(1)+1) = &
             RHSEPerNB(ParallelInfo % NeighbourList(i) % Neighbours(1)+1) + 1
        END IF
     END IF
  END DO

  !----------------------------------------------------------------------
  !
  ! Allocate some buffers for communication
  !
  !----------------------------------------------------------------------

  RHSEperNB = RHSEperNB
  ResEPerNB = ResEPerNB

  DO i = 1, ParEnv % PEs
     NULLIFY( SplittedMatrix % RHS(i) % RHSVec )
     NULLIFY( SplittedMatrix % RHS(i) % RHSInd )
     IF ( RHSEPerNB(i) /= 0 ) THEN
        ALLOCATE( SplittedMatrix % RHS(i) % RHSVec( RHSEperNB(i) ) )
        ALLOCATE( SplittedMatrix % RHS(i) % RHSind( RHSEperNB(i) ) )
     END IF

!    NULLIFY( SplittedMatrix % ResBuf(i) % ResVal )
!    NULLIFY( SplittedMatrix % ResBuf(i) % ResInd )
     IF ( ResEPerNB(i) /= 0 ) THEN
        ALLOCATE( SplittedMatrix % ResBuf(i) % ResVal( ResEPerNB(i) ) )
        ALLOCATE( SplittedMatrix % ResBuf(i) % ResInd( ResEPerNB(i) ) )
     END IF

  END DO

  DEALLOCATE( ResEPerNB, RHSEPerNB )
!*********************************************************************
END SUBROUTINE CountNeighbourConns
!*********************************************************************


!-----------------------------------------------------------------------
!> This subroutine combines indices of two matrices in CRS format.
!> Row and column indices are assumed to be in sorted order.
!-----------------------------------------------------------------------
SUBROUTINE CombineCRSMatIndices ( SMat1, SMat2, DMat )
!-----------------------------------------------------------------------
  USE Types
  IMPLICIT NONE

  TYPE (BasicMatrix_t), TARGET :: SMat1, SMat2, DMat
! External routines

! EXTERNAL SearchIAItem
! INTEGER :: SearchIAItem

  ! Local variables

  INTEGER :: i, j, k, i1, i2, j1, j2, ind, ind1, ind2, DRows, DCols, row, col

  INTEGER, POINTER :: cols(:)
  LOGICAL, POINTER :: done(:)

  !*******************************************************************

  IF ( SMat1 % NumberOfRows == 0 .AND. SMat2 % NumberOfRows == 0 ) THEN

     RETURN

  ELSE IF ( SMat1 % NumberOfRows == 0 ) THEN

     ALLOCATE( DMat % Rows(  SMat2 % NumberOfRows + 1) )
     ALLOCATE( DMat % GRows( SMat2 % NumberOfRows ) )
     ALLOCATE( DMat % RowOwner( SMat2 % NumberOfRows ) )
     ALLOCATE( DMat % Cols( SMat2 % Rows(SMat2 % NumberOfRows + 1)-1 ) )

     DMat % NumberOfRows = SMat2 % NumberOfRows
     DMat % Rows = SMat2 % Rows(1:SMat2 % NumberOfRows+1)
     DMat % GRows = SMat2 % GRows(1:SMat2 % NumberOfRows)
     DMat % Cols = SMat2 % Cols(1:SIZE(DMat % Cols))
     DMat % RowOwner = SMat2 % RowOwner(1:SMat2 % NumberOfRows)

     RETURN

  ELSE IF ( SMat2 % NumberOfRows == 0 ) THEN

     ALLOCATE( DMat % Rows( SMat1 % NumberOfRows + 1) )
     ALLOCATE( DMat % GRows( SMat1 % NumberOfRows ) )
     ALLOCATE( DMat % RowOwner( SMat1 % NumberOfRows ) )
     ALLOCATE( DMat % Cols( SMat1 % Rows(SMat1 % NumberOfRows + 1)-1 ) )

     DMat % NumberOfRows = SMat1 % NumberOfRows
     DMat % Rows = SMat1 % Rows(1:SMat1 % NumberOfRows+1)
     DMat % GRows = SMat1 % GRows(1:SMat1 % NumberOfRows)
     DMat % Cols = SMat1 % Cols(1:SIZE(DMat % Cols))
     DMat % RowOwner = SMat1 % RowOwner(1:SMat1 % NumberOfRows)
     RETURN

  END IF
        
  !----------------------------------------------------------------------
  !
  ! First we have to compute the strorage allocations
  !
  !----------------------------------------------------------------------
  CALL CRS_SortBasicMatrix( SMat1 )
  CALL CRS_SortBasicMatrix( SMat2 )

  DRows = SMat2 % NumberOfRows;
  DO i = 1, SMat1 % NumberOfRows
     ind1 = SearchIAItem( SMat2 % NumberOfRows, &
          SMat2 % GRows, SMat1 % GRows(i) )
     IF ( Ind1 == -1 ) DRows = DRows + 1
  END DO

  DCols = SIZE(SMat1 % Cols) + SIZE(SMat2 % Cols)
  DMat % NumberOfRows = DRows
  ALLOCATE( DMat % Rows(DRows + 1) )
  ALLOCATE( DMat % Cols(DCols) )
  ALLOCATE( DMat % GRows(DRows) )
  ALLOCATE( DMat % RowOwner( DRows ) )

  !----------------------------------------------------------------------
  !
  ! Then we combine the index structures of the two CRS matrices
  !
  !----------------------------------------------------------------------

  row = 1; col = 1; i1 = 1; i2 = 1
  ALLOCATE( Done( Smat2 % NumberOfRows ) )
  done = .FALSE.

  DO WHILE ( i1 <= SMat1 % NumberOfRows .OR. &
             i2 <= SMat2 % NumberOfRows )

     ind = -1
     IF ( i1 <= SMat1 % NumberOfRows ) THEN
        ind = SearchIAItem( SMat2 % NumberOfRows, &
          SMat2 % GRows, SMat1 % GRows(i1) )
     END IF

     IF ( i1 > SMat1 % NumberOfRows ) THEN

        DO i=1,Smat2 % NumberOfRows
           IF ( .NOT.done(i) ) THEN
              DMat % Rows(row) = Col
              DMat % GRows(row)    = SMat2 % GRows(i)
              DMat % RowOwner(row) = SMat2 % RowOwner(i)
              Row = Row + 1
              DO k = SMat2 % Rows(i), SMat2 % Rows(i+1) - 1
                 DMat % Cols(col) = SMat2 % Cols(k)
                 Col = Col + 1
              END DO
              i2 = i2 + 1
           END IF
        END DO

     ELSE IF ( i2 > SMat2 % NumberOfRows .OR. Ind == -1 ) THEN

        DMat % Rows(Row)     = Col
        DMat % GRows(Row)    = SMat1 % GRows(i1)
        DMat % RowOwner(Row) = SMat1 % RowOwner(i1)
        Row = Row + 1
        DO k = SMat1 % Rows(i1), SMat1 % Rows(i1+1) - 1
           DMat % Cols(Col) = SMat1 % Cols(k)
           Col = Col + 1
        END DO
        i1 = i1 + 1

     ELSE IF ( Ind /= -1 ) THEN
              
        DMat % Rows(Row)  = Col
        DMat % GRows(Row) = SMat1 % GRows(i1)
        DMat % RowOwner(Row) = SMat1 % RowOwner(i1)
        Row = Row + 1
        j1 = SMat1 % Rows(i1)
        j2 = SMat2 % Rows(Ind)

        DO WHILE ( j1 < SMat1 % Rows(i1+1) .OR.  j2 < SMat2 % Rows(Ind+1) )
          IF ( j1 <  SMat1 % Rows(i1+1)) THEN
            IF (j2 >= SMat2 % Rows(Ind+1) ) THEN

              DMat % Cols(col) = SMat1 % Cols(j1)
              Col = Col + 1
              j1 = j1 + 1
              CYCLE
            ELSE IF ( SMat1 % Cols(j1) < SMat2 % Cols(j2) ) THEN

              DMat % Cols(col) = SMat1 % Cols(j1)
              Col = Col + 1
              j1 = j1 + 1
              CYCLE
            ELSE IF (j2 >= SMat2 % Rows(Ind+1) ) THEN
              
              DMat % Cols(col) = SMat1 % Cols(j1)
              Col = Col + 1
              j1 = j1 + 1
              CYCLE
            ELSE IF ( SMat1 % Cols(j1) < SMat2 % Cols(j2) ) THEN

              DMat % Cols(col) = SMat1 % Cols(j1)
              Col = Col + 1
              j1 = j1 + 1
              CYCLE
            END IF
          END IF
          IF ( j2 <  SMat2 % Rows(Ind+1) ) THEN
            IF ( j1 >= SMat1 % Rows(i1+1) ) THEN

              DMat % Cols(col) = SMat2 % Cols(j2)
              Col = Col + 1
              j2 = j2 + 1
              CYCLE
            ELSE IF ( SMat2 % Cols(j2) < SMat1 % Cols(j1) ) THEN

              DMat % Cols(col) = SMat2 % Cols(j2)
              Col = Col + 1
              j2 = j2 + 1
              CYCLE
            END IF
          END IF
          IF ( SMat1 % Cols(j1) == SMat2 % Cols(j2) ) THEN
            
            DMat % Cols(col) = SMat1 % Cols(j1)
            Col = Col + 1
            j1 = j1 + 1
            j2 = j2 + 1
            CYCLE
          END IF

          ! Should not happen
          CALL Fatal('CombineCRSMatIndices','Internal error while merging matrix rows')
        END DO

        Done(Ind) = .TRUE.
        i1 = i1 + 1
        i2 = i2 + 1
     END IF
  END DO
  DMat % Rows(Row) = Col

  ALLOCATE(Cols(col-1))
  Cols(1:col-1) = DMat % Cols(1:col-1)
  DEALLOCATE(DMAT % Cols)
  ALLOCATE(DMAT % Cols(col-1))
  DMAT % Cols = Cols
  DEALLOCATE( Done, Cols )
!*********************************************************************
END SUBROUTINE CombineCRSMatIndices
!*********************************************************************



!----------------------------------------------------------------------
!> Update all necessary global structures after the local matrix has
!> been built.
!----------------------------------------------------------------------
SUBROUTINE GlueFinalize( SourceMatrix, SplittedMatrix, ParallelInfo )
  USE Types
  IMPLICIT NONE

  TYPE (ParallelInfo_t) :: ParallelInfo
  TYPE(Matrix_t) :: SourceMatrix
  TYPE (SplittedMatrixT) :: SplittedMatrix

  ! External routines

! EXTERNAL SearchIAItem
! INTEGER :: SearchIAItem

  ! Local variables

  INTEGER, ALLOCATABLE :: RevDOFList(:)

  TYPE (BasicMatrix_t), POINTER :: CurrIf
  TYPE (Matrix_t), POINTER :: InsideMatrix
  INTEGER :: i, j, k, l, RowInd, Rows, ColInd, ColIndA
  TYPE (BasicMatrix_t), DIMENSION(:), ALLOCATABLE :: RecvdIfMatrix

  LOGICAL :: Found, NeedMass, NeedDamp, NeedPrec, NeedILU

  !*******************************************************************

  !----------------------------------------------------------------------
  !
  ! Exchange interface block values stored in GlueLocalMatrix routine
  !
  !----------------------------------------------------------------------
  InsideMatrix => SplittedMatrix % InsideMatrix

  NeedILU = .FALSE.
  DO i=1,Parenv % PEs
    IF ( SplittedMatrix % IfMatrix(i) % NumberOfRows>0 ) THEN
      NeedILU = ALLOCATED(SplittedMatrix % IfMatrix(i) % ILUValues)
      EXIT
    END IF
  END DO

  NeedMass = ASSOCIATED(InsideMatrix % MassValues)
  NeedDamp = ASSOCIATED(InsideMatrix % DampValues)
  NeedPrec = ASSOCIATED(InsideMatrix % PrecValues)

  ALLOCATE( RecvdIfMatrix(ParEnv % PEs) )
  RecvdIfMatrix(:) % NumberOfRows = 0

  CALL ExchangeIfValues(SplittedMatrix % NbsIfMatrix, &
     RecvdIfMatrix, NeedMass, NeedDamp, NeedPrec, NeedILU )

  !----------------------------------------------------------------------
  !
  ! Add received interface block values into our own interface blocks or
  ! into InsideMatrix.
  !
  !----------------------------------------------------------------------
  ALLOCATE(RevDofList(SourceMatrix % NumberOfRows))
  RevDOFList = 0
  j = 0
  DO i = 1, SourceMatrix % NumberOfRows
     IF (ParallelInfo % NeighbourList(i) % Neighbours(1)==ParEnv % MyPE) THEN
       j = j + 1
       RevDofList(i) = j
     END IF
  END DO

  DO i = 1, ParEnv % PEs
     IF ( RecvdIfMatrix(i) % NumberOfRows == 0 ) CYCLE

     CurrIf => SplittedMatrix % IfMatrix(i)
     DO j = 1, RecvdIfMatrix(i) % NumberOfRows

        RowInd = SearchIAItem( CurrIf % NumberOfRows, &
             CurrIf % GRows, RecvdIfMatrix(i) % GRows(j) )

        Found=.FALSE.
        IF ( RowInd>0 ) THEN
           ! Value to be added to IfMatrix

           IF (CurrIf % RowOwner(RowInd)/=ParEnv % MyPE) THEN
             DO  k=RecvdIfMatrix(i) % Rows(j),RecvdIfMatrix(i) % Rows(j+1)-1
               DO l = CurrIf % Rows(RowInd), CurrIf % Rows(RowInd+1)-1
                  IF ( RecvdIfMatrix(i) % Cols(k) == CurrIf % Cols(l) ) THEN
                     CurrIf % Values(l) = CurrIf % Values(l) + &
                           RecvdIfMatrix(i) % Values(k)
                     IF ( NeedPrec ) &
                        CurrIf % PrecValues(l) = CurrIf % PrecValues(l) + &
                               RecvdIfMatrix(i) % PrecValues(k)
                     IF ( NeedMass ) &
                        CurrIf % MassValues(l) = CurrIf % MassValues(l) + &
                               RecvdIfMatrix(i) % MassValues(k)
                     IF ( NeedDamp ) &
                        CurrIf % DampValues(l) = CurrIf % DampValues(l) + &
                                 RecvdIfMatrix(i) % DampValues(k)
                     IF ( NeedILU) &
                        CurrIf % ILUValues(l) = CurrIf % ILUValues(l) + &
                               RecvdIfMatrix(i) % ILUValues(k)
                      EXIT
                   END IF
               END DO
             END DO
             Found=.TRUE.
           END IF
        END IF

        IF (.NOT.Found) THEN
           ! Value to be added to InsideMatrix
           RowInd = SearchIAItem( InsideMatrix % NumberOfRows, &
               InsideMatrix % GRows, RecvdIfMatrix(i) % GRows(j), &
                   InsideMatrix % Gorder )

           IF ( RowInd > 0 ) THEN
              DO k = RecvdIfMatrix(i) % Rows(j), RecvdIfMatrix(i) % Rows(j+1) - 1
                 l = RecvdIfMatrix(i) % Cols(k)

!                ColIndA = SearchNode(ParallelInfo,l,Order=SourceMatrix % Perm)
!XYXY
                 ColIndA =  SearchNode(ParallelInfo,l, Order=ParallelInfo % Gorder )
                 IF (ColIndA>0) ColIndA = RevDOFList(ColIndA)

                 IF ( ColIndA  <= 0 ) CYCLE

                 Found = .FALSE.
                 DO l = InsideMatrix % Rows(RowInd), InsideMatrix % Rows(RowInd+1) - 1
                     IF ( ColIndA == InsideMatrix % Cols(l) ) THEN
                        InsideMatrix % Values(l) = InsideMatrix % Values(l) + &
                                RecvdIfMatrix(i) % Values(k)
                       IF ( NeedPrec ) &
                          InsideMatrix % PrecValues(l) = InsideMatrix % PrecValues(l) + &
                                   RecvdIfMatrix(i) % PrecValues(k)
                       IF ( NeedMass ) &
                          InsideMatrix % MassValues(l) = InsideMatrix % MassValues(l) + &
                                   RecvdIfMatrix(i) % MassValues(k)
                       IF ( NeedDamp ) &
                          InsideMatrix % DampValues(l) = InsideMatrix % DampValues(l) + &
                                   RecvdIfMatrix(i) % DampValues(k)
                       IF ( NeedILU ) &
                          InsideMatrix % ILUValues(l) = InsideMatrix % ILUValues(l) + &
                                   RecvdIfMatrix(i) % ILUValues(k)
                       Found = .TRUE.
                       EXIT
                     END IF
                 END DO

                 IF ( .NOT. Found ) THEN
                   IF(RecvdIfMatrix(i) % Values(k) /= 0._dp ) THEN
                     PRINT*,ParEnv % MyPE,  'GlueFinalize: This should not happen 0'
                     PRINT*,ParEnv % MyPE, i-1, RecvdIfMatrix(i) % GRows(j), RowInd, &
                                 RecvdIfMatrix(i) % Cols(k), ColIndA, recvdifmatrix(i) % values(k)
                   END IF
                 END IF
              END DO
           ELSE
              PRINT*,ParEnv % MyPE, 'GlueFinalize: This should not happen 1', RowInd
              CALL FLUSH(6)
           END IF
        END IF
     END DO
  END DO

  DEALLOCATE( RevDOFList )
  DO i=1,ParEnv % PEs
     IF ( RecvdIfMatrix(i) % NumberOfRows > 0 ) THEN
       DEALLOCATE( RecvdIfMatrix(i) % Rows,  RecvdIfMatrix(i) % Cols,  &
            RecvdIfMatrix(i) % GRows, RecvdIfMatrix(i) % Values )

       IF (ALLOCATED(RecvdIfMatrix(i) % ILUValues)) &
         DEALLOCATE(RecvdIfMatrix(i) % ILUValues)

       IF (ALLOCATED(RecvdIfMatrix(i) % MassValues)) &
         DEALLOCATE(RecvdIfMatrix(i) % MassValues)

       IF (ALLOCATED(RecvdIfMatrix(i) % DampValues)) &
         DEALLOCATE(RecvdIfMatrix(i) % DampValues)

       IF (ALLOCATED(RecvdIfMatrix(i) % PrecValues)) &
         DEALLOCATE(RecvdIfMatrix(i) % PrecValues)
     END IF
  END DO
  DEALLOCATE( RecvdIfMatrix )
END SUBROUTINE GlueFinalize
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!> Compress the given interface matrix deleting the inside connections.
!> Note: modifies the structure of the input matrix RecvdIfMatrix.
!----------------------------------------------------------------------
SUBROUTINE ClearInsideC( SourceMatrix, InsideMatrix, &
            RecvdIfMatrix, ParallelInfo )

  USE Types
  IMPLICIT NONE

  ! Parameters

  TYPE (Matrix_t) :: SourceMatrix, InsideMatrix
  TYPE (ParallelInfo_t) :: ParallelInfo
  TYPE (BasicMatrix_t), DIMENSION(:) :: RecvdIfMatrix

  ! External routines

! EXTERNAL SearchIAItem
! INTEGER :: SearchIAItem

  ! Local variables

  INTEGER :: NewRow, NewCol,old_nv,nc
  INTEGER :: p,i,j,k,l,RowInd,ColInd,GCol
  
  !*********************************************************************
  !
  ! Compression of the matrix is done in place and lengths are
  ! stored in New* variables
  !
  DO p = 1, ParEnv % PEs
     IF ( RecvdIfMatrix(p) % NumberOfRows <= 0 ) CYCLE

     NewRow = 1; NewCol = 1; old_nv = 1
     DO i = 1, RecvdIfMatrix(p) % NumberOfRows

        NC = 0
        RowInd = SearchIAItem( InsideMatrix % NumberOfRows,  &
          InsideMatrix % GRows, RecvdIfMatrix(p) % GRows(i), &
              InsideMatrix % Gorder )

        IF ( RowInd /= -1 ) THEN
           DO j = RecvdIfMatrix(p) % Rows(i),RecvdIfMatrix(p) % Rows(i+1) - 1

!             GCol = SearchNode( ParallelInfo,  RecvdIfMatrix(p) % Cols(j),&
!                       Order = SourceMatrix % Perm )
!XYXY
              GCol = SearchNode( ParallelInfo,  RecvdIfMatrix(p) % Cols(j), Order=ParallelInfo % Gorder )

!             GCol = SourceMatrix % Perm( Gcol )

              ColInd = -1
              IF ( Gcol  > 0 ) THEN
                DO k = InsideMatrix % Rows(RowInd), &
                       InsideMatrix % Rows(RowInd+1) - 1
                    IF ( GCol == InsideMatrix % Cols(k) ) THEN
                       ColInd = GCol
                       EXIT
                    END IF
                END DO
              END IF

              IF ( ColInd == -1 ) THEN
                 RecvdIfMatrix(p) % Cols(NewCol) = RecvdIfMatrix(p) % Cols(j)
                 NewCol = NewCol + 1
                 NC = 1
              END IF
           END DO
        ELSE
           DO j = RecvdIfMatrix(p) % Rows(i), RecvdIfMatrix(p) % Rows(i+1)-1
              RecvdIfMatrix(p) % Cols(NewCol) = RecvdIfMatrix(p) % Cols(j)
              NewCol = NewCol + 1
              NC = 1
           END DO

        END IF

        IF ( NC /= 0 ) THEN
           RecvdIfMatrix(p) % GRows(NewRow)    = RecvdIfMatrix(p) % GRows(i)
           RecvdIfMatrix(p) % RowOwner(NewRow) = RecvdIfMatrix(p) % RowOwner(i)
           RecvdIfMatrix(p) % Rows(NewRow)  = old_nv
           NewRow = NewRow + 1
        END IF
        old_nv = NewCol
           
     END DO
     RecvdIfMatrix(p) % Rows(NewRow) = NewCol
     RecvdIfMatrix(p) % NumberOfRows  = NewRow - 1
  END DO
END SUBROUTINE ClearInsideC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Convert the original local DoF numbering to a compressed one to
!> enable direct use of column indices in matrix operations.
!-----------------------------------------------------------------------
SUBROUTINE RenumberDOFs( SourceMatrix, SplittedMatrix, ParallelInfo )

  USE Types

  IMPLICIT NONE

  TYPE (SplittedMatrixT) :: SplittedMatrix
  TYPE (ParallelInfo_t) :: ParallelInfo
  TYPE( Matrix_t) :: SourceMatrix

  ! Local variables

  INTEGER, DIMENSION(:), ALLOCATABLE :: RevDofList
  INTEGER :: i, j, k, Inside
  TYPE (Matrix_t), POINTER :: InsideMatrix

  !*******************************************************************

  ! Construct a list to convert original DOF (Degrees of Freedom)
  ! numbering to truncated one (for InsideMatrix).

  ALLOCATE( RevDofList( SourceMatrix % NumberOfRows ) )

  Inside = 0
  DO i = 1, SourceMatrix % NumberOfRows
     IF ( ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
        Inside = Inside + 1
        RevDofList(i) = Inside
     ELSE
        RevDofList(i) = -1
     END IF
  END DO

  ! Scan the InsideMatrix and change the numbering

  InsideMatrix => SplittedMatrix % InsideMatrix
  DO i = 1, InsideMatrix % NumberOfRows
     DO j = InsideMatrix % Rows(i), InsideMatrix % Rows(i+1)-1
        InsideMatrix % Cols(j) = RevDofList( InsideMatrix % Cols(j) )
     END DO
  END DO

  DEALLOCATE( RevDofList )
!********************************************************************
END SUBROUTINE RenumberDOFs
!********************************************************************

END MODULE SParIterSolve

!> \}
