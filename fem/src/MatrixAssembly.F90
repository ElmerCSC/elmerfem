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
! *  Utilities for *Solver - routines
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 28 Sep 1998
! *
! *****************************************************************************/

!> Basic utilities for creating and manipulating the global matrix.
!> These routines should only have dependence on the CRS, List, and Band matrix
!> routines. 
!------------------------------------------------------------------------------

!> \ingroup ElmerLib
!> \{


MODULE MatrixAssembly

#include "../config.h"


   USE ListMatrix
   USE CRSMatrix
   USE BandMatrix
   USE PElementMaps, ONLY : isActivePElement

   
   IMPLICIT NONE

   INTERFACE CondensateP
     MODULE PROCEDURE CondensatePR, CondensatePC
   END INTERFACE CondensateP

CONTAINS
   
!> Sets the matrix element to a desired value. 
!------------------------------------------------------------------------------
   SUBROUTINE SetMatrixElement( A, i, j, VALUE )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A  !< Structure holding the matrix
     INTEGER :: i                            !< Row index
     INTEGER :: j                            !< Column index
     REAL(KIND=dp) :: VALUE                  !< Value to be obtained
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_SetMatrixElement( A, i, j, VALUE )
         IF(A % FORMAT == MATRIX_LIST) THEN
           CALL List_toListMatrix(A)
           CALL List_SetMatrixElement( A % ListMatrix, i, j, VALUE )
         END IF

       CASE( MATRIX_LIST )
         CALL List_SetMatrixElement( A % ListMatrix, i, j, VALUE )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_SetMatrixElement( A, i, j, VALUE )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE SetMatrixElement
!------------------------------------------------------------------------------

!> Gets a matrix element. 
!------------------------------------------------------------------------------
   FUNCTION GetMatrixElement( A, i, j ) RESULT ( VALUE )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A  !< Structure holding the matrix
     INTEGER :: i                            !< Row index
     INTEGER :: j                            !< Column index
     REAL(KIND=dp) :: VALUE                  !< Value to be obtained
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         VALUE = CRS_GetMatrixElement( A, i, j )

      CASE( MATRIX_LIST )
         VALUE = List_GetMatrixElement( A % ListMatrix, i, j )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         VALUE = Band_GetMatrixElement( A, i, j )
     END SELECT
!------------------------------------------------------------------------------
   END FUNCTION GetMatrixElement
!------------------------------------------------------------------------------

!> Changes the value of a given matrix element.
!------------------------------------------------------------------------------
   FUNCTION ChangeMatrixElement( A, i, j, NewValue ) RESULT ( OldValue )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i,j
     REAL(KIND=dp) :: NewValue, OldValue
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         OldValue = CRS_ChangeMatrixElement( A, i, j, NewValue )

       CASE DEFAULT
         CALL Warn('ChangeMatrixElement','Not implemented for this type')

     END SELECT
!------------------------------------------------------------------------------
   END FUNCTION ChangeMatrixElement
!------------------------------------------------------------------------------


!> Adds to the value of a given matrix element.
!------------------------------------------------------------------------------
   SUBROUTINE AddToMatrixElement( A, i, j,VALUE )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i,j
     REAL(KIND=dp) :: VALUE
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_AddToMatrixElement( A, i, j, VALUE )
         IF(A % FORMAT == MATRIX_LIST) THEN
           CALL List_toListMatrix(A)
           CALL List_AddToMatrixElement( A % ListMatrix, i, j, VALUE )
         END IF

      CASE( MATRIX_LIST )
         CALL List_AddToMatrixElement( A % ListMatrix, i, j, VALUE )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_AddToMatrixElement( A, i, j, VALUE )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE AddToMatrixElement
!------------------------------------------------------------------------------

!> Adds CMPLX value to the value of a given CMPLX matrix element. -ettaka
!------------------------------------------------------------------------------
  SUBROUTINE AddToCmplxMatrixElement(CM, RowId, ColId, Re, Im)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: CM
    INTEGER :: RowId, ColId
    REAL(KIND=dp) :: Re, Im

    CALL AddToMatrixElement(CM, RowId, ColId, Re)
    CALL AddToMatrixElement(CM, RowId, ColId+1, -Im)
    CALL AddToMatrixElement(CM, RowId+1, ColId, Im)
    CALL AddToMatrixElement(CM, RowId+1, ColId+1, Re)

!------------------------------------------------------------------------------
  END SUBROUTINE AddToCmplxMatrixElement
!------------------------------------------------------------------------------

!> Moves a matrix element from one position adding it to the value of another one.
!------------------------------------------------------------------------------
   SUBROUTINE MoveMatrixElement( A, i1, j1, i2, j2 )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i1,j1,i2,j2
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: VALUE

     VALUE = ChangeMatrixElement(A, i1, j1, 0.0_dp)
     CALL AddToMatrixElement(A, i2, j2, VALUE )
     
!------------------------------------------------------------------------------
   END SUBROUTINE MoveMatrixElement
!------------------------------------------------------------------------------


!> Zeros a row in matrix.
!------------------------------------------------------------------------------
   SUBROUTINE ZeroRow( A, n )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A  !< Structure holding the matrix 
      INTEGER :: n                           !< Row to be zerored.
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_ZeroRow( A,n )

       CASE( MATRIX_LIST )
         CALL List_ZeroRow( A % ListMatrix, n )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_ZeroRow( A,n )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE ZeroRow
!------------------------------------------------------------------------------

!> Moves a row and and sums it with the values of a second one, optionally 
!> multiplying with a constant.
!------------------------------------------------------------------------------
   SUBROUTINE MoveRow( A, n1, n2, Coeff, StayCoeff, MoveCoeff )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n1, n2
     REAL(KIND=dp), OPTIONAL :: Coeff, StayCoeff, MoveCoeff
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_MoveRow( A,n1,n2,Coeff,StayCoeff )

         ! If entries are not found the format is changed on-the-fly
         IF( A % FORMAT == MATRIX_LIST ) THEN
           CALL CRS_MoveRow( A,n1,n2,Coeff,StayCoeff ) ! does this make sense?
         END IF
         
       CASE( MATRIX_LIST )
         CALL List_MoveRow( A % ListMatrix,n1,n2,Coeff,StayCoeff )

       CASE DEFAULT
         CALL Warn('MoveRow','Not implemented for this type')
         
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE MoveRow
!------------------------------------------------------------------------------


!---------------------------------------------------------------------------
!> If we have antiperiodic DOFs in periodic system and want to do elimination
!> for conforming mesh, then we need to flip entries in stiffness/mass matrix.
!---------------------------------------------------------------------------
   SUBROUTINE FlipPeriodicLocalMatrix( Solver, n, Indexes, dofs, A )
     TYPE(Solver_t), POINTER :: Solver
     INTEGER :: n, dofs
     INTEGER :: Indexes(:)
     REAL(KIND=dp) :: A(:,:)

     LOGICAL, POINTER :: PerFlip(:)
     INTEGER :: i,j,k,l

     IF( .NOT. Solver % PeriodicFlipActive ) RETURN

     PerFlip => Solver % Mesh % PeriodicFlip           

     IF( .NOT. ANY( PerFlip( Indexes(1:n) ) ) ) RETURN
     
     IF( dofs == 1 ) THEN
       DO i=1,n
         DO j=1,n
           IF( XOR(PerFlip(Indexes(i)),PerFlip(Indexes(j))) ) THEN
             A(i,j) = -A(i,j)
           END IF
         END DO
       END DO
     ELSE
       DO i=1,n
         DO j=1,n
           IF( XOR(PerFlip(Indexes(i)),PerFlip(Indexes(j))) ) THEN
             DO k=1,dofs
               DO l=1,dofs
                 A(dofs*(i-1)+k,dofs*(j-1)+l) = -A(dofs*(i-1)+k,dofs*(j-1)+l)
               END DO
             END DO
           END IF
         END DO
       END DO       
     END IF
              
   END SUBROUTINE FlipPeriodicLocalMatrix


!---------------------------------------------------------------------------
!> If we have antiperiodic DOFs in periodic system and want to do elimination
!> for conforming mesh, then we need to flip entries in local force.
!---------------------------------------------------------------------------
   SUBROUTINE FlipPeriodicLocalForce( Solver, n, Indexes, dofs, F )
     TYPE(Solver_t), POINTER :: Solver
     INTEGER :: n, dofs
     INTEGER :: Indexes(:)
     REAL(KIND=dp) :: F(:)
     
     LOGICAL, POINTER :: PerFlip(:)
     INTEGER :: i,j

     IF( .NOT. Solver % PeriodicFlipActive ) RETURN

     PerFlip => Solver % Mesh % PeriodicFlip           
     
     IF( .NOT. ANY( PerFlip( Indexes(1:n) ) ) ) RETURN
     
     IF( dofs == 1 ) THEN
       DO i=1,n
         IF( PerFlip(Indexes(i))) F(i) = -F(i)
       END DO
     ELSE
       DO i=1,n
         IF( PerFlip(Indexes(i))) THEN
           DO j=1,dofs
             F(dofs*(i-1)+j) = -F(dofs*(i-1)+j)
           END DO
         END IF
       END DO
     END IF
          
   END SUBROUTINE FlipPeriodicLocalForce


!---------------------------------------------------------------------------
!> Check if there is something to flip.
!---------------------------------------------------------------------------
   FUNCTION AnyFlipPeriodic( Solver, n, Indexes ) RESULT ( DoFlip ) 
     TYPE(Solver_t), POINTER :: Solver
     INTEGER :: n
     INTEGER :: Indexes(:)
     LOGICAL :: DoFlip 
     
     LOGICAL, POINTER :: PerFlip(:)

     DoFlip = .FALSE.
     IF( .NOT. Solver % PeriodicFlipActive ) RETURN
    
     PerFlip => Solver % Mesh % PeriodicFlip                
     DoFlip = ANY( PerFlip(Indexes(1:n)))

   END FUNCTION AnyFlipPeriodic

   
   
!> Glues a local matrix to the global one.
!------------------------------------------------------------------------------
   SUBROUTINE GlueLocalSubMatrix( A,row0,col0,Nrow,Ncol,RowInds,ColInds,&
       RowDofs,ColDofs,LocalMatrix )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: LocalMatrix(:,:)
     TYPE(Matrix_t) :: A
     INTEGER :: Nrow,Ncol,RowDofs,ColDofs,Col0,Row0,RowInds(:),ColInds(:)
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )

       CASE( MATRIX_CRS )       
         CALL CRS_GlueLocalSubMatrix( A,row0,col0,Nrow,Ncol,RowInds,ColInds,&
             RowDofs,ColDofs,LocalMatrix )
      
       CASE( MATRIX_LIST )
         CALL List_GlueLocalSubMatrix( A % ListMatrix,row0,col0,Nrow,Ncol,RowInds,ColInds,&
             RowDofs,ColDofs,LocalMatrix )
        
       CASE DEFAULT
         CALL Warn('GlueLocalSubMatrix','Not implemented for this type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE GlueLocalSubMatrix
!------------------------------------------------------------------------------


!> Create a copy of the matrix entries A % Values to A % BulkValues. 
!> Optionally the entries of the RHS vector, the mass matrix and the damping matrix
!> may also be copied. The RHS vector is copied by default, while the mass and
!> damping matrices are copied only if asked. 
!------------------------------------------------------------------------------
   SUBROUTINE CopyBulkMatrix( A, BulkMass, BulkDamp, BulkRHS )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     LOGICAL, OPTIONAL :: BulkMass, BulkDamp, BulkRHS

     INTEGER :: i,n
     LOGICAL :: CopyRHS

     IF (PRESENT(BulkRHS)) THEN
       CopyRHS = BulkRHS
     ELSE
       CopyRHS = .TRUE.
     END IF
     
     IF (CopyRHS) THEN
       n = SIZE( A % Rhs )
       IF( ASSOCIATED( A % BulkRhs ) ) THEN
         IF( SIZE( A % BulkRhs ) /= n ) THEN
           DEALLOCATE( A % BulkRhs ) 
           A % BulkRHS => NULL()
         END IF
       END IF
       IF ( .NOT. ASSOCIATED( A % BulkRHS ) ) THEN
         ALLOCATE( A % BulkRHS( n ) )
       END IF
       DO i=1,n
         A % BulkRHS(i) = A % Rhs(i)
       END DO
     END IF
     
     n = SIZE( A % Values )
     IF( ASSOCIATED( A % BulkValues ) ) THEN
       IF( SIZE( A % BulkValues ) /= n ) THEN
          DEALLOCATE( A % BulkValues ) 
          A % BulkValues => NULL()
       END IF
     END IF
     IF ( .NOT. ASSOCIATED( A % BulkValues ) ) THEN
       ALLOCATE( A % BulkValues( n ) )
     END IF

     DO i=1,n
       A % BulkValues(i) = A % Values(i)
     END DO

     IF( PRESENT( BulkMass ) .AND. ASSOCIATED( A % MassValues) ) THEN
       IF( BulkMass ) THEN
         n = SIZE( A % MassValues )
         IF( ASSOCIATED( A % BulkMassValues ) ) THEN
           IF( SIZE( A % BulkMassValues ) /= n ) THEN
             DEALLOCATE( A % BulkMassValues ) 
             A % BulkMassValues => NULL()
           END IF
         END IF
         IF ( .NOT. ASSOCIATED( A % BulkMassValues ) ) THEN
           ALLOCATE( A % BulkMassValues( n ) )
         END IF

         DO i=1,n
           A % BulkMassValues(i) = A % MassValues(i)
         END DO
       END IF
     END IF

     IF( PRESENT( BulkDamp ) .AND. ASSOCIATED( A % DampValues) ) THEN
       IF( BulkDamp ) THEN
         n = SIZE( A % DampValues )
         IF( ASSOCIATED( A % BulkDampValues ) ) THEN
           IF( SIZE( A % BulkDampValues ) /= n ) THEN
             DEALLOCATE( A % BulkDampValues ) 
             A % BulkDampValues => NULL()
           END IF
         END IF
         IF ( .NOT. ASSOCIATED( A % BulkDampValues ) ) THEN
           ALLOCATE( A % BulkDampValues( n ) )
         END IF

         DO i=1,n
           A % BulkDampValues(i) = A % DampValues(i)
         END DO
       END IF
     END IF
     
   END SUBROUTINE CopyBulkMatrix
!------------------------------------------------------------------------------


!> Restores the RHS vector and the stiffness, mass and damping matrices by
!> using the data objects A % Bulk* of the argument A
!------------------------------------------------------------------------------
   SUBROUTINE RestoreBulkMatrix( A )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i,n
     
     IF( ASSOCIATED( A % BulkRhs ) ) THEN
       n = SIZE( A % Rhs )
       IF( SIZE( A % BulkRhs ) /= n ) THEN
         CALL Fatal('RestoreBulkMatrix','Cannot restore rhs of different size!')
       END IF
       A % Rhs(1:n) = A % BulkRhs(1:n)
     END IF
     
     IF( ASSOCIATED( A % BulkValues ) ) THEN
       n = SIZE( A % Values )
       IF( SIZE( A % BulkValues ) /= n ) THEN
         CALL Fatal('RestoreBulkMatrix','Cannot restore matrix of different size!')
       END IF
       DO i=1,n
         A % Values(i) = A % BulkValues(i)
       END DO
     END IF

     IF( ASSOCIATED( A % BulkMassValues ) ) THEN
       n = SIZE( A % MassValues )
       IF( SIZE( A % BulkMassValues ) /= n ) THEN
         CALL Fatal('RestoreBulkMatrix','Cannot restore mass matrix of different size!')
       END IF
       DO i=1,n
         A % MassValues(i) = A % BulkMassValues(i)
       END DO
     END IF

     IF( ASSOCIATED( A % BulkDampValues ) ) THEN
       n = SIZE( A % DampValues )
       IF( SIZE( A % BulkDampValues ) /= n ) THEN
         CALL Fatal('RestoreBulkMatrix','Cannot restore damp matrix of different size!')
       END IF
       DO i=1,n
         A % DampValues(i) = A % BulkDampValues(i)
       END DO
     END IF
     
   END SUBROUTINE RestoreBulkMatrix
!------------------------------------------------------------------------------


   

!> Create a child matrix of same toopology but optioanally different size than the
!> parent matrix.
!------------------------------------------------------------------------------

   FUNCTION CreateChildMatrix( ParentMat, ParentDofs, Dofs, ColDofs, CreateRhs, &
       NoReuse, Diagonal ) RESULT ( ChildMat )
     TYPE(Matrix_t) :: ParentMat
     INTEGER :: ParentDofs
     INTEGER :: Dofs
     TYPE(Matrix_t), POINTER :: ChildMat
     INTEGER, OPTIONAL :: ColDofs
     LOGICAL, OPTIONAL :: CreateRhs
     LOGICAL, OPTIONAL :: NoReuse
     LOGICAL, OPTIONAL :: Diagonal
     INTEGER :: i,j,ii,jj,k,l,m,n,nn,Cdofs
     LOGICAL :: ReuseMatrix

     IF( ParentMat % FORMAT /= MATRIX_CRS ) THEN
       CALL Fatal('CreateChildMatrix','Only available for CRS matrix format!')
     END IF

     ChildMat => AllocateMatrix()

     CALL CRS_CreateChildMatrix( ParentMat, ParentDofs, ChildMat, Dofs, ColDofs, CreateRhs, &
         NoReuse, Diagonal )

   END FUNCTION CreateChildMatrix

   
   !> Return number of degrees of freedom and their indexes.
   !------------------------------------------------------------------------------
   FUNCTION mGetElementDOFs( Indexes, UElement, USolver, NotDG )  RESULT(NB)
   !------------------------------------------------------------------------------
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     INTEGER :: Indexes(:)
     LOGICAL, OPTIONAL  ::  NotDG

     TYPE(Solver_t),  POINTER :: Solver
     TYPE(Element_t), POINTER :: Element, Parent, Edge, Face

     LOGICAL :: Found, GB, DGdisable, NeedEdges
     INTEGER :: nb,i,j,k,id, NDOFs, EDOFs, FDOFs, BDOFs,FaceDOFs, EdgeDOFs, BubbleDOFs
     INTEGER :: Ind, ElemFamily, DOFsPerNode

     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF

     NB = 0

     IF (.NOT. ASSOCIATED(Solver)) THEN
       CALL Warn('mGetElementDOFS', 'Cannot return DOFs data without knowing solver')
       RETURN
     END IF

     IF ( PRESENT( UElement ) ) THEN
       Element => UElement
     ELSE
       Element => CurrentModel % CurrentElement
     END IF
     ElemFamily = Element % TYPE % ElementCode / 100

     DGDisable=.FALSE.
     IF (PRESENT(NotDG)) DGDisable=NotDG

     IF ( .NOT. DGDisable .AND. Solver % DG ) THEN
       DO i=1,Element % DGDOFs
         NB = NB + 1
         Indexes(NB) = Element % DGIndexes(i)
       END DO

       IF ( ASSOCIATED( Element % BoundaryInfo ) ) THEN
         IF ( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
           DO i=1,Element % BoundaryInfo % Left % DGDOFs
             NB = NB + 1
             Indexes(NB) = Element % BoundaryInfo % Left % DGIndexes(i)
           END DO
         END IF
         IF ( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
           DO i=1,Element % BoundaryInfo % Right % DGDOFs
             NB = NB + 1
             Indexes(NB) = Element % BoundaryInfo % Right % DGIndexes(i)
           END DO
         END IF
       END IF

       IF ( NB > 0 ) RETURN
     END IF

     id = Element % BodyId
     IF ( Id==0 .AND. ASSOCIATED(Element % BoundaryInfo) ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo % Left) ) &
           id = Element % BoundaryInfo % Left % BodyId

       IF ( ASSOCIATED(Element % BoundaryInfo % Right) ) &
           id = Element % BoundaryInfo % Right % BodyId
     END IF
     IF ( id == 0 ) id=1

     IF (.NOT.ASSOCIATED(Solver % Mesh)) THEN
       IF ( Solver % Def_Dofs(ElemFamily,id,1)>0 ) THEN  
         CALL Warn('mGetElementDOFS', &
             'Solver mesh unknown, the node indices are returned')
         NDOFs = 1
       ELSE
         CALL Warn('mGetElementDOFS', &
             'Solver mesh unknown, no indices returned')
       END IF
     ELSE
       NDOFs = Solver % Mesh % MaxNDOFs
     END IF

     IF ( Solver % Def_Dofs(ElemFamily,id,1)>0 ) THEN
       DOFsPerNode = Element % NDOFs / Element % TYPE % NumberOfNodes
       DO i=1,Element % TYPE % NumberOfNodes
         DO j=1,DOFsPerNode
           NB = NB + 1
           Indexes(NB) = NDOFs * (Element % NodeIndexes(i)-1) + j
         END DO
       END DO
     END IF

     ! The DOFs of advanced elements cannot be returned without knowing mesh
     ! ---------------------------------------------------------------------
     IF (.NOT.ASSOCIATED(Solver % Mesh)) RETURN

     NeedEdges = .FALSE.
     DO i=2,SIZE(Solver % Def_Dofs,3)
       IF (Solver % Def_Dofs(ElemFamily, id, i)>=0) THEN
         NeedEdges = .TRUE.
         EXIT
       END IF
     END DO

     IF (.NOT. NeedEdges) THEN
       !
       ! Check whether face DOFs have been generated by "-quad_face b: ..." or
       ! "-tri_face b: ..."
       !
       IF (ElemFamily == 3 .OR. ElemFamily == 4) THEN
         IF (Solver % Def_Dofs(6+ElemFamily, id, 5)>=0) NeedEdges = .TRUE.
       END IF
     END IF

     IF ( .NOT. NeedEdges ) RETURN

     FaceDOFs   = Solver % Mesh % MaxFaceDOFs
     EdgeDOFs   = Solver % Mesh % MaxEdgeDOFs
     BubbleDOFs = Solver % Mesh % MaxBDOFs

     IF ( ASSOCIATED(Element % EdgeIndexes) ) THEN
       DO j=1,Element % TYPE % NumberOFEdges
         EDOFs = Solver % Mesh % Edges(Element % EdgeIndexes(j)) % BDOFs
         DO i=1,EDOFs
           NB = NB + 1
           Indexes(NB) = EdgeDOFs*(Element % EdgeIndexes(j)-1) + &
               i + NDOFs * Solver % Mesh % NumberOfNodes
         END DO
       END DO
     END IF

     IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
       DO j=1,Element % TYPE % NumberOFFaces
         FDOFs = Solver % Mesh % Faces( Element % FaceIndexes(j) ) % BDOFs
         DO i=1,FDOFs
           NB = NB + 1
           Indexes(NB) = FaceDOFs*(Element % FaceIndexes(j)-1) + i + &
               NDOFs * Solver % Mesh % NumberOfNodes + &
               EdgeDOFs*Solver % Mesh % NumberOfEdges
         END DO
       END DO
     END IF

     GB = Solver % GlobalBubbles
     !ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
     !IF (.NOT. Found) GB = .TRUE.

     IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
       Parent => Element % BoundaryInfo % Left
       IF (.NOT.ASSOCIATED(Parent) ) &
           Parent => Element % BoundaryInfo % Right
       IF (.NOT.ASSOCIATED(Parent) ) RETURN

       SELECT CASE(ElemFamily)
       CASE(2)
         IF ( ASSOCIATED(Parent % EdgeIndexes) ) THEN
           IF ( isActivePElement(Element) ) THEN
             Ind=Element % PDefs % LocalNumber
           ELSE
             DO Ind=1,Parent % TYPE % NumberOfEdges
               Edge => Solver % Mesh % Edges(Parent % EdgeIndexes(ind))
               k = 0
               DO i=1,Edge % TYPE % NumberOfNodes
                 DO j=1,Element % TYPE % NumberOfNodes
                   IF ( Edge % NodeIndexes(i)==Element % NodeIndexes(j) ) k=k+1
                 END DO
               END DO
               IF ( k==Element % TYPE % NumberOfNodes) EXIT
             END DO
           END IF

           EDOFs = Element % BDOFs
           DO i=1,EDOFs
             NB = NB + 1
             Indexes(NB) = EdgeDOFs*(Parent % EdgeIndexes(Ind)-1) + &
                 i + NDOFs * Solver % Mesh % NumberOfNodes
           END DO
         END IF

       CASE(3,4)
         IF ( ASSOCIATED( Parent % FaceIndexes ) ) THEN
           IF ( isActivePElement(Element) ) THEN
             Ind=Element % PDefs % LocalNumber
           ELSE
             DO Ind=1,Parent % TYPE % NumberOfFaces
               Face => Solver % Mesh % Faces(Parent % FaceIndexes(ind))
               k = 0
               DO i=1,Face % TYPE % NumberOfNodes
                 DO j=1,Element % TYPE % NumberOfNodes
                   IF ( Face % NodeIndexes(i)==Element % NodeIndexes(j)) k=k+1
                 END DO
               END DO
               IF ( k==Face % TYPE % NumberOfNodes) EXIT
             END DO
           END IF

           FDOFs = Element % BDOFs
           DO i=1,FDOFs
             NB = NB + 1
             Indexes(NB) = FaceDOFs*(Parent % FaceIndexes(Ind)-1) + i + &
                 NDOFs * Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges
           END DO
         END IF
       END SELECT
     ELSE IF ( GB ) THEN
       IF ( ASSOCIATED(Element % BubbleIndexes) ) THEN
         DO i=1,Element % BDOFs
           NB = NB + 1
           Indexes(NB) = FaceDOFs*Solver % Mesh % NumberOfFaces + &
               NDOFs * Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges + &
               Element % BubbleIndexes(i)
         END DO
       END IF
     END IF
!------------------------------------------------------------------------------
   END FUNCTION mGetElementDOFs
!------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
!> Eliminates bubble degrees of freedom from a local linear system.
!> This version is suitable for flow models with velocity and pressure as 
!> unknowns.
!------------------------------------------------------------------------------
SUBROUTINE NSCondensate( N, Nb, dim, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb, dim
    REAL(KIND=dp) :: K(:,:), F(:)
    REAL(KIND=dp), OPTIONAL :: F1(:)

    REAL(KIND=dp) :: Kbb(nb*dim,nb*dim)
    REAL(KIND=dp) :: Kbl(nb*dim,n*(dim+1)), Klb(n*(dim+1),nb*dim), Fb(nb*dim)

    INTEGER :: m, i, j, l, p, Cdofs((dim+1)*n), Bdofs(dim*nb)

    m = 0
    DO p = 1,n
      DO i = 1,dim+1
        m = m + 1
        Cdofs(m) = (dim+1)*(p-1) + i
      END DO
    END DO

    m = 0
    DO p = 1,nb
      DO i = 1,dim
        m = m + 1
        Bdofs(m) = (dim+1)*(p-1) + i + n*(dim+1)
      END DO
    END DO

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Cdofs)
    Klb = K(Cdofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb*dim )

    F(1:(dim+1)*n) = F(1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
    K(1:(dim+1)*n,1:(dim+1)*n) = &
    K(1:(dim+1)*n,1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb,Kbl ) )

    IF (PRESENT(F1)) THEN
      Fb  = F1(Bdofs)
      F1(1:(dim+1)*n) = F1(1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
    END IF
!------------------------------------------------------------------------------
END SUBROUTINE NSCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Subroutine for the static condensation of element bubbles when there are
!> as many bubbles as DOFs left in the matrix (historically this convention
!> was used; now the count of elementwise bubble functions can be chosen
!> flexibly and then the subroutine CondensateP should be called instead).
!------------------------------------------------------------------------------
SUBROUTINE Condensate( N, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N
    REAL(KIND=dp) :: K(:,:),F(:)
    REAL(KIND=dp), OPTIONAL :: F1(:)
!------------------------------------------------------------------------------    
    IF ( PRESENT(F1) ) THEN
      CALL CondensateP( N, N, K, F, F1 )
    ELSE
      CALL CondensateP( N, N, K, F )
    END IF
!------------------------------------------------------------------------------
END SUBROUTINE Condensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Subroutine for condensation of p element bubbles from linear problem.
!> Modifies given stiffness matrix and force vector(s) 
!------------------------------------------------------------------------------
SUBROUTINE CondensatePR( N, Nb, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N               !< Sum of nodal, edge and face degrees of freedom.
    INTEGER :: Nb              !< Sum of internal (bubble) degrees of freedom.
    REAL(KIND=dp) :: K(:,:)    !< Local stiffness matrix.
    REAL(KIND=dp) :: F(:)      !< Local force vector.
    REAL(KIND=dp), OPTIONAL :: F1(:)  !< Local second force vector.
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Kbb(Nb,Nb), Kbl(Nb,N), Klb(N,Nb), Fb(Nb)
    INTEGER :: i, Ldofs(N), Bdofs(Nb)

    IF ( nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    IF (PRESENT(F1)) THEN
      Fb  = F1(Bdofs)
      F1(1:n) = F1(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    END IF

    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
END SUBROUTINE CondensatePR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Subroutine for condensation of p element bubbles from complex-valued linear 
!> problem. Modifies given stiffness matrix and force vector(s) 
!------------------------------------------------------------------------------
SUBROUTINE CondensatePC( N, Nb, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N               !< Sum of nodal, edge and face degrees of freedom.
    INTEGER :: Nb              !< Sum of internal (bubble) degrees of freedom.
    COMPLEX(KIND=dp) :: K(:,:)    !< Local stiffness matrix.
    COMPLEX(KIND=dp) :: F(:)      !< Local force vector.
    COMPLEX(KIND=dp), OPTIONAL :: F1(:)  !< Local second force vector.
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: Kbb(Nb,Nb), Kbl(Nb,N), Klb(N,Nb), Fb(Nb)
    INTEGER :: i, Ldofs(N), Bdofs(Nb)

    IF ( nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL ComplexInvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    IF (PRESENT(F1)) THEN
      Fb  = F1(Bdofs)
      F1(1:n) = F1(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    END IF

    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE CondensatePC
!------------------------------------------------------------------------------


END MODULE MatrixAssembly
