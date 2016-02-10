!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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

!> Basic utilities used by individual solvers. 
!------------------------------------------------------------------------------

!> \ingroup ElmerLib
!> \{


MODULE SolverUtils

#include "../config.h"

#ifdef USE_ISO_C_BINDINGS
   USE LoadMod
#endif
   USE DirectSolve
   USE Multigrid
   USE IterSolve
   USE ElementUtils
   USE TimeIntegrate
   USE ModelDescription
   USE MeshUtils
   USE ParallelUtils
   USE ParallelEigenSolve
   USE ListMatrix
   USE CRSMatrix

   IMPLICIT NONE

   CHARACTER(LEN=MAX_NAME_LEN), PRIVATE :: NormalTangentialName
   INTEGER, PRIVATE :: NormalTangentialNOFNodes
   INTEGER, POINTER, PRIVATE :: NTelement(:,:)
   LOGICAL, POINTER, PRIVATE :: NTzeroing_done(:,:)
   INTEGER, POINTER, PRIVATE :: BoundaryReorder(:)
   REAL(KIND=dp), POINTER, PRIVATE :: BoundaryNormals(:,:),  &
                                      BoundaryTangent1(:,:), &
                                      BoundaryTangent2(:,:)

   SAVE BoundaryReorder, NormalTangentialNOFNodes, BoundaryNormals, &
              BoundaryTangent1, BoundaryTangent2, NormalTangentialName


CONTAINS

!> Initialize matrix structure and vector to zero initial value.
!------------------------------------------------------------------------------
   SUBROUTINE InitializeToZero( A, ForceVector )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: A  !< Matrix to be initialized
     REAL(KIND=dp) :: ForceVector(:)         !< vector to be initialized
!------------------------------------------------------------------------------
     INTEGER :: i,dim
     LOGICAL :: Found, AnyNT, AnyProj, DoDisplaceMesh
     TYPE(Solver_t), POINTER :: Solver
!------------------------------------------------------------------------------

     IF ( ASSOCIATED( A ) ) THEN
       SELECT CASE( A % FORMAT )
         CASE( MATRIX_CRS )
           CALL CRS_ZeroMatrix( A )

         CASE( MATRIX_BAND,MATRIX_SBAND )
           CALL Band_ZeroMatrix( A )
       END SELECT

       IF ( ASSOCIATED(A % PrecValues) ) THEN
         A % PrecValues(:) = 0._dp 
       END IF

       IF ( ASSOCIATED( A % MassValues ) ) THEN
         A % MassValues(:) = 0.d0
       END IF

       IF ( ASSOCIATED( A % DampValues ) ) THEN
         A % DampValues(:) = 0.d0
       END IF

       IF ( ASSOCIATED( A % Force ) ) THEN
         A % Force(:,1) = 0.0d0
       END IF

       IF ( ASSOCIATED( A % RHS_im ) )  THEN
         A % RHS_im(:) = 0.0d0
       END IF
     END IF

     ForceVector = 0.0d0
     Solver => CurrentModel % Solver

     NormalTangentialNOFNodes = 0
     IF ( Solver % Variable % DOFs <= 1 ) RETURN

     NormalTangentialName = 'Normal-Tangential'
     IF ( SEQL(Solver % Variable % Name, 'flow solution') ) THEN
       NormalTangentialName = TRIM(NormalTangentialName) // ' Velocity'
     ELSE
       NormalTangentialName = TRIM(NormalTangentialName) // ' ' // &
                   GetVarName(Solver % Variable)
     END IF

     AnyNT = ListGetLogicalAnyBC( CurrentModel, NormalTangentialName ) 
     AnyProj =  ListGetLogicalAnyBC( CurrentModel, 'Mortar BC Nonlinear')
     IF( .NOT. (AnyNT .OR. AnyProj ) ) RETURN

     DoDisplaceMesh = ListGetLogical( Solver % Values,'Displace Mesh At Init',Found )
     IF( DoDisplaceMesh ) THEN
       CALL Info('InitializeToZero','Displacing mesh for nonlinear projectors',Level=8)
       CALL DisplaceMesh( Solver % Mesh, Solver % variable % Values, 1, &
           Solver % Variable % Perm, Solver % variable % Dofs )
     END IF

     IF( AnyNT ) THEN
       dim = CoordinateSystemDimension()
       CALL CheckNormalTangentialBoundary( CurrentModel, NormalTangentialName, &
           NormalTangentialNOFNodes, BoundaryReorder, &
           BoundaryNormals, BoundaryTangent1, BoundaryTangent2, dim )
       
       CALL AverageBoundaryNormals( CurrentModel, NormalTangentialName, &
           NormalTangentialNOFNodes, BoundaryReorder, &
           BoundaryNormals, BoundaryTangent1, BoundaryTangent2, &
           dim )
     END IF

     IF( AnyProj ) THEN
       CALL GenerateProjectors(CurrentModel,Solver,Nonlinear = .TRUE. )
     END IF

     IF( DoDisplaceMesh ) THEN
       CALL DisplaceMesh( Solver % Mesh, Solver % variable % Values, -1, &
           Solver % Variable % Perm, Solver % variable % Dofs )
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE InitializeToZero
!------------------------------------------------------------------------------


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

!> Moves a row and and sumes it with the values of a second one, optionally 
!> multiplying with a constant.
!------------------------------------------------------------------------------
   SUBROUTINE MoveRow( A, n1, n2, Coeff )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n1, n2
     REAL(KIND=dp), OPTIONAL :: Coeff
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         IF( PRESENT( Coeff ) ) THEN
           CALL CRS_MoveRow( A,n1,n2,Coeff )
         ELSE
           CALL CRS_MoveRow( A,n1,n2 )  
         END IF

         ! If entries are not found the format is changed on-the-fly
         IF( A % FORMAT == MATRIX_LIST ) THEN
           IF( PRESENT( Coeff ) ) THEN
             CALL CRS_MoveRow( A,n1,n2,Coeff )
           ELSE
             CALL CRS_MoveRow( A,n1,n2 )  
         END IF
       END IF

       CASE( MATRIX_LIST )
         IF( PRESENT( Coeff ) ) THEN
           CALL List_MoveRow( A % ListMatrix,n1,n2,Coeff )
         ELSE
           CALL List_MoveRow( A % ListMatrix,n1,n2 )
         END IF

       CASE DEFAULT
         CALL Warn('MoveRow','Not implemented for this type')
         
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE MoveRow
!------------------------------------------------------------------------------


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


!> Matrix vector multiplication of sparse matrices.
!------------------------------------------------------------------------------
   SUBROUTINE MatrixVectorMultiply( A,u,v )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n
     REAL(KIND=dp), DIMENSION(:) CONTIG :: u,v
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_MatrixVectorMultiply( A,u,v )

     CASE( MATRIX_BAND,MATRIX_SBAND )
       CALL Band_MatrixVectorMultiply( A,u,v )

     CASE( MATRIX_LIST )
       CALL Warn('MatrixVectorMultiply','Not implemented for List matrix type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE MatrixVectorMultiply
!------------------------------------------------------------------------------


!> Matrix vector multiplication of sparse matrices.
!------------------------------------------------------------------------------
   SUBROUTINE TransposeMatrixVectorMultiply( A,u,v )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n
     REAL(KIND=dp), DIMENSION(:) CONTIG :: u,v
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_TransposeMatrixVectorMultiply( A,u,v )

     CASE DEFAULT 
       CALL Fatal('TransposeMatrixVectorMultiply','Not implemented for other than CRS type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE TransposeMatrixVectorMultiply
!------------------------------------------------------------------------------


!> Create a copy of the linear system (Values,Rhs) to (BulkValues,BulkRhs).
!------------------------------------------------------------------------------
   SUBROUTINE CopyBulkMatrix( A, BulkMass )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i,n
     LOGICAL, OPTIONAL :: BulkMass
     
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
     


   END SUBROUTINE CopyBulkMatrix
!------------------------------------------------------------------------------



!> Search faces between passive / non-passive domains; add to boundary
!> elements with given bc-id.
!------------------------------------------------------------------------------
  SUBROUTINE GetPassiveBoundary(Model,Mesh,BcId)
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER :: BcId
    TYPE(Mesh_t) :: Mesh 

    INTEGER, ALLOCATABLE :: arr(:)
    INTEGER :: i,j,n,cnt,ind, sz
    LOGICAL :: L1,L2
    TYPE(Element_t), POINTER :: Faces(:), Telems(:), Face, P1, P2

    CALL FindMeshEdges(Mesh,.FALSE.)
    SELECT CASE(Mesh % MeshDim)
    CASE(2)
      Faces => Mesh % Edges
      n = Mesh % NumberOfEdges
    CASE(3)
      Faces => Mesh % Faces
      n = Mesh % NumberOfFaces
    END SELECT

    ALLOCATE(arr(n)); cnt=0
    DO i=1,n
      P1 => Faces(i) % BoundaryInfo % Right
      P2 => Faces(i) % BoundaryInfo % Left
      IF ( .NOT. ASSOCIATED(P1) .OR. .NOT. ASSOCIATED(P2) ) CYCLE

      L1 = CheckPassiveElement(P1)
      L2 = CheckPassiveElement(P2)

      IF ( L1.NEQV.L2) THEN
        cnt = cnt+1
        arr(cnt) = i
      END IF
    END DO

    sz = Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements - &
             Mesh % PassBCcnt
    IF ( sz+cnt>SIZE(Mesh % Elements) ) THEN
      Telems => Mesh % Elements
      ALLOCATE(Mesh % Elements(sz+cnt))
      IF ( ASSOCIATED(Model % Elements,Telems) ) &
        Model % Elements => Mesh % Elements

      Mesh % Elements(1:sz) = Telems

      ! fix boundary element parent pointers to use new array ...
      ! --------------------------------------------------------
      DO i=1,Mesh % NumberOfBoundaryElements-Mesh % PassBCcnt
        ind = i+Mesh % NumberOfBulkElements
        Face => Mesh % Elements(ind)
        IF ( ASSOCIATED(Face % BoundaryInfo % Left) ) &
          Face % BoundaryInfo % Left  => &
             Mesh % Elements(Face % BoundaryInfo % Left % ElementIndex)
        IF ( ASSOCIATED(Face % BoundaryInfo % Right ) ) &
          Face % BoundaryInfo % Right => &
             Mesh % Elements(Face % BoundaryInfo % Right % ElementIndex)
      END DO

      ! ...likewise for  faces (edges).
      ! -------------------------------
      DO i=1,n
        Face => Faces(i)
        IF ( ASSOCIATED(Face % BoundaryInfo % Left) ) &
          Face % BoundaryInfo % Left  => &
             Mesh % Elements(Face % BoundaryInfo % Left % ElementIndex)
        IF ( ASSOCIATED(Face % BoundaryInfo % Right ) ) &
          Face % BoundaryInfo % Right => &
             Mesh % Elements(Face % BoundaryInfo % Right % ElementIndex)
      END DO

      DEALLOCATE(Telems)
    END IF

    DO i=1,cnt
      sz = sz+1
      Mesh % Elements(sz) = Faces(arr(i))
      Mesh % Elements(sz) % Copy = .TRUE.
      Mesh % Elements(sz) % ElementIndex = sz
      Mesh % Elements(sz) % BoundaryInfo % Constraint = BcId
    END DO
    Mesh % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements - &
                Mesh % PassBCcnt + cnt
    Mesh % PassBCcnt = cnt
    IF ( ASSOCIATED(Model % Elements,Mesh % Elements) ) &
      Model % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
  END SUBROUTINE GetPassiveBoundary
!------------------------------------------------------------------------------

   
!------------------------------------------------------------------------------
!> Check if the current element has been defined passive.
!> This is done by inspecting a looking an the values of "varname Passive"
!> in the Body Force section. It is determined to be passive if it has 
!> more positive than negative hits in an element.
!------------------------------------------------------------------------------
   FUNCTION CheckPassiveElement( UElement )  RESULT( IsPassive )
!------------------------------------------------------------------------------
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     LOGICAL :: IsPassive
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     REAL(KIND=dp), ALLOCATABLE :: Passive(:)
     INTEGER :: body_id, bf_id, nlen, NbrNodes,PassNodes, LimitNodes
     LOGICAL :: Found
     CHARACTER(LEN=MAX_NAME_LEN) :: PassName

     SAVE Passive
!------------------------------------------------------------------------------
     IsPassive = .FALSE.

     IF ( PRESENT( UElement ) ) THEN
       Element => UElement
     ELSE
       Element => CurrentModel % CurrentElement
     END IF

     body_id = Element % BodyId 
     IF ( body_id <= 0 )  RETURN   ! body_id == 0 for boundary elements

     bf_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
         'Body Force', Found, minv=1,maxv=CurrentModel % NumberOfBodyForces )
     IF ( .NOT. Found )  RETURN

     nlen = CurrentModel % Solver % Variable % NameLen
     PassName = GetVarName(CurrentModel % Solver % Variable) // ' Passive'

     IF ( ListCheckPresent(CurrentModel % BodyForces(bf_id) % Values, PassName) ) THEN
       NbrNodes = Element % TYPE % NumberOfNodes
       IF ( ALLOCATED(Passive) ) THEN
         IF ( SIZE(Passive) < NbrNodes ) THEN
           DEALLOCATE(Passive)
           ALLOCATE( Passive(NbrNodes) )
         END IF
       ELSE
         ALLOCATE( Passive(NbrNodes) )
       END IF
       Passive(1:NbrNodes) = ListGetReal( CurrentModel % BodyForces(bf_id) % Values, &
           PassName, NbrNodes, Element % NodeIndexes )
       PassNodes = COUNT(Passive(1:NbrNodes)>0)

       ! Go through the extremum cases first, and if the element is not either fully 
       ! active or passive, then check for some possible given criteria for determining 
       ! the element active / passive. 
       !------------------------------------------------------------------------------
       IF( PassNodes == 0 ) THEN
         CONTINUE
       ELSE IF( PassNodes == NbrNodes ) THEN
         IsPassive = .TRUE.
       ELSE
         LimitNodes = ListGetInteger( CurrentModel % BodyForces(bf_id) % Values, &
             'Passive Element Min Nodes',Found )
         IF( Found ) THEN
           IsPassive = ( PassNodes >= LimitNodes )
         ELSE
           LimitNodes = ListGetInteger( CurrentModel % BodyForces(bf_id) % Values, &
               'Active Element Min Nodes',Found )
           IF( Found ) THEN
             IsPassive = ( PassNodes > NbrNodes - LimitNodes )
           ELSE
             IsPassive = ( 2*PassNodes > NbrNodes )
           END IF
         END IF
       END IF
     END IF
   
!------------------------------------------------------------------------------
   END FUNCTION CheckPassiveElement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  For time dependent simulations add the time derivative coefficient terms
!>  to the local matrix containing other coefficients.
!------------------------------------------------------------------------------
   SUBROUTINE Add1stOrderTime( MassMatrix, StiffMatrix,  &
          Force, dt, n, DOFs, NodeIndexes, Solver, UElement )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: MassMatrix(:,:)   !< Local mass matrix.
     REAL(KIND=dp) :: StiffMatrix(:,:)  !< Local stiffness matrix.
     REAL(KIND=dp) :: Force(:)          !< Local right-hand-side vector.
     REAL(KIND=dp) :: dt                !< Simulation timestep size
     INTEGER :: n                       !< number of element nodes
     INTEGER :: DOFs                    !< variable degrees of freedom
     INTEGER :: NodeIndexes(:)          !< element nodes
     TYPE(Solver_t) :: Solver           !< Solver structure.
     TYPE(Element_t), TARGET, OPTIONAL :: UElement !< Element structure
!------------------------------------------------------------------------------
     LOGICAL :: GotIt
     INTEGER :: i,j,k,l,m,Order
     REAL(KIND=dp) :: s, t
     CHARACTER(LEN=MAX_NAME_LEN) :: Method
     REAL(KIND=dp) :: PrevSol(DOFs*n,Solver % Order), CurSol(DOFs*n), LForce(n*DOFs)
     TYPE(Variable_t), POINTER :: DtVar
     REAL(KIND=dp) :: Dts(Solver % Order)
     LOGICAL :: ConstantDt
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------

     IF ( PRESENT(UElement) ) THEN
       Element => UElement
     ELSE
       Element => CurrentModel % CurrentElement
     END IF

     IF ( Solver % Matrix % Lumped ) THEN
#ifndef OLD_LUMPING
       s = 0.d0
       t = 0.d0
       DO i=1,n*DOFs
         DO j=1,n*DOFs
           s = s + MassMatrix(i,j)
           IF (i /= j) THEN
             MassMatrix(i,j) = 0.d0
           END IF
         END DO
         t = t + MassMatrix(i,i)
       END DO
  
       DO i=1,n
         DO j=1,DOFs
           K = DOFs * (i-1) + j
           L = DOFs * (NodeIndexes(i)-1) + j
           IF ( t /= 0.d0 ) THEN
             MassMatrix(K,K) = MassMatrix(K,K) * s / t
           END IF
         END DO
       END DO
#else
       DO i=1,n*DOFs
         s = 0.0d0
         DO j = 1,n*DOFs
           s = s + MassMatrix(i,j)
           MassMatrix(i,j) = 0.0d0
         END DO
         MassMatrix(i,i) = s
       END DO

       DO i=1,n
         DO j=1,DOFs
           K = DOFs * (i-1) + j
           L = DOFs * (NodeIndexes(i)-1) + j
         END DO
       END DO
#endif
     END IF
!------------------------------------------------------------------------------
     Order = MIN(Solver % DoneTime, Solver % Order)

     DO i=1,n
       DO j=1,DOFs
         K = DOFs * (i-1) + j
         L = DOFs * (NodeIndexes(i)-1) + j
         DO m=1, Order
           PrevSol(K,m) = Solver % Variable % PrevValues(L,m)
         END DO
         CurSol(K) = Solver % Variable % Values(L)
       END DO
     END DO

     LForce(1:n*DOFs) = Force(1:n*DOFs)
     CALL UpdateGlobalForce( Solver % Matrix % Force(:,1), LForce, &
                  n, DOFs, NodeIndexes, UElement=Element )
!------------------------------------------------------------------------------
!PrevSol(:,Order) needed for BDF
     Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )

     SELECT CASE( Method )
     CASE( 'fs' ) 
       CALL FractionalStep( n*DOFs, dt, MassMatrix, StiffMatrix, Force, &
                   PrevSol(:,1), Solver % Beta, Solver )
     CASE('bdf')
       Dts(1) = Dt
       ConstantDt = .TRUE.
       IF(Order > 1) THEN
         DtVar => VariableGet( Solver % Mesh % Variables, 'Timestep size' )
         DO i=2,Order
           Dts(i) = DtVar % PrevValues(1,i-1)
           IF(ABS(Dts(i)-Dts(1)) > 1.0d-6 * Dts(1)) ConstantDt = .FALSE.
         END DO
       END IF

       IF(ConstantDt) THEN
         CALL BDFLocal( n*DOFs, dt, MassMatrix, StiffMatrix, Force, PrevSol, &
                          Order )
       ELSE     
         CALL VBDFLocal( n*DOFs, dts, MassMatrix, StiffMatrix, Force, PrevSol, &
                         Order )
       END IF

     CASE('runge-kutta')
       CALL RungeKutta( n*DOFs, dt, MassMatrix, StiffMatrix, Force, &
                PrevSol(:,1), CurSol )

     CASE DEFAULT
       CALL NewmarkBeta( n*DOFs, dt, MassMatrix, StiffMatrix, Force, &
                 PrevSol(:,1), Solver % Beta )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE Add1stOrderTime
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  For time dependent simulations add the time derivative coefficient terms
!>  to the global matrix containing other coefficients.
!------------------------------------------------------------------------------
   SUBROUTINE Add1stOrderTime_CRS( Matrix, Force, dt, Solver )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: Matrix  !< Global matrix (including stiffness and mass)
     REAL(KIND=dp) :: Force(:)          !< Global right-hand-side vector.
     REAL(KIND=dp) :: dt                !< Simulation timestep size
     TYPE(Solver_t) :: Solver           !< Solver structure.
!------------------------------------------------------------------------------
     LOGICAL :: GotIt
     INTEGER :: i,j,k,l,m,n,Order
     REAL(KIND=dp) :: s, t, msum
     CHARACTER(LEN=MAX_NAME_LEN) :: Method
     TYPE(Variable_t), POINTER :: DtVar
     REAL(KIND=dp) :: Dts(Solver % Order)
     REAL(KIND=dp), POINTER :: PrevSol(:,:), ML(:), CurrSol(:)
     INTEGER, POINTER :: Rows(:), Cols(:)
     LOGICAL :: ConstantDt, Lumped, Found
!------------------------------------------------------------------------------

     CALL Info('Add1stOrderTime_CRS','Adding time discretization to CRS matrix')

!------------------------------------------------------------------------------
     Order = MIN(Solver % DoneTime, Solver % Order)
     Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )
     CurrSol => Solver % Variable % Values
     PrevSol => Solver % Variable % PrevValues

     SELECT CASE( Method )

     CASE( 'fs' ) 
       CALL FractionalStep_CRS( dt, Matrix, Force, PrevSol(:,1), Solver )

     CASE('bdf')
       Dts(1) = Dt
       ConstantDt = .TRUE.
       IF(Order > 1) THEN
         DtVar => VariableGet( Solver % Mesh % Variables, 'Timestep size' )
         DO i=2,Order
           Dts(i) = DtVar % PrevValues(1,i-1)
           IF(ABS(Dts(i)-Dts(1)) > 1.0d-6 * Dts(1)) ConstantDt = .FALSE.
         END DO
       END IF

       IF(ConstantDt) THEN
         CALL BDF_CRS( dt, Matrix, Force, PrevSol, Order )
       ELSE     
         CALL VBDF_CRS( dts, Matrix, Force, PrevSol, Order )
       END IF

     CASE('runge-kutta')
       CALL RungeKutta_CRS( dt, Matrix, Force, PrevSol(:,1), CurrSol )

     CASE DEFAULT
       CALL NewmarkBeta_CRS( dt, Matrix, Force, PrevSol(:,1), &
             Solver % Beta )

     END SELECT

!------------------------------------------------------------------------------
   END SUBROUTINE Add1stOrderTime_CRS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  For time dependent simulations add the time derivative coefficient terms
!>  to the matrix containing other coefficients.
!------------------------------------------------------------------------------
   SUBROUTINE Add2ndOrderTime( MassMatrix, DampMatrix, StiffMatrix,  &
         Force, dt, n, DOFs, NodeIndexes, Solver )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: MassMatrix(:,:)   !< Local mass matrix.
     REAL(KIND=dp) :: DampMatrix(:,:)   !< Local damping matrix.
     REAL(KIND=dp) :: StiffMatrix(:,:)  !< Local stiffness matrix.
     REAL(KIND=dp) :: Force(:)          !< Local right-hand-side vector.
     REAL(KIND=dp) :: dt                !< Simulation timestep size
     INTEGER :: n                       !< number of element nodes
     INTEGER :: DOFs                    !< variable degrees of freedom
     INTEGER :: NodeIndexes(:)          !< element nodes
     TYPE(Solver_t) :: Solver           !< Solver structure.
!------------------------------------------------------------------------------
     LOGICAL :: GotIt
     INTEGER :: i,j,k,l
     CHARACTER(LEN=MAX_NAME_LEN) :: Method
     REAL(KIND=dp) :: s,t
     REAL(KIND=dp) :: X(DOFs*n),V(DOFs*N),A(DOFs*N),LForce(n*DOFs)

!------------------------------------------------------------------------------

     IF ( Solver % Matrix % Lumped ) THEN
!------------------------------------------------------------------------------
#ifndef OLD_LUMPING
       s = 0.d0
       t = 0.d0
       DO i=1,n*DOFs
         DO j=1,n*DOFs
           s = s + MassMatrix(i,j)
           IF (i /= j) THEN
             MassMatrix(i,j) = 0.d0
           END IF
         END DO
         t = t + MassMatrix(i,i)
       END DO

       DO i=1,n
         DO j=1,DOFs
           K = DOFs * (i-1) + j
           IF ( t /= 0.d0 ) THEN
             MassMatrix(K,K) = MassMatrix(K,K) * s / t
           END IF
         END DO
       END DO

       s = 0.d0
       t = 0.d0
       DO i=1,n*DOFs
         DO j=1,n*DOFs
           s = s + DampMatrix(i,j)
           IF (i /= j) THEN
             DampMatrix(i,j) = 0.d0
           END IF
         END DO
         t = t + DampMatrix(i,i)
       END DO

       DO i=1,n
         DO j=1,DOFs
           K = DOFs * (i-1) + j
           IF ( t /= 0.d0 ) THEN
             DampMatrix(K,K) = DampMatrix(K,K) * s / t
           END IF
         END DO
       END DO
#else
!------------------------------------------------------------------------------
!      Lump the second order time derivative terms ...
!------------------------------------------------------------------------------
       DO i=1,n*DOFs
         s = 0.0D0
         DO j=1,n*DOFs
           s = s + MassMatrix(i,j)
           MassMatrix(i,j) = 0.0d0
         END DO
         MassMatrix(i,i) = s
       END DO

!------------------------------------------------------------------------------
!      ... and the first order terms.
!------------------------------------------------------------------------------
       DO i=1,n*DOFs
         s = 0.0D0
         DO j=1,n*DOFs
           s = s + DampMatrix(i,j)
           DampMatrix(i,j) = 0.0d0
         END DO
         DampMatrix(i,i) = s
       END DO
#endif
!------------------------------------------------------------------------------
     END IF
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get previous solution vectors and update current force
!-----------------------------------------------------------------------------
     DO i=1,n
       DO j=1,DOFs
         K = DOFs * (i-1) + j
         IF ( NodeIndexes(i) > 0 ) THEN
           L = DOFs * (NodeIndexes(i)-1) + j
           SELECT CASE(Method)
           CASE DEFAULT
             X(K) = Solver % Variable % PrevValues(L,3)
             V(K) = Solver % Variable % PrevValues(L,4)
             A(K) = Solver % Variable % PrevValues(L,5)
           END SELECT
         END IF
       END DO
     END DO

     LForce(1:n*DOFs) = Force(1:n*DOFs)
     CALL UpdateGlobalForce( Solver % Matrix % Force(:,1), LForce, &
                  n, DOFs, NodeIndexes )
!------------------------------------------------------------------------------
     Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )
     SELECT CASE(Method)
     CASE DEFAULT
       CALL Bossak2ndOrder( n*DOFs, dt, MassMatrix, DampMatrix, StiffMatrix, &
                    Force, X, V, A, Solver % Alpha )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE Add2ndOrderTime
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Update the right-hand-side of the global equation by adding the local entry. 
!------------------------------------------------------------------------------
   SUBROUTINE UpdateTimeForce( StiffMatrix, &
           ForceVector, LocalForce, n, NDOFs, NodeIndexes )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix  !< Global stiffness matrix.
     REAL(KIND=dp) :: LocalForce(:)     !< Local right-hand-side vector.
     REAL(KIND=dp) :: ForceVector(:)    !< Global right-hand-side vector.
     INTEGER :: n                       !< number of element nodes
     INTEGER :: nDOFs                   !< variable degrees of freedom
     INTEGER :: NodeIndexes(:)          !< Element node to global node numbering mapping.
!------------------------------------------------------------------------------
     INTEGER :: i,j,k
!------------------------------------------------------------------------------
     CALL UpdateGlobalForce( StiffMatrix % Force(:,1), LocalForce, &
                     n, NDOFs, NodeIndexes )
     LocalForce = 0.0d0
!------------------------------------------------------------------------------
   END SUBROUTINE UpdateTimeForce
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Add element local matrices & vectors to global matrices and vectors.
!------------------------------------------------------------------------------
   SUBROUTINE UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
      ForceVector, LocalForce, n, NDOFs, NodeIndexes, RotateNT, UElement )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix  !< The global matrix
     REAL(KIND=dp) :: LocalStiffMatrix(:,:)  !< Local matrix to be added to the global matrix.
     REAL(KIND=dp) :: LocalForce(:)          !< Element local force vector.
     REAL(KIND=dp) :: ForceVector(:)         !< The global RHS vector.
     INTEGER :: n                            !< Number of nodes.
     INTEGER :: NDOFs                        !< Number of element nodes. 
     INTEGER :: NodeIndexes(:)               !< Element node to global node numbering mapping.
     LOGICAL, OPTIONAL :: RotateNT           !< Should the global equation be done in local normal-tangential coordinates.
     TYPE(Element_t), OPTIONAL, TARGET :: UElement !< Element to be updated
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,dim, Indexes(n)
     LOGICAL :: Rotate
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
!    Update global matrix and rhs vector....
!------------------------------------------------------------------------------
     IF (PRESENT(UElement)) THEN
        Element => UElement
     ELSE
        Element => CurrentModel % CurrentElement
     END IF
!------------------------------------------------------------------------------
!    Check first if this element has been defined passive
!------------------------------------------------------------------------------
     IF ( CheckPassiveElement(Element) )  RETURN

!------------------------------------------------------------------------------
     Rotate = .TRUE.
     IF ( PRESENT(RotateNT) ) Rotate = RotateNT

     dim = CoordinateSystemDimension()	
     IF ( Rotate .AND. NormalTangentialNOFNodes > 0 .AND. ndofs>=dim) THEN
       Indexes = 0
       Indexes(1:Element % TYPE % NumberOfNodes) = &
             BoundaryReorder(Element % NodeIndexes)
       CALL RotateMatrix( LocalStiffMatrix, LocalForce, n, dim, NDOFs, &
          Indexes, BoundaryNormals, BoundaryTangent1, BoundaryTangent2 )
     END IF
!------------------------------------------------------------------------------
     IF ( ASSOCIATED( StiffMatrix ) ) THEN
       SELECT CASE( StiffMatrix % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_GlueLocalMatrix( StiffMatrix,n,NDOFs,NodeIndexes, &
                          LocalStiffMatrix )

       CASE( MATRIX_LIST )
         CALL List_GlueLocalMatrix( StiffMatrix % ListMatrix,n,NDOFs,NodeIndexes, &
                          LocalStiffMatrix )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_GlueLocalMatrix( StiffMatrix,n,NDOFs,NodeIndexes, &
                          LocalStiffMatrix )
       END SELECT
     END IF

     DO i=1,n
       IF ( Nodeindexes(i) > 0 ) THEN
         DO j=1,NDOFs
           k = NDOFs * (NodeIndexes(i)-1) + j
!$omp atomic
           ForceVector(k) = ForceVector(k) + LocalForce(NDOFs*(i-1)+j)
         END DO
       END IF
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE UpdateGlobalEquations
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Update the global vector with the local vector entry.
!------------------------------------------------------------------------------
   SUBROUTINE UpdateGlobalForce(ForceVector, LocalForce, n, &
             NDOFs, NodeIndexes, RotateNT, UElement )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: LocalForce(:)          !< Element local force vector.
     REAL(KIND=dp) :: ForceVector(:)         !< The global RHS vector.
     INTEGER :: n                            !< Number of nodes.
     INTEGER :: NDOFs                        !< Number of element nodes. 
     INTEGER :: NodeIndexes(:)               !< Element node to global node numbering mapping.
     LOGICAL, OPTIONAL :: RotateNT           !< Should the global equation be done in local normal-tangential coordinates.
     TYPE(Element_t), OPTIONAL, TARGET :: UElement !< Element to be updated
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k, dim,indexes(n)
     LOGICAL :: Rotate
     REAL(KIND=dp) :: LocalStiffMatrix(n*NDOFs,n*NDOFs), LForce(n*NDOFs)
!------------------------------------------------------------------------------
!    Update global matrix and rhs vector....
!------------------------------------------------------------------------------
     IF (PRESENT(UElement)) THEN
        Element => UElement
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

     IF ( CheckPassiveElement( Element ) )  RETURN

     Rotate = .TRUE.
     IF ( PRESENT(RotateNT) ) Rotate=RotateNT

     IF ( Rotate .AND. NormalTangentialNOFNodes>0 ) THEN
       dim = CoordinateSystemDimension()
       Indexes = 0
       ! Element => CurrentModel % CurrentElement
       Indexes(1:Element % TYPE % NumberOfNodes) = &
             BoundaryReorder(Element % NodeIndexes)
       CALL RotateMatrix( LocalStiffMatrix, LocalForce, n, dim, NDOFs, &
          Indexes, BoundaryNormals, BoundaryTangent1, BoundaryTangent2 )
     END IF

     DO i=1,n
       IF ( NodeIndexes(i) > 0 ) THEN
         DO j=1,NDOFs
           k = NDOFs * (NodeIndexes(i)-1) + j
!$omp atomic
           ForceVector(k) = ForceVector(k) + LocalForce(NDOFs*(i-1)+j)
         END DO
       END IF
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE UpdateGlobalForce
!------------------------------------------------------------------------------


!> Updates the mass matrix only.
!------------------------------------------------------------------------------
   SUBROUTINE UpdateMassMatrix( StiffMatrix, LocalMassMatrix, &
                  n, NDOFs, NodeIndexes )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix  !< The global matrix
     REAL(KIND=dp) :: LocalMassMatrix(:,:)   !< Local matrix to be added to the global matrix
     INTEGER :: n                            !<  number of nodes in element
     INTEGER :: NDOFs                        !< number of DOFs per node
     INTEGER :: NodeIndexes(:)               !< Element node to global node numbering mapping
!------------------------------------------------------------------------------
     INTEGER :: i,j,k
     REAL(KIND=dp) :: s,t
     REAL(KIND=dp), POINTER  :: SaveValues(:)
!------------------------------------------------------------------------------
!    Check first if this element has been defined passive
!------------------------------------------------------------------------------
     IF ( CheckPassiveElement() )  RETURN

!------------------------------------------------------------------------------
!    Update global matrix and rhs vector....
!------------------------------------------------------------------------------

     IF ( StiffMatrix % Lumped ) THEN
       s = 0.d0
       t = 0.d0
       DO i=1,n*NDOFs
          DO j=1,n*NDOFs
             s = s + LocalMassMatrix(i,j)
             IF (i /= j) LocalMassMatrix(i,j) = 0.0d0
          END DO
          t = t + LocalMassMatrix(i,i)
       END DO

        DO i=1,n*NDOFs
           LocalMassMatrix(i,i) = LocalMassMatrix(i,i) * s / t
        END DO
     END IF

     SaveValues => StiffMatrix % Values
     StiffMatrix % Values => StiffMatrix % MassValues 

     SELECT CASE( StiffMatrix % FORMAT )
        CASE( MATRIX_CRS )
           CALL CRS_GlueLocalMatrix( StiffMatrix, &
                n, NDOFs, NodeIndexes, LocalMassMatrix )

        CASE( MATRIX_LIST )
           CALL List_GlueLocalMatrix( StiffMatrix % ListMatrix, &
                n, NDOFs, NodeIndexes, LocalMassMatrix )

       CASE( MATRIX_BAND,MATRIX_SBAND )
           CALL Band_GlueLocalMatrix( StiffMatrix, &
                n, NDOFs, NodeIndexes, LocalMassMatrix )
     END SELECT

     StiffMatrix % Values => SaveValues
!------------------------------------------------------------------------------
   END SUBROUTINE UpdateMassMatrix
!------------------------------------------------------------------------------
 

!------------------------------------------------------------------------------
!> Determine soft limiters set. This is called after the solution.
!> and can therefore be active only on the 2nd nonlinear iteration round.
!------------------------------------------------------------------------------
   SUBROUTINE DetermineSoftLimiter( Solver )
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver
!-----------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(variable_t), POINTER :: Var, LoadVar, IterV, LimitVar
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k,n,t,ind,dofs, dof, bf, bc, Upper, Removed, Added, &
         ElemFirst, ElemLast, totsize, i2, j2, ind2
     REAL(KIND=dp), POINTER :: FieldValues(:), LoadValues(:), &
         ElemLimit(:),ElemInit(:), ElemActive(:)
     REAL(KIND=dp) :: LimitSign, EqSign, ValEps, LoadEps, val
     INTEGER, POINTER :: FieldPerm(:), NodeIndexes(:)
     LOGICAL :: Found,AnyLimitBC, AnyLimitBF, GotInit, GotActive
     LOGICAL, ALLOCATABLE :: LimitDone(:)
     LOGICAL, POINTER :: LimitActive(:)
     TYPE(ValueList_t), POINTER :: Params, Entity
     CHARACTER(LEN=MAX_NAME_LEN) :: Name, LimitName, InitName, ActiveName
     LOGICAL, ALLOCATABLE :: InterfaceDof(:)
     INTEGER :: ConservativeAfterIters, NonlinIter, CoupledIter
     LOGICAL :: Conservative, ConservativeAdd, ConservativeRemove, &
         DoAdd, DoRemove, DirectionActive, FirstTime

     Model => CurrentModel
     Var => Solver % Variable
     

     ! Check the iterations counts and determine whether this is the first 
     ! time with this solver. 
     !------------------------------------------------------------------------
     FirstTime = .TRUE.
     iterV => VariableGet( Solver % Mesh % Variables,'nonlin iter')
     IF( ASSOCIATED( iterV ) ) THEN
       NonlinIter =  NINT( iterV % Values(1) ) 
       IF( NonlinIter > 1 ) FirstTime = .FALSE.
     END IF

     iterV => VariableGet( Solver % Mesh % Variables,'coupled iter')
     IF( ASSOCIATED( iterV ) ) THEN
       CoupledIter = NINT( iterV % Values(1) )
       IF( CoupledIter > 1 ) FirstTime = .FALSE.
     END IF
          
     ! Determine variable for computing the contact load used to determine the 
     ! soft limit set.
     !------------------------------------------------------------------------
     CALL Info('DetermineSoftLimiter','Determining soft limiter problems',Level=8)
     LoadVar => VariableGet( Model % Variables, &
         GetVarName(Var) // ' Contact Load',ThisOnly = .TRUE. )
     CALL CalculateLoads( Solver, Solver % Matrix, Var % Values, Var % DOFs, .FALSE., LoadVar ) 

     IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
       CALL Fatal('DetermineSoftLimiter', &
           'No Loads associated with variable '//GetVarName(Var) )
       RETURN
     END IF
     LoadValues => LoadVar % Values


     ! The variable to be constrained by the soft limiters
     FieldValues => Var % Values
     FieldPerm => Var % Perm
     totsize = SIZE( FieldValues )
     dofs = Var % Dofs
     Params => Solver % Values

     ConservativeAdd = .FALSE.
     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Add After Iterations',Conservative ) 
     IF( Conservative ) THEN
       ConservativeAdd = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeAdd ) THEN
         CALL Info('DetermineSoftLimiter','Adding dofs in conservative fashion',Level=8)
       END IF
     END IF

     ConservativeRemove = .FALSE.
     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Remove After Iterations',Found ) 
     IF( Found ) THEN
       Conservative = .TRUE.  
       ConservativeRemove = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeRemove ) THEN
         CALL Info('DetermineSoftLimiter','Removing dofs in conservative fashion',Level=8)
       END IF
     END IF

     LoadEps = ListGetConstReal(Params,'Limiter Load Tolerance',Found ) 
     IF(.NOT. Found ) LoadEps = EPSILON( LoadEps )
         
     ValEps = ListGetConstReal(Params,'Limiter Value Tolerance',Found ) 
     IF(.NOT. Found ) ValEps = EPSILON( ValEps )

     ! The user may want to toggle the sign for various kinds of equations
     ! The default sign that come from standard formulation of Laplace equation.
     !---------------------------------------------------------------------------       
     IF( ListGetLogical( Params,'Limiter Load Sign Negative',Found) ) THEN
       EqSign = -1.0_dp
     ELSE
       EqSign = 1.0_dp
     END IF

     ! Loop through upper and lower limits     
     !------------------------------------------------------------------------
     DO Upper=0,1

       DirectionActive = .FALSE.

       ! If we have both upper and lower limiter then these logical vectors need to be 
       ! reinitialized for the 2nd sweep.
       IF( ALLOCATED( LimitDone) ) LimitDone = .FALSE.
       IF( ALLOCATED( InterfaceDof ) ) InterfaceDof = .FALSE. 

       ! Upper and lower limits have different sign for testing
       !----------------------------------------------------------------------       
       IF( Upper == 0 ) THEN
         LimitSign = -EqSign
       ELSE
         LimitSign = EqSign
       END IF       
       
       ! Go through the components of the field, if many
       !-------------------------------------------------
       DO DOF = 1,dofs

         name = Var % name
         IF ( Var % DOFs > 1 ) name = ComponentName(name,DOF)

         ! The keywords for the correct lower or upper limit of the variable
         !------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           LimitName = TRIM(name)//' Lower Limit'           
           InitName = TRIM(name)//' Lower Initial'
           ActiveName = TRIM(name)//' Lower Active'
         ELSE
           LimitName = TRIM(name)//' Upper Limit' 
           InitName = TRIM(name)//' Upper Initial' 
           ActiveName = TRIM(name)//' Upper Active' 
         END IF

         AnyLimitBC = ListCheckPresentAnyBC( Model, LimitName )
         AnyLimitBF = ListCheckPresentAnyBodyForce( Model, LimitName )

         ! If there is no active keyword then there really is nothing to do
         !----------------------------------------------------------------
         IF( .NOT. ( AnyLimitBC .OR. AnyLimitBF ) ) CYCLE
         DirectionActive = .TRUE.
         
         CALL Info('DetermineSoftLimiter','Applying limit: '//TRIM(LimitName),Level=8)

         ! OK: Do contact for a particular dof and only upper or lower limit
         !------------------------------------------------------------------------

         ! Define the range of elements for which the limiters are active
         !---------------------------------------------------------------
         ElemFirst = Model % NumberOfBulkElements + 1           
         ElemLast = Model % NumberOfBulkElements 
        
         IF( AnyLimitBF ) ElemFirst = 1
         IF( AnyLimitBC ) ElemLast = Model % NumberOfBulkElements + &
             Model % NumberOfBoundaryElements 
         
         IF(.NOT. ALLOCATED( LimitDone) ) THEN
           n = Model % MaxElementNodes
           ALLOCATE( LimitDone( totsize ), ElemLimit(n), ElemInit(n), ElemActive(n) )
           LimitDone = .FALSE.
         END IF

         ! Check that active set vectors for limiters exist, otherwise allocate
         !---------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           IF( .NOT. ASSOCIATED(Var % LowerLimitActive ) ) THEN
             ALLOCATE( Var % LowerLimitActive( totsize ) )
             Var % LowerLimitActive = .FALSE.
           END IF
           LimitActive => Var % LowerLimitActive
         ELSE
           IF( .NOT. ASSOCIATED( Var % UpperLimitActive ) ) THEN
             ALLOCATE( Var % UpperLimitActive( totsize ) )
             Var % UpperLimitActive = .FALSE.
           END IF
           LimitActive => Var % UpperLimitActive
         END IF
 
         Removed = 0
         Added = 0        
         IF(.NOT. ALLOCATED( LimitDone) ) THEN
           n = Model % MaxElementNodes
           ALLOCATE( LimitDone( totsize ), ElemLimit(n), ElemInit(n), ElemActive(n) )
           LimitDone = .FALSE.
         END IF


         IF( FirstTime ) THEN
           ! In the first time set the initial set 
           !----------------------------------------------------------------------
           DO t = ElemFirst, ElemLast
             
             Element => Model % Elements(t)
             Model % CurrentElement => Element

             n = Element % TYPE % NumberOfNodes
             NodeIndexes => Element % NodeIndexes
             
             Found = .FALSE.
             IF( t > Model % NumberOfBulkElements ) THEN
               DO bc = 1,Model % NumberOfBCs
                 IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
                   Found = .TRUE.
                   Entity => Model % BCs(bc) % Values
                   EXIT
                 END IF
               END DO
               IF(.NOT. Found ) CYCLE
             ELSE             
               bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                   'Body Force', Found)
               IF(.NOT. Found ) CYCLE
               Entity => Model % BodyForces(bf) % Values
             END IF
             
             ElemLimit(1:n) = ListGetReal( Entity, &
                 LimitName, n, NodeIndexes, Found)             
             IF(.NOT. Found) CYCLE

             ElemInit(1:n) = ListGetReal( Entity, &
                 InitName, n, NodeIndexes, GotInit)
             ElemActive(1:n) = ListGetReal( Entity, &
                 ActiveName, n, NodeIndexes, GotActive)
             IF(.NOT. ( GotInit .OR. GotActive ) ) CYCLE


             DO i=1,n
               j = FieldPerm( NodeIndexes(i) )
               IF( j == 0 ) CYCLE
               ind = Dofs * ( j - 1) + Dof

               IF( LimitDone(ind) ) CYCLE
             
               ! Go through the active set and free nodes with wrong sign in contact force
               !--------------------------------------------------------------------------       
               IF( GotInit .AND. ElemInit(i) > 0.0_dp ) THEN
                 added = added + 1
                 LimitActive(ind) = .TRUE.
               ELSE IF( GotActive .AND. ElemActive(i) > 0.0_dp ) THEN
                 added = added + 1
                 LimitActive(ind) = .TRUE.
               ELSE
                 LimitActive(ind) = .FALSE.
               END IF

               ! Enforce the values to limits because nonlinear material models
               ! may otherwise lead to divergence of the iteration
               !--------------------------------------------------------------
               IF( LimitActive(ind) ) THEN
                 IF( Upper == 0 ) THEN
                   Var % Values(ind) = MAX( val, ElemLimit(i) )
                 ELSE
                   Var % Values(ind) = MIN( val, ElemLimit(i) )
                 END IF
               END IF
               
               LimitDone(ind) = .TRUE.             
             END DO
           END DO

           CYCLE
         END IF


         IF( Conservative ) THEN
           IF(.NOT. ALLOCATED( InterfaceDof ) ) THEN
             ALLOCATE( InterfaceDof( totsize ) )
             InterfaceDof = .FALSE. 
           END IF

           
           ! Mark limited and unlimited neighbours and thereby make a 
           ! list of interface dofs. 
           !----------------------------------------------------------------------
           DO t = ElemFirst, ElemLast
             
             Element => Model % Elements(t)
             Model % CurrentElement => Element
             n = Element % TYPE % NumberOfNodes
             NodeIndexes => Element % NodeIndexes
             
             Found = .FALSE.
             IF( t > Model % NumberOfBulkElements ) THEN
               DO bc = 1,Model % NumberOfBCs
                 IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
                   Found = .TRUE.
                   Entity => Model % BCs(bc) % Values
                   EXIT
                 END IF
               END DO
               IF(.NOT. Found ) CYCLE
             ELSE             
               bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                   'Body Force', Found)
               IF(.NOT. Found ) CYCLE
               Entity => Model % BodyForces(bf) % Values
             END IF

             ElemLimit(1:n) = ListGetReal( Entity, &
                 LimitName, n, NodeIndexes, Found)
             IF(.NOT. Found) CYCLE

             DO i=1,n
               j = FieldPerm( NodeIndexes(i) )
               IF( j == 0 ) CYCLE
               ind = Dofs * ( j - 1) + Dof
               
               DO i2 = i+1,n
                 j2 = FieldPerm( NodeIndexes(i2) )
                 IF( j2 == 0 ) CYCLE
                 ind2 = Dofs * ( j2 - 1) + Dof
                 
                 IF( LimitActive(ind) .NEQV. LimitActive(ind2) ) THEN
                   InterfaceDof(ind) = .TRUE.
                   InterfaceDof(ind2) = .TRUE.
                 END IF
               END DO
             END DO
           END DO

           CALL Info('DetermineSoftLimiter',&
               'Number of interface dofs: '//TRIM(I2S(COUNT(InterfaceDof))),Level=8)
         END IF

       
         ! Add and release dofs from the contact set:
         ! If it is removed it cannot be added. 
         !----------------------------------------------------------------------
         DO t = ElemFirst, ElemLast

           Element => Model % Elements(t)
           Model % CurrentElement => Element
           n = Element % TYPE % NumberOfNodes
           NodeIndexes => Element % NodeIndexes
           
           Found = .FALSE.
           IF( t > Model % NumberOfBulkElements ) THEN
             DO bc = 1,Model % NumberOfBCs
               IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
                 Found = .TRUE.
                 Entity => Model % BCs(bc) % Values
                 EXIT
               END IF
             END DO
             IF(.NOT. Found ) CYCLE
           ELSE             
             bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                 'Body Force', Found)
             IF(.NOT. Found ) CYCLE
             Entity => Model % BodyForces(bf) % Values
           END IF
           
           ElemLimit(1:n) = ListGetReal( Entity, &
               LimitName, n, NodeIndexes, Found)             
           IF(.NOT. Found) CYCLE
           
           ElemActive(1:n) = ListGetReal( Entity, &
               ActiveName, n, NodeIndexes, GotActive)

           DO i=1,n
             j = FieldPerm( NodeIndexes(i) )
             IF( j == 0 ) CYCLE
             ind = Dofs * ( j - 1) + Dof

             IF( LimitDone(ind) ) CYCLE
             
             ! Go through the active set and free nodes with wrong sign in contact force
             !--------------------------------------------------------------------------       
             IF( GotActive .AND. ElemActive(i) > 0.0_dp ) THEN
               IF(.NOT. LimitActive( ind ) ) THEN
                 added = added + 1
                 LimitActive(ind) = .TRUE. 
               END IF
             ELSE IF( LimitActive( ind ) ) THEN
               DoRemove = ( LimitSign * LoadValues(ind) > LimitSign * LoadEps ) 
               IF( DoRemove ) THEN
                 ! In the conservative mode only release nodes from contact set 
                 ! when they are adjacent to dofs that previously was not in the set.
                 ! This means that set is released only at the boundaries. 
                 IF( ConservativeRemove ) DoRemove = InterfaceDof( ind ) 
                 IF( DoRemove ) THEN
                   removed = removed + 1
                   LimitActive(ind) = .FALSE.
                   CYCLE
                 END IF
               END IF
             ELSE
               ! Go through the dofs that are beyond the contact surface.
               !-----------------------------------------------------------
               val = Var % Values(ind) 
               IF( Upper == 0 ) THEN
                 DoAdd = ( val < ElemLimit(i) - ValEps )
               ELSE
                 DoAdd = ( val > ElemLimit(i) + ValEps )
               END IF
               
               IF( DoAdd ) THEN
                 IF( ConservativeAdd ) DoAdd = InterfaceDof( ind ) 
                 IF( DoAdd ) THEN
                   IF( .NOT. LimitActive(ind) ) THEN
                     added = added + 1
                     LimitActive(ind) = .TRUE.
                   END IF
                 END IF
               END IF
             END IF

             ! Enforce the values to limits because nonlinear material models
             ! may otherwise lead to divergence of the iteration
             !--------------------------------------------------------------
             IF( LimitActive(ind) ) THEN
               IF( Upper == 0 ) THEN
                 Var % Values(ind) = MAX( val, ElemLimit(i) )
               ELSE
                 Var % Values(ind) = MIN( val, ElemLimit(i) )
               END IF
             END IF

             LimitDone(ind) = .TRUE.             
           END DO
         END DO
       END DO

       IF( DirectionActive ) THEN      
         ! Output some information before exiting
         !---------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           CALL Info('DetermineSoftLimiter','Determined lower soft limit set',Level=5)
         ELSE
           CALL Info('DetermineSoftLimiter','Determined upper soft limit set',Level=5)
         END IF

         WRITE(Message,'(A,I0)') 'Number of limited dofs for '&
             //TRIM(GetVarName(Var))//': ',COUNT( LimitActive )
         CALL Info('DetermineSoftLimiter',Message,Level=5)
         
         IF(added >= 0) THEN
           WRITE(Message,'(A,I0,A)') 'Added ',added,' dofs to the set'
           CALL Info('DetermineSoftLimiter',Message,Level=5)
         END IF
         
         IF(removed >= 0) THEN
           WRITE(Message,'(A,I0,A)') 'Removed ',removed,' dofs from the set'
           CALL Info('DetermineSoftLimiter',Message,Level=5)
         END IF
       END IF
     END DO

     ! Optionally save the limiters as a field variable so that 
     ! lower limit is given value -1.0 and upper limit value +1.0.
     IF( ListGetLogical( Params,'Save Limiter',Found ) ) THEN
       LimitVar => VariableGet( Model % Variables, &
           GetVarName(Var) // ' Contact Active',ThisOnly = .TRUE. )
       IF(.NOT. ASSOCIATED( LimitVar ) ) THEN
         CALL VariableAddVector( Model % Variables, Model % Mesh, Solver,&
             GetVarName(Var) //' Contact Active', Perm = FieldPerm )
         LimitVar => VariableGet( Model % Variables, &
             GetVarName(Var) // ' Contact Active',ThisOnly = .TRUE. )
       END IF

       ! Currently the visulized limit is always scalar even though the limited field could be a vector!
       DO i = 1, SIZE( LimitVar % Values ) 
         LimitVar % Values(i) = 0.0_dp
         DO j=1,Var % Dofs
           IF( ASSOCIATED( Var % LowerLimitActive ) ) THEN
             IF( Var % LowerLimitActive(Var%Dofs*(i-1)+j) ) LimitVar % Values(i) = -1.0_dp
           END IF
           IF( ASSOCIATED( Var % UpperLimitActive ) ) THEN
             IF( Var % UpperLimitActive(Var%Dofs*(i-1)+j) ) LimitVar % Values(i) = 1.0_dp
           END IF
         END DO
       END DO
     END IF

     IF( ALLOCATED( LimitDone ) ) THEN
       DEALLOCATE( LimitDone, ElemLimit, ElemInit, ElemActive ) 
     END IF
     
     IF( ALLOCATED( InterfaceDof ) ) THEN
       DEALLOCATE( InterfaceDof )
     END IF

     CALL Info('DetermineSoftLimiter','All done',Level=12)

  END SUBROUTINE DetermineSoftLimiter
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Subroutine for determine the contact set and create the necessary data
!> for setting up the contact conditions. As input the mortar projectors,
!> the current solution, and the stiffness matrix are used.  
!------------------------------------------------------------------------------
   SUBROUTINE DetermineContact( Solver )
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver
!-----------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(variable_t), POINTER :: Var, LoadVar, IterVar
     TYPE(Variable_t), POINTER :: DistVar, NormalLoadVar, SlipLoadVar, VeloVar, &
         WeightVar, NormalActiveVar, StickActiveVar, GapVar
     TYPE(Element_t), POINTER :: Element
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: i,j,k,l,n,m,t,ind,dofs, bf, Upper, &
         ElemFirst, ElemLast, totsize, i2, j2, ind2, bc_ind, master_ind, &
         DistSign, LimitSign, DofN, DofT1, DofT2, Limited, LimitedMin, TimeStep
     REAL(KIND=dp), POINTER :: FieldValues(:), LoadValues(:), ElemLimit(:),pNormal(:,:),&
         RotatedField(:)
     REAL(KIND=dp) :: ValEps, LoadEps, val, ContactNormal(3), &
         ContactT1(3), ContactT2(3), LocalT1(3), LocalT2(3), &
         LocalNormal(3), NodalForce(3), wsum, coeff, &
         Dist, DistN, DistT1, DistT2, NTT(3,3), RotVec(3), dt
     INTEGER, POINTER :: FieldPerm(:), NodeIndexes(:)
     LOGICAL :: Found,AnyLimitBC, AnyLimitBF
     LOGICAL, ALLOCATABLE :: LimitDone(:),InterfaceDof(:)
     LOGICAL, POINTER :: LimitActive(:)
     TYPE(ValueList_t), POINTER :: Params
     CHARACTER(LEN=MAX_NAME_LEN) :: Str, LimitName, VarName, ContactType
     INTEGER :: ConservativeAfterIters, ActiveDirection, NonlinIter, CoupledIter
     LOGICAL :: ConservativeAdd, ConservativeRemove, &
         DoAdd, DoRemove, DirectionActive, Rotated, FlatProjector, PlaneProjector, &
         RotationalProjector, FirstTime = .TRUE., &
         AnyRotatedContact, ThisRotatedContact, StickContact, TieContact, FrictionContact, SlipContact, &
         CalculateVelocity, NodalNormal, ResidualMode, AddDiag, SkipFriction, DoIt
     TYPE(MortarBC_t), POINTER :: MortarBC
     TYPE(Matrix_t), POINTER :: Projector, DualProjector
     TYPE(ValueList_t), POINTER :: BC, MasterBC
     REAL(KIND=dp), POINTER :: nWrk(:,:)
     LOGICAL :: CreateDual

     
     SAVE FirstTime

     Model => CurrentModel
     Var => Solver % Variable
     VarName = GetVarName( Var ) 
     Mesh => Solver % Mesh

     ! Is any boundary rotated or not
     AnyRotatedContact = ( NormalTangentialNOFNodes > 0 ) 

     ! The variable to be constrained by the contact algorithm
     ! Here it is assumed to be some "displacement" i.e. a vector quantity
     FieldValues => Var % Values
     FieldPerm => Var % Perm
     totsize = SIZE( FieldValues )
     dofs = Var % Dofs
     Params => Solver % Values

     IterVar => VariableGet( Model % Variables,'coupled iter')
     CoupledIter = NINT( IterVar % Values(1) )

     IterVar => VariableGet( Model % Variables,'nonlin iter')
     NonlinIter = NINT( IterVar % Values(1) )
     
     IterVar => VariableGet( Model % Variables,'timestep')
     Timestep = NINT( IterVar % Values(1) )

     IterVar => VariableGet( Mesh % Variables,'timestep size')
     IF( ASSOCIATED( IterVar ) ) THEN
       dt = IterVar % Values(1)
     ELSE
       dt = 1.0_dp
     END IF

     !FirstTime = ( NonlinIter == 1 .AND. CoupledIter == 1 ) 

     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Add After Iterations',ConservativeAdd ) 
     IF( ConservativeAdd ) THEN
       IF( CoupledIter == 1 ) ConservativeAdd = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeAdd ) THEN
         CALL Info('DetermineContact','Adding dofs in conservative fashion',Level=8)
       END IF
     END IF

     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Remove After Iterations',ConservativeRemove ) 
     IF( ConservativeRemove ) THEN
       IF( CoupledIter == 1 ) ConservativeRemove = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeRemove ) THEN
         CALL Info('DetermineContact','Removing dofs in conservative fashion',Level=8)
       END IF
     END IF
         
     ResidualMode = ListGetLogical(Params,&
         'Linear System Residual Mode',Found )

     CalculateVelocity = ListGetLogical(Params,&
         'Apply Contact Velocity',Found )
     IF(.NOT. Found ) THEN
       Str = ListGetString( CurrentModel % Simulation, 'Simulation Type' )
       CalculateVelocity =  ( Str == 'transient' ) 
     END IF

     NodalNormal = ListGetLogical(Params,&
         'Use Nodal Normal',Found )

     LoadEps = ListGetConstReal(Params,'Limiter Load Tolerance',Found ) 
     IF(.NOT. Found ) LoadEps = EPSILON( LoadEps )
         
     ValEps = ListGetConstReal(Params,'Limiter Value Tolerance',Found ) 
     IF(.NOT. Found ) ValEps = EPSILON( ValEps )

     IF( .NOT. ASSOCIATED( Model % Solver % MortarBCs ) ) THEN
       CALL Fatal('DetermineContact','Cannot apply contact without projectors!')
     END IF

     ! a) Create rotateted contact if needed
     CALL RotatedDisplacementField() 

     ! b) Create and/or obtain pointers to boundary variables 
     CALL GetContactFields( FirstTime )

     ! c) Calculate the contact loads to the normal direction
     LoadVar => CalculateContactLoad() 
     LoadValues => LoadVar % Values

     
     ! Loop over each contact pair
     !--------------------------------------------------------------
     DO bc_ind = 1, Model % NumberOfBCs
 
       MortarBC => Model % Solver % MortarBCs(bc_ind)  
       IF( .NOT. ASSOCIATED( MortarBC ) ) CYCLE

       Projector => MortarBC % Projector
       IF(.NOT. ASSOCIATED(Projector) ) CYCLE

       BC => Model % BCs(bc_ind) % Values

       CALL Info('DetermineContact','Set contact for boundary: '&
           //TRIM(I2S(bc_ind)),Level=8)
       Model % Solver % MortarBCsChanged = .TRUE.

       FlatProjector = ListGetLogical( BC, 'Flat Projector',Found ) 
       PlaneProjector = ListGetLogical( BC, 'Plane Projector',Found )
       RotationalProjector = ListGetLogical( BC, 'Rotational Projector',Found ) .OR. &
           ListGetLogical( BC, 'Cylindrical Projector',Found )

       ! Is the current boundary rotated or not
       ThisRotatedContact = ListGetLogical( BC,'Normal-Tangential '//TRIM(VarName),Found)

       IF( FlatProjector ) THEN
         ActiveDirection = ListGetInteger( BC, 'Flat Projector Coordinate',Found )
         IF( .NOT. Found ) ActiveDirection = dofs       
       ELSE IF( PlaneProjector ) THEN
         pNormal => ListGetConstRealArray( BC,'Plane Projector Normal',Found)
         IF( ThisRotatedContact ) THEN
           ActiveDirection = 1
         ELSE
           ActiveDirection = 1
           DO i=2,3
             IF( ABS( pnormal(i,1) ) > ABS( pnormal(ActiveDirection,1) ) ) THEN
               ActiveDirection = i
             END IF
           END DO
           CALL Info('DetermineContact','Active direction set to: '//TRIM(I2S(ActiveDirection)))
         END IF
       ELSE IF( RotationalProjector ) THEN
         ActiveDirection = 1
         IF( .NOT. ThisRotatedContact ) THEN
           CALL Warn('DetermineContact','Rotational projector may not work without N-T coordinates')
         END IF
       ELSE
         CALL Fatal('DetermineContact','Projector must be current either flat, plane, cylinder or rotational!')
       END IF
      

       ! Get the pointer to the other side i.e. master boundary  
       master_ind = ListGetInteger( BC,'Mortar BC',Found )
       IF( .NOT. Found ) master_ind = ListGetInteger( BC,'Contact BC',Found )
       MasterBC => Model % BCs(master_ind) % Values
       
       ! If we have dual projector we may use it to map certain quantities directly to master nodes
       DualProjector => Projector % Ematrix
       CreateDual = ASSOCIATED( DualProjector )
       IF( CreateDual ) THEN
         CALL Info('DetermineContact','Using also the dual projector',Level=8)
       END IF
      
       ! If we have N-T system then the mortar condition for the master side
       ! should have reverse sign as both normal displacement diminish the gap.
       IF( ThisRotatedContact ) THEN
         IF( master_ind > 0 ) THEN
           IF( .NOT. ListGetLogical( MasterBC, &
               'Normal-Tangential '//TRIM(VarName),Found) ) THEN
             CALL Fatal('DetermineContact','Master boundary '//TRIM(I2S(master_ind))//&
                 ' should also have N-T coordinates!')
           END IF
         END IF

         CALL Info('DetermineContact','We have a normal-tangential system',Level=6)
         MortarBC % MasterScale = -1.0_dp
         DofN = 1
       ELSE                 
         DofN = ActiveDirection 
       END IF

       ! Get the degrees of freedom related to the normal and tangential directions
       DofT1 = 0; DofT2 = 0
       DO i=1,dofs
         IF( i == DofN ) CYCLE
         IF( DofT1 == 0 ) THEN
           DofT1 = i 
           CYCLE
         END IF
         IF( DofT2 == 0 ) THEN
           DofT2 = i
           CYCLE
         END IF
       END DO

       ! This is the normal that is used to detect the signed distance
       ! and tangent vectors used to detect surface velocity
       IF( PlaneProjector ) THEN
         ContactNormal = pNormal(1:3,1)
       ELSE
         ContactNormal = 0.0_dp
         ContactNormal(ActiveDirection) = 1.0_dp
       END IF
       ContactT1 = 0.0_dp
       ContactT1(DofT1) = 1.0_dp
       ContactT2 = 0.0_dp
       IF(DofT2>0) ContactT2(DofT2) = 1.0_dp

       ! Get the contact type. There are four possibilities currently. 
       ! Only one is active at a time while others are false. 
       StickContact = .FALSE.; TieContact = .FALSE.
       FrictionContact = .FALSE.; SlipContact = .FALSE.

       ContactType = ListGetString( BC,'Contact Type',Found ) 
       IF( Found ) THEN
         SELECT CASE ( ContactType )
         CASE('stick')
           StickContact = .TRUE.
         CASE('tie')
           TieContact = .TRUE.
         CASE('friction')
           FrictionContact = .TRUE.
         CASE('slide')
           SlipContact = .TRUE.
         CASE Default
           CALL Fatal('DetermineContact','Unknown contact type: '//TRIM(ContactType))
         END SELECT
       ELSE
         StickContact = ListGetLogical( BC,'Stick Contact',Found )
         IF(.NOT. Found ) TieContact = ListGetLogical( BC,'Tie Contact',Found )
         IF(.NOT. Found ) FrictionContact = ListGetLogical( BC,'Friction Contact',Found )
         IF(.NOT. Found ) SlipContact = ListGetLogical( BC,'Slide Contact',Found )
         IF(.NOT. Found ) THEN 
           CALL Warn('DetermineContact','No contact type given, assuming > Slip Contact <')
           SlipContact = .TRUE.
         END IF
       END IF

       IF( StickContact ) CALL Info('DetermineContact','Using stick contact for displacement',Level=10)
       IF( TieContact ) CALL Info('DetermineContact','Using tie contact for displacement',Level=10)
       IF( FrictionContact ) CALL Info('DetermineContact','Using friction contact for displacement',Level=10)
       IF( SlipContact ) CALL Info('DetermineContact','Using slip contact for displacement',Level=10)
       

       ! At the first time it may be beneficial to assume frictionless initial contact.
       SkipFriction = .FALSE.
       IF( (FrictionContact .OR. StickContact .OR. SlipContact ) .AND. TimeStep == 1 ) THEN
         DoIt = .NOT. ListGetLogical(BC,'Initial Contact Friction',Found )
         IF( DoIt ) THEN
           FrictionContact = .FALSE.; StickContact = .FALSE.
           SlipContact = .TRUE.
           SkipFriction = .TRUE.
           CALL Info('DetermineContact','Assuming frictionless initial contact',Level=10)
         END IF
       ELSE IF( ( FrictionContact .OR. SlipContact) .AND. NonlinIter == 1 ) THEN
         DoIt = ListGetLogical(BC,'Nonlinear System Initial Stick',Found )
         IF(.NOT. Found ) THEN
           ! If contact velocity is not given then it is difficult to determine the direction at 
           ! start of nonlinear iteration when the initial guess still reflects the old displacements. 
           DoIt = .NOT. ListCheckPresent( BC,'Contact Velocity') .AND. &
               ListCheckPresent( BC,'Dynamic Friction Coefficient')
         END IF
         IF( DoIt ) THEN
           FrictionContact = .FALSE.
           SlipContact = .FALSE.
           StickContact = .TRUE.
           CALL Info('DetermineContact','Assuming sticking in first iteration initial contact',Level=10)
         END IF        
       END IF

       ! If we have stick contact then create a diagonal entry to the projection matrix.
       IF( StickContact .OR. FrictionContact ) THEN
         AddDiag = ListCheckPresent( BC,'Stick Contact Coefficient')      
       ELSE
         AddDiag = .FALSE.
       END IF

       ! d) allocate and initialize all necessary vectors for the contact 
       !------------------------------------------------------------------
       CALL InitializeMortarVectors()
     
       ! e) If the contact set is set up in a conservative fashion we need to mark interface nodes
       !------------------------------------------------------------------
       IF( ConservativeAdd .OR. ConservativeRemove ) THEN
         CALL MarkInterfaceDofs()
       END IF

       ! f) Compute the normal load used to determine whether contact should be released.
       !    Also check the direction to which the signed distance should be computed
       !------------------------------------------------------------------
       CALL CalculateContactPressure()

       ! g) Calculate the distance used to determine whether contact should be added
       !------------------------------------------------------------------
       CALL CalculateMortarDistance()
        
       ! h) Determine the contact set in normal direction
       !------------------------------------------------------------------
       CALL NormalContactSet()

       ! i) If requested ensure a minumum number of contact nodes
       !-------------------------------------------------------------------
       LimitedMin = ListGetInteger( BC,'Contact Active Set Minimum',Found)
       IF( Found ) CALL IncreaseContactSet( LimitedMin )

       ! j) Determine the stick set in tangent direction
       !------------------------------------------------------------------
       CALL TangentContactSet()
       
       ! k) Add the stick coefficient if present
       !------------------------------------------------------------------
       IF( AddDiag ) THEN
         CALL StickCoefficientSet()
       END IF

       ! l) We can map information from slave to master either by creating a dual projector
       !    or using the transpose of the original projector to map field from slave to master.
       !-----------------------------------------------------------------------------------
       IF(.NOT. CreateDual ) THEN
         CALL ProjectFromSlaveToMaster()
       END IF

       ! m) If we have dynamic friction then add it 
       IF( .NOT. SkipFriction .AND. ( SlipContact .OR. FrictionContact ) ) THEN
         CALL SetSlideFriction()
       END IF

       IF( ConservativeAdd .OR. ConservativeRemove ) THEN
         DEALLOCATE( InterfaceDof )
       END IF
     END DO
     
     ! Use N-T coordinate system for the initial guess
     ! This is mandatory if using the residual mode linear solvers 
     IF( AnyRotatedContact ) THEN
       DEALLOCATE( RotatedField ) 
     END IF
     

     FirstTime = .FALSE.
     CALL Info('DetermineContact','All done',Level=10)

   CONTAINS


     ! Given the cartesian solution compute the rotated solution.
     !-------------------------------------------------------------------------
     SUBROUTINE RotatedDisplacementField( ) 

       REAL(KIND=dp) :: RotVec(3)
       INTEGER :: i,j,k,m

       IF( .NOT. AnyRotatedContact ) RETURN

       CALL Info('DetermineContact','Rotating displacement field',Level=8)
       ALLOCATE( RotatedField(Solver % Matrix % NumberOfRows ) )
       RotatedField = Var % Values
       
       DO i=1,Solver % Mesh % NumberOfNodes
         j = Solver % Variable % Perm(i)
         IF( j == 0 ) CYCLE
         m = BoundaryReorder(i)
         IF( m == 0 ) CYCLE
         
         RotVec = 0._dp
         DO k=1,Var % DOFs
           RotVec(k) = RotatedField(Var % DOfs*(j-1)+k)
         END DO
         CALL RotateNTSystem( RotVec, i )
         DO k=1,Var % DOFs
           RotatedField(Var % Dofs*(j-1)+k) = RotVec( k )
         END DO
       END DO

     END SUBROUTINE RotatedDisplacementField


     ! Given the previous solution and the current stiffness matrix 
     ! computes the load normal to the surface i.e. the contact load.
     ! If we have normal-tangential coordinate system then also the load is in 
     ! the same coordinate system. 
     !-------------------------------------------------------------------------
     FUNCTION CalculateContactLoad( ) RESULT ( LoadVar )

       TYPE(Variable_t), POINTER :: LoadVar
       REAL(KIND=dp), POINTER :: TempX(:)
       REAL(KIND=dp) :: RotVec(3)
       INTEGER :: i,j,k,m


       CALL Info('DetermineContact','Determining contact load for contact problems',Level=10)

       LoadVar => VariableGet( Model % Variables, &
           TRIM(VarName) // ' Contact Load',ThisOnly = .TRUE. )
       IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
         CALL Fatal('DetermineContact', &
             'No Loads associated with variable: '//GetVarName(Var) )
       END IF

       IF( AnyRotatedContact ) THEN
         TempX => RotatedField 
       ELSE
         TempX => FieldValues
       END IF

       CALL CalculateLoads( Solver, Solver % Matrix, TempX, Var % DOFs, .FALSE., LoadVar ) 

     END FUNCTION CalculateContactLoad


     ! Create fields where the contact information will be saved.
     ! Create the fields both for slave and master nodes at each 
     ! contact pair.
     !--------------------------------------------------------------
     SUBROUTINE GetContactFields( DoAllocate )

       LOGICAL :: DoAllocate
       INTEGER, POINTER :: BoundaryPerm(:), Indexes(:)
       INTEGER :: i,j,k,t,n
       TYPE(Element_t), POINTER :: Element
       LOGICAL, ALLOCATABLE :: ActiveBCs(:)


       IF( DoAllocate ) THEN
         CALL Info('DetermineContact','Creating contact fields',Level=8)
        
         ALLOCATE( BoundaryPerm(Mesh % NumberOfNodes) )
         BoundaryPerm = 0

         ALLOCATE( ActiveBCs(Model % NumberOfBcs ) )
         ActiveBCs = .FALSE.

         DO i=1,Model % NumberOfBCs 
           j = ListGetInteger( Model % BCs(i) % Values,'Mortar BC',Found ) 
           IF(.NOT. Found ) THEN
             j = ListGetInteger( Model % BCs(i) % Values,'Contact BC',Found ) 
           END IF
           IF( j > 0 ) THEN
             ActiveBCs(i) = .TRUE.
             ActiveBCs(j) = .TRUE. 
           END IF
         END DO

         DO t=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

           Element => Mesh % Elements( t )                 
           DO i = 1, Model % NumberOfBCs
             IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
               IF( ActiveBCs(i) ) THEN
                 BoundaryPerm( Element % NodeIndexes ) = 1
               END IF
             END IF
           END DO
         END DO

         DEALLOCATE( ActiveBCs )

         j = 0
         DO i=1,Mesh % NumberOfNodes
           IF( BoundaryPerm(i) > 0 ) THEN
             j = j + 1
             BoundaryPerm(i) = j
           END IF
         END DO

         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Distance',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Gap',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Normalload',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Slipload',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Weight',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Active',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Stick',1,Perm = BoundaryPerm )
         IF( CalculateVelocity ) THEN
           CALL VariableAddVector( Model % Variables,Mesh,Solver,&
               TRIM(VarName)//' Contact Velocity',Dofs,Perm = BoundaryPerm )
         END IF
       END IF

       DistVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Distance')
       GapVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Gap')
       NormalLoadVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Normalload')
       SlipLoadVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Slipload')
       WeightVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Weight')
       NormalActiveVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Active')
       StickActiveVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Stick') 
       IF( CalculateVelocity ) THEN
         VeloVar => VariableGet( Model % Variables,&
             TRIM(VarName)//' Contact Velocity')
       END IF

       NormalActiveVar % Values = -1.0_dp
       StickActiveVar % Values = -1.0_dp

     END SUBROUTINE GetContactFields



     ! Allocates the vectors related to the mortar contact surface, if needed.
     ! Initialize the mortar vectors and mortar permutation future use. 
     ! As the geometry changes the size of the projectors may also change. 
     !----------------------------------------------------------------------------
     SUBROUTINE InitializeMortarVectors()

       INTEGER :: onesize, totsize
       INTEGER, POINTER :: Perm(:)
       LOGICAL, POINTER :: Active(:)
       REAL(KIND=dp), POINTER :: Diag(:)
       LOGICAL :: SamePerm, SameSize

       
       onesize = Projector % NumberOfRows
       totsize = Dofs * onesize

       IF( .NOT. AddDiag .AND. ASSOCIATED(MortarBC % Diag) ) THEN
         DEALLOCATE( MortarBC % Diag ) 
       END IF


       ! Create the permutation that is later need in putting the diag and rhs to correct position
       ALLOCATE( Perm( SIZE( FieldPerm ) ) )
       Perm = 0
       DO i=1,SIZE( Projector % InvPerm )
         j = Projector % InvPerm(i) 
         IF( j == 0 ) CYCLE
         Perm( j ) = i
       END DO

       ! First time nothing is allocated
       IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN
         CALL Info('DetermineContact','Allocating projector mortar vectors',Level=10)
         ALLOCATE( MortarBC % Active( totsize ), MortarBC % Rhs( totsize) )
         MortarBC % Active = .FALSE.
         MortarBC % Rhs = 0.0_dp
         MortarBC % Perm => Perm 

         IF( AddDiag ) THEN
           ALLOCATE( MortarBC % Diag( totsize ) )
           MortarBC % Diag = 0.0_dp
         END IF

         RETURN
       END IF

       
       ! If permutation has changed we need to change the vectors also
       SamePerm = ANY( Perm /= MortarBC % Perm )
       SameSize = ( SIZE(MortarBC % Rhs) == totsize )

       ! Permutation unchanged, just return
       IF( SamePerm ) THEN
         DEALLOCATE( Perm ) 
         RETURN
       END IF
       
       ! Permutation changes, and also sizes changed
       IF(.NOT. SameSize ) THEN
         DEALLOCATE( MortarBC % Rhs )
         ALLOCATE( MortarBC % Rhs( totsize ) )
         MortarBC % Rhs = 0.0_dp
       END IF

       ! .NOT. SamePerm
       ALLOCATE(Active(totsize))
       Active = .FALSE.

       IF( AddDiag ) THEN
         ALLOCATE( Diag(totsize) )
         Diag = 0.0_dp
       END IF


       DO i=1,SIZE( Perm ) 
         j = Perm(i)
         IF( j == 0 ) CYCLE

         k = MortarBC % Perm(i)
         IF( k == 0 ) CYCLE

         DO l=1,Dofs
           Active(Dofs*(j-1)+l) = MortarBC % Active(Dofs*(k-1)+l)
         END DO
       END DO

       DEALLOCATE( MortarBC % Active ) 
       DEALLOCATE( MortarBC % Perm ) 
       MortarBC % Active => Active 
       MortarBC % Perm => Perm 

       IF( AddDiag ) THEN
         IF( ASSOCIATED( MortarBC % Diag ) ) THEN
           DEALLOCATE( MortarBC % Diag ) 
         END IF
         MortarBC % Diag => Diag 
       END IF

       CALL Info('DetermineContact','Copied > Active < flag to changed projector',Level=8)

     END SUBROUTINE InitializeMortarVectors

     

     ! Make a list of interface dofs to allow conservative algorithms. 
     ! There only nodes that are at the interface are added or removed from the set.
     !------------------------------------------------------------------------------
     SUBROUTINE MarkInterfaceDofs()
       
       INTEGER :: i,j,i2,j2,k,k2,l,n,ind,ind2,elem
       INTEGER, POINTER :: Indexes(:)
       TYPE(Element_t), POINTER :: Element
       
       CALL Info('DetermineContact','Marking interface dofs for conservative adding/removal',Level=8)

       IF(.NOT. ALLOCATED( InterfaceDof ) ) THEN
         ALLOCATE( InterfaceDof( SIZE(MortarBC % Active) ) )
       END IF
       InterfaceDof = .FALSE. 


       DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         
         Element => Mesh % Elements( elem )         
         IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE
         
         n = Element % TYPE % NumberOfNodes
         Indexes => Element % NodeIndexes
         
         DO i=1,n
           j = FieldPerm( Indexes(i) )
           IF( j == 0 ) CYCLE
           k = MortarBC % Perm( Indexes(i) )
           
           DO i2 = i+1,n
             j2 = FieldPerm( Indexes(i2) )
             IF( j2 == 0 ) CYCLE
             k2 = MortarBC % perm( Indexes(i2) )
             
             DO l=1,Dofs             
               ind = Dofs * ( k - 1 ) + l
               ind2 = Dofs * ( k2 - 1) + l
               
               IF( MortarBC % Active(ind) .NEQV. MortarBC % Active(ind2) ) THEN
                 InterfaceDof(ind) = .TRUE.
                 InterfaceDof(ind2) = .TRUE.
               END IF
             END DO
           END DO
         END DO
       END DO

       n = COUNT(InterfaceDof)
       CALL Info('DetermineContact',&
           'Number of interface dofs: '//TRIM(I2S(n)),Level=8)
     END SUBROUTINE MarkInterfaceDofs
     

     ! Calculates the signed distance that is used to define whether we have contact or not.
     ! If distance is negative then we can later add the corresponding node to the contact set
     ! Also computes the right-hand-side of the mortar equality constrained which is the 
     ! desired distance in the active direction. Works also for residual mode which greatly 
     ! improves the convergence for large displacements.  
     !----------------------------------------------------------------------------------------
     SUBROUTINE CalculateMortarDistance()

       REAL(KIND=dp) :: Disp(3), Coord(3), PrevDisp(3), Velo(3), ContactDist(3), ContactVelo(3), &
           LocalNormal0(3), SlipCoord(3)
       REAL(KIND=dp), POINTER :: DispVals(:), PrevDispVals(:) 
       REAL(KIND=dp) :: MinDist, MaxDist
       TYPE(Matrix_t), POINTER :: ActiveProjector
       LOGICAL :: IsSlave, IsMaster, DistanceSet
       LOGICAL, ALLOCATABLE :: SlaveNode(:), MasterNode(:), NodeDone(:)
       INTEGER, POINTER :: Indexes(:)
       INTEGER :: elemcode, CoeffSign
       REAL(KIND=dp), ALLOCATABLE :: CoeffTable(:)
       INTEGER :: l2
       LOGICAL :: DebugNormals

       CALL Info('DetermineContact','Computing distance between mortar boundaries',Level=14)

       DispVals => Solver % Variable % Values
       IF( .NOT. ASSOCIATED( DispVals ) ) THEN
         CALL Fatal('CalculateMortarDistance','Displacement variable not associated!')
       END IF

       IF( CalculateVelocity ) THEN
         IF( .NOT. ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
           CALL Fatal('CalculateMortarDistance','Displacement PrevValues not associated!')         
         END IF
         IF( Solver % TimeOrder == 1 ) THEN
           PrevDispVals => Solver % Variable % PrevValues(:,1)
         ELSE
           PrevDispVals => Solver % Variable % PrevValues(:,3)
         END IF
         IF(.NOT. ASSOCIATED( PrevDispVals ) ) CALL Fatal('DetermineContact',&
             'Previous displacement field required!')
       END IF


       ALLOCATE( SlaveNode( Mesh % NumberOfNodes ) ) 
       SlaveNode = .FALSE.

       IF( CreateDual ) THEN
         ALLOCATE( MasterNode( Mesh % NumberOfNodes ) ) 
         MasterNode = .FALSE.
       END IF

       DO i=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( i )                  
         IF( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) THEN
           SlaveNode( Element % NodeIndexes ) = .TRUE.
         END IF
         IF( CreateDual ) THEN
           IF ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) THEN
             MasterNode( Element % NodeIndexes ) = .TRUE.
           END IF
         END IF
       END DO

       ! First create the master, then the slave if needed
       IsSlave = .TRUE.
       IsMaster = .NOT. IsSlave
       ActiveProjector => Projector
       MinDist = HUGE(MinDist)
       MaxDist = -HUGE(MaxDist)

       IF( .NOT. ASSOCIATED( ActiveProjector ) ) THEN
         CALL Fatal('CalculateMortarDistance','Projector not associated!')
       END IF

       DebugNormals = ListGetLogical( Params,'Debug Normals',Found ) 

100    CONTINUE

       DO i = 1,ActiveProjector % NumberOfRows

         j = ActiveProjector % InvPerm(i)

         IF( j == 0 ) CYCLE

         wsum = 0.0_dp
         Dist = 0.0_dp
         DistN = 0.0_dp
         DistT1 = 0.0_dp
         DistT2 = 0.0_dp
         ContactVelo = 0.0_dp
         ContactDist = 0.0_dp
         DistanceSet = .FALSE.

         ! This is the most simple contact condition. We just want no slip on the contact.
         IF( TieContact .AND. .NOT. ResidualMode ) GOTO 200


         IF( ThisRotatedContact ) THEN
           Rotated = GetSolutionRotation(NTT, j )
           LocalNormal = NTT(:,1)
           LocalNormal0 = LocalNormal
           LocalT1 = NTT(:,2)
           IF( Dofs == 3 ) LocalT2 = NTT(:,3)
         ELSE
           LocalNormal = ContactNormal
           LocalT1 = ContactT1
           IF( Dofs == 3 ) LocalT2 = ContactT2 
         END IF

         ! Compute normal direction from the average sum of normals
         IF( NodalNormal ) THEN
           LocalNormal = 0.0_dp
           LocalT1 = 0.0_dp
           LocalT2 = 0.0_dp

           DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
             k = ActiveProjector % Cols(j)

             l = FieldPerm( k ) 
             IF( l == 0 ) CYCLE

             coeff = ActiveProjector % Values(j)             
             Rotated = GetSolutionRotation(NTT, k )

             ! Weighted direction for the unit vectors
             LocalNormal = LocalNormal + coeff * NTT(:,1)
             LocalT1 = LocalT1 + coeff * NTT(:,2)
             IF( Dofs == 3 ) LocalT2 = LocalT2 + coeff * NTT(:,3)
           END DO

           ! Normalize the unit vector length to one
           LocalNormal = LocalNormal / SQRT( SUM( LocalNormal**2 ) )
           LocalT1 = LocalT1 / SQRT( SUM( LocalT1**2 ) )
           IF( Dofs == 3 ) LocalT2 = LocalT2 / SQRT( SUM( LocalT1**2 ) )

           !PRINT *,'NodalNormal:',i,j,LocalNormal0,LocalNormal
         END IF

         ! For debugging reason, check that normals are roughly opposite
         IF( DebugNormals ) THEN
           DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
             k = ActiveProjector % Cols(j)
             
             l = FieldPerm( k ) 
             IF( l == 0 ) CYCLE

             Rotated = GetSolutionRotation(NTT, k )
             coeff = SUM( LocalNormal * NTT(:,1) ) + SUM( LocalT1*NTT(:,2)) + SUM(LocalT2*NTT(:,3))
             IF( SlaveNode(k) .AND. coeff < 2.5_dp ) THEN
               Found = .TRUE.
               !PRINT *,'Slave Normal:',i,j,k,Rotated,coeff
             ELSE IF( .NOT. SlaveNode(k) .AND. coeff > -2.5_dp ) THEN
               Found = .TRUE.
               !PRINT *,'Master Normal:',i,j,k,Rotated,coeff
             ELSE
               Found = .FALSE.
             END IF
             IF( Found ) THEN
               !PRINT *,'Prod:',SUM( LocalNormal * NTT(:,1) ), SUM( LocalT1*NTT(:,2)), SUM(LocalT2*NTT(:,3))
               !PRINT *,'N:',LocalNormal,NTT(:,1)
               !PRINT *,'T1:',LocalT1,NTT(:,2)
               !PRINT *,'T2:',LocalT2,NTT(:,3)
             END IF
           END DO
         END IF


         DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
           k = ActiveProjector % Cols(j)

           l = FieldPerm( k ) 
           IF( l == 0 ) CYCLE

           ! This includes only the coordinate since the displacement
           ! is added to the coordinate!
           coeff = ActiveProjector % Values(j)
           CoeffSign = 1

           ! Only compute the sum related to the active projector
           IF( SlaveNode(k) ) THEN
             wsum = wsum + coeff
           ELSE IF( ThisRotatedContact ) THEN
             CoeffSign = -1
           END IF
             
           IF( dofs == 2 ) THEN
             disp(1) = DispVals( 2 * l - 1)
             disp(2) = DispVals( 2 * l )
             disp(3) = 0.0_dp
           ELSE
             disp(1) = DispVals( 3 * l - 2)
             disp(2) = DispVals( 3 * l - 1 )
             disp(3) = DispVals( 3 * l )
           END IF

           ! If nonliear analysis is used we may need to cancel the introduced gap due to numerical errors 
           IF( TieContact .AND. ResidualMode ) THEN
             IF( ThisRotatedContact ) THEN
               ContactDist(1) = ContactDist(1) + coeff * SUM( LocalNormal * Disp )
               ContactDist(2) = ContactDist(2) + coeff * SUM( LocalT1 * Disp )
               IF( Dofs == 3) ContactDist(3) = ContactDist(3) + coeff * SUM( LocalT2 * Disp )
             ELSE
               ContactDist(1) = ContactDist(1) + coeff * SUM( ContactNormal * Disp )
               ContactDist(2) = ContactDist(2) + coeff * SUM( ContactT1 * Disp )
               IF( Dofs == 3 ) ContactDist(3) = ContactDist(3) + coeff * SUM( ContactT2 * Disp ) 
             END IF
             CYCLE
           END IF

           coord(1) = Mesh % Nodes % x( k ) 
           coord(2) = Mesh % Nodes % y( k ) 
           coord(3) = Mesh % Nodes % z( k ) 

           IF( CalculateVelocity ) THEN
             IF( dofs == 2 ) THEN
               PrevDisp(1) = PrevDispVals( 2 * l - 1)
               PrevDisp(2) = PrevDispVals( 2 * l )
               PrevDisp(3) = 0.0_dp
             ELSE
               PrevDisp(1) = PrevDispVals( 3 * l - 2)
               PrevDisp(2) = PrevDispVals( 3 * l - 1 )
               PrevDisp(3) = PrevDispVals( 3 * l )
             END IF
           END IF

           ! If the linear system is in residual mode also set the ContactDist in residual mode too!
           IF( ResidualMode ) Coord = Coord + Disp

           ! DistN is used to give the distance that we need to move the original coordinates
           ! in the wanted direction in order to have contact.
           IF( ThisRotatedContact ) THEN
             ContactDist(1) = ContactDist(1) + coeff * SUM( LocalNormal * Coord )
           ELSE
             ContactDist(1) = ContactDist(1) + coeff * SUM( ContactNormal * Coord )
           END IF

           ! Tangential distances needed to move the original coordinates to the contact position
           ! If stick is required then we want to keep the tangential slip zero. 
           IF( StickContact ) THEN             
             SlipCoord = -PrevDisp 
             IF( ResidualMode ) SlipCoord = SlipCoord + Disp 

             IF( ThisRotatedContact ) THEN
               ContactDist(2) = ContactDist(2) + coeff * SUM( LocalT1 * SlipCoord )
               IF( Dofs == 3) ContactDist(3) = ContactDist(3) + coeff * SUM( LocalT2 * SlipCoord )
             ELSE
               ContactDist(2) = ContactDist(2) + coeff * SUM( ContactT1 * SlipCoord )
               IF( Dofs == 3 ) ContactDist(3) = ContactDist(3) + coeff * SUM( ContactT2 * SlipCoord )
             END IF
           END IF

           ! If not in the residual mode still take into account the displacement for the condition
           IF( .NOT. ResidualMode ) Coord = Coord + Disp

           ! Dist is used to compute the current signed distance that is used to determine
           ! whether we have contact or not. 
           IF( RotationalProjector ) THEN
             Dist = Dist + coeff * SQRT( SUM( Coord**2 ) )
           ELSE
             Dist = Dist + coeff * SUM( ContactNormal * Coord )
           END IF

           IF( CalculateVelocity ) THEN
             Velo = ( Disp - PrevDisp ) !/ dt
             ContactVelo(1) = ContactVelo(1) + coeff * SUM( Velo * LocalNormal ) 
             ContactVelo(2) = ContactVelo(2) + coeff * SUM( Velo * LocalT1 )
             ContactVelo(3) = ContactVelo(3) + coeff * SUM( Velo * LocalT2 ) 
           END IF
           DistanceSet = .TRUE.
         END DO

         ! Divide by weight to get back to real distance in the direction of the normal
         IF( ABS( wsum ) > EPSILON( wsum )  ) THEN
           ContactDist = ContactDist / wsum 
           Dist = DistSign * Dist / wsum
           IF( CalculateVelocity ) THEN
             ContactVelo = ContactVelo / wsum
           END IF
         ELSE
           ContactDist = 0.0_dp
           Dist = 1.0_dp
           ContactVelo = 0.0_dp
         END IF

         ! PRINT *,'ContactVelo:',i, ContactVelo(1:Dofs), wsum, IsSlave

200      IF( IsSlave ) THEN
           MortarBC % Rhs(Dofs*(i-1)+DofN) = -ContactDist(1)
           IF( StickContact .OR. TieContact ) THEN
             MortarBC % Rhs(Dofs*(i-1)+DofT1) = -ContactDist(2) 
             IF( Dofs == 3 ) THEN
               MortarBC % Rhs(Dofs*(i-1)+DofT2) = -ContactDist(3)
             END IF
           END IF

           MinDist = MIN( Dist, MinDist ) 
           MaxDist = MAX( Dist, MaxDist )
         END IF

         IF( IsMaster ) THEN
           Dist = -Dist
           ContactVelo = -ContactVelo
         END IF

         ! We use the same permutation for all boundary variables
         IF(ActiveProjector % InvPerm(i) <= 0 ) CYCLE
         j = DistVar % Perm( ActiveProjector % InvPerm(i) )

         DistVar % Values( j ) = Dist

         GapVar % Values( j ) = ContactDist(1)

         IF( CalculateVelocity ) THEN
           DO k=1,Dofs             
             VeloVar % Values( Dofs*(j-1)+k ) = ContactVelo(k) 
           END DO
         END IF
       END DO

       IF( IsSlave ) THEN
         IF( CreateDual ) THEN
           IsSlave = .FALSE.
           IsMaster = .NOT. IsSlave
           ActiveProjector => DualProjector
           GOTO 100
         END IF
       END IF

       IF( CalculateVelocity ) THEN
         !PRINT *,'Velo range:',MINVAL( VeloVar % Values), MAXVAL( VeloVar % Values)
       END IF

       DEALLOCATE( SlaveNode )
       IF( CreateDual ) DEALLOCATE( MasterNode )

       !PRINT *,'Distance Range:',MinDist, MaxDist
       !PRINT *,'Distance Offset:',MINVAL( MortarBC % Rhs ), MAXVAL( MortarBC % Rhs )
       
     END SUBROUTINE CalculateMortarDistance



     ! Calculates the contact pressure in the normal direction from the nodal loads.
     ! The nodal loads may be given either in cartesian or n-t coordinate system. 
     !-------------------------------------------------------------------------------
     SUBROUTINE CalculateContactPressure()
       
       INTEGER :: elem
       INTEGER, POINTER :: Indexes(:)
       TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
       TYPE(Nodes_t) :: Nodes
       INTEGER :: i,j,k,t,CoordSys, NormalSign0, NormalSign, NormalCount
       REAL(KIND=dp) :: s, x, DetJ, u, v, w, Normal(3),NodalForce(3),DotProd, &
           NormalForce, SlipForce
       REAL(KIND=dp), ALLOCATABLE :: Basis(:)
       LOGICAL :: Stat, IsSlave, IsMaster
       TYPE(Matrix_t), POINTER :: ActiveProjector
       LOGICAL, ALLOCATABLE :: NodeDone(:)

       n = Mesh % MaxElementNodes
       ALLOCATE(Basis(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

       CoordSys = CurrentCoordinateSystem()
       NodalForce = 0.0_dp

       NormalSign0 = 0
       NormalCount = 0
       
       ALLOCATE( NodeDone( Mesh % NumberOfNodes ) )
       NodeDone = .FALSE.


100    DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         
         Element => Mesh % Elements( elem )         

         IsSlave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) 
         IsMaster = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) 

         IF( .NOT. ( IsSlave .OR. ( CreateDual .AND. IsMaster ) ) ) CYCLE
         
         Indexes => Element % NodeIndexes
         
         n = Element % TYPE % NumberOfNodes
         Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
         Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
         Nodes % z(1:n) = Mesh % Nodes % z(Indexes)
         
         IntegStuff = GaussPoints( Element, n )

         DO t=1,IntegStuff % n        
           U = IntegStuff % u(t)
           V = IntegStuff % v(t)
           W = IntegStuff % w(t)
           
           stat = ElementInfo( Element, Nodes, U, V, W, detJ, Basis )
           S = DetJ * IntegStuff % s(t)
           
           IF ( CoordSys /= Cartesian ) THEN
             X = SUM( Nodes % X(1:n) * Basis(1:n) )
             s = s * x
           END IF
           
           Normal = NormalVector( Element,Nodes,u,v,.TRUE. )

           ! Check the consistency of sign in the projector
           IF( IsSlave .AND. ( FlatProjector .OR. PlaneProjector ) ) THEN
             DotProd = SUM( Normal * ContactNormal ) 
             IF( DotProd < 0.0 ) THEN
               NormalSign = 1
             ELSE
               NormalSign = -1 
             END IF
             IF( NormalSign0 == 0 ) THEN
               NormalSign0 = NormalSign
             ELSE
               IF( NormalSign0 /= NormalSign ) NormalCount = NormalCount + 1
             END IF
           END IF

           DO i=1,n
             j = NormalLoadVar % Perm( Indexes(i) ) 

             IF( .NOT. NodeDone( Indexes(i) ) ) THEN             
               NodeDone( Indexes(i) ) = .TRUE.
               WeightVar % Values(j) = 0.0_dp
               NormalLoadVar % Values(j) = 0.0_dp
               SlipLoadVar % Values(j) = 0.0_dp
             END IF

             k = FieldPerm( Indexes(i) )
             IF( k == 0 .OR. j == 0 ) CYCLE
             DO l=1,dofs
               NodalForce(l) = LoadValues(dofs*(k-1)+l)
             END DO

             IF( ThisRotatedContact ) THEN
               NormalForce = NodalForce(1)
             ELSE
               NormalForce = SUM( NodalForce * Normal ) 
             END IF
             SlipForce = SQRT( SUM( NodalForce**2 ) - NormalForce**2 )

             NormalLoadVar % Values(j) = NormalLoadVar % Values(j) - &
                 s * Basis(i) * NormalForce
             SlipLoadVar % Values(j) = SlipLoadVar % Values(j) + &
                 s * Basis(i) * SlipForce
             
             WeightVar % Values(j) = WeightVar % Values(j) + s * Basis(i)
           END DO
           
         END DO
       END DO

       ! Normalize the computed normal loads such that the unit will be that of pressure
       DO i=1,Mesh % NumberOfNodes
         IF( NodeDone( i ) ) THEN             
           j = WeightVar % Perm(i)
           SlipLoadVar % Values(j) = SlipLoadVar % Values(j) / WeightVar % Values(j)**2
           NormalLoadVar % Values(j) = NormalLoadVar % Values(j) / WeightVar % Values(j)**2
         END IF
       END DO

       IF( FlatProjector .OR. PlaneProjector ) THEN
         IF( NormalCount == 0 ) THEN
           CALL Info('DetermineContact','All normals are consistently signed',Level=10)
         ELSE
           CALL Warn('DetermineContact','There are normals with conflicting signs: '&
               //TRIM(I2S(NormalCount) ) )
           NormalSign = 1
         END IF
         CALL Info('DetermineContact','Normal direction for distance measure: '&
             //TRIM(I2S(NormalSign)),Level=8)
         DistSign = NormalSign 
       END IF

       ! Check whether the normal sign has been enforced
       IF( ListGetLogical( BC,'Normal Sign Negative',Found ) ) DistSign = -1
       IF( ListGetLogical( BC,'Normal Sign Positive',Found ) ) DistSign = 1


       DEALLOCATE( Basis, Nodes % x, Nodes % y, Nodes % z, NodeDone )

     END SUBROUTINE CalculateContactPressure

  

     ! Sets the contact in the normal direction by looking at the signed distance and 
     ! contact force. The initial contact set may be enlarged to eliminate null-space
     ! related to rigid-body motion.
     !----------------------------------------------------------------------------------
     SUBROUTINE NormalContactSet() 
       
       INTEGER :: LimitSign, Removed, Added
       REAL(KIND=dp) :: DistOffSet, MinLoad, MaxLoad, NodeLoad, MinDist, MaxDist, NodeDist
       INTEGER :: i,j,k,ind
       LOGICAL :: Found

       ! This is related to the formulation of the PDE and is probably fixed for all elasticity solvers
       LimitSign = -1

       Removed = 0
       Added = 0        
       MinLoad = HUGE(MinLoad)
       MaxLoad = -HUGE(MaxLoad)
       MinDist = HUGE(MinDist)
       MaxDist = -HUGE(MaxDist)

       Found = .FALSE.
       IF( FirstTime ) THEN
         DistOffset = ListGetCReal( BC,&
             'Mortar BC Initial Contact Depth',Found)
         IF(.NOT. Found ) DistOffset = ListGetCReal( BC,&
             'Contact Depth Offset Initial',Found)
       END IF
       IF( .NOT. Found ) DistOffset = ListGetCReal( BC,&
           'Contact Depth Offset',Found)

       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = NormalLoadVar % Perm(j)

         ind = Dofs * (i-1) + DofN

         ! Tie contact should always be in contact - if we have found a counterpart
         IF( TieContact ) THEN
           MortarBC % Active(ind) = .TRUE.
           CYCLE
         END IF

         ! Enforce contact 
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Contact Active Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           MortarBC % Active(ind) = .TRUE.
           CYCLE
         END IF

         ! Enforce no contact
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Contact Passive Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           MortarBC % Active(ind) = .FALSE.
           CYCLE
         END IF

         ! Free nodes with wrong sign in contact force
         !--------------------------------------------------------------------------       
         IF( MortarBC % Active( ind ) ) THEN
           NodeLoad = NormalLoadVar % Values(k)
           MaxLoad = MAX( MaxLoad, NodeLoad )
           MinLoad = MIN( MinLoad, NodeLoad )
           DoRemove = ( LimitSign * NodeLoad > LimitSign * LoadEps ) 
           IF( DoRemove .AND. ConservativeRemove ) THEN
             DoRemove = InterfaceDof(ind) 
           END IF
           IF( DoRemove ) THEN
             removed = removed + 1
             MortarBC % Active(ind) = .FALSE.
           END IF
         ELSE 
           NodeDist = DistVar % Values(k)
           MaxDist = MAX( MaxDist, NodeDist ) 
           MinDist = MIN( MinDist, NodeDist )

           DoAdd = ( NodeDist < -ValEps + DistOffset )            
           IF( DoAdd .AND. ConservativeAdd ) THEN
             DoAdd = InterfaceDof(ind)
           END IF
           IF( DoAdd ) THEN
             added = added + 1
             MortarBC % Active(ind) = .TRUE.
           END IF
         END IF
       END DO

       IF( MaxDist - MinDist >= 0.0_dp ) THEN
         PRINT *,'NormalContactSet Dist:',MinDist,MaxDist
       END IF
       IF( MaxLoad - MinLoad >= 0.0_dp ) THEN
         PRINT *,'NormalContactSet Load:',MinLoad,MaxLoad
       END IF

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' nodes to the set'
         CALL Info('DetermineContactSet',Message,Level=5)
       END IF
       
       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' nodes from the set'
         CALL Info('DetermineContactSet',Message,Level=5)
       END IF

     END SUBROUTINE NormalContactSet



     ! If requested add new nodes to the contact set
     ! This would be typically done in order to make the elastic problem well defined
     ! Without any contact the bodies may float around.
     !---------------------------------------------------------------------------------
     SUBROUTINE IncreaseContactSet( LimitedMin )
       INTEGER :: LimitedMin

       REAL(KIND=dp), ALLOCATABLE :: DistArray(:)
       INTEGER, ALLOCATABLE :: IndArray(:)
       REAL(KIND=dp) :: Dist
       INTEGER :: i,j,ind,LimitedNow,NewNodes

       ! Nothing to do 
       IF( LimitedMin <= 0 ) RETURN

       LimitedNow = COUNT( MortarBC % active(DofN::Dofs) )      
       NewNodes = LimitedMin - LimitedNow
       IF( NewNodes <= 0 ) RETURN

       WRITE(Message,'(A,I0)') 'Initial number of contact nodes for '&
           //TRIM(VarName)//': ',LimitedNow 
       CALL Info('DetermineContactSet',Message,Level=5)

       CALL Info('DetermineContact',&
           'Setting '//TRIM(I2S(NewNodes))//' additional contact nodes',Level=5)

       ALLOCATE( DistArray( NewNodes ), IndArray( NewNodes ) ) 
       DistArray = HUGE( DistArray ) 
       IndArray = 0

       ! Find additional contact nodes from the closest non-contact nodes
       DO i = 1,Projector % NumberOfRows
         ind = Dofs * (i-1) + DofN
         IF( MortarBC % Active(ind)  ) CYCLE

         IF( Projector % InvPerm(i) == 0 ) CYCLE
         j = DistVar % Perm(Projector % InvPerm(i))
         Dist = DistVar % Values(j)
 
         IF( Dist < DistArray(NewNodes) ) THEN
           DistArray(NewNodes) = Dist
           IndArray(NewNodes) = i

           ! Order the new nodes such that the last node always has the largest distance
           ! This way we only need to compare to the one distance when adding new nodes.
           DO j=1,NewNodes-1
             IF( DistArray(j) > DistArray(NewNodes) ) THEN
               Dist = DistArray(NewNodes)
               DistArray(NewNodes) = DistArray(j)
               DistArray(j) = Dist                
               ind = IndArray(NewNodes)
               IndArray(NewNodes) = IndArray(j)
               IndArray(j) = ind                
             END IF
           END DO
         END IF
       END DO

       IF( ANY( IndArray == 0 ) ) THEN
         CALL Fatal('DetermineContact','Could not define sufficient number of new nodes!')
       END IF

       WRITE(Message,'(A,ES12.4)') 'Maximum distance needed for new nodes:',DistArray(NewNodes)
       CALL Info('DetermineContact',Message,Level=8)

       MortarBC % Active( Dofs*(IndArray-1)+DofN ) = .TRUE.

       DEALLOCATE( DistArray, IndArray ) 

     END SUBROUTINE IncreaseContactSet



     ! Sets the contact in the tangent direction(s) i.e. the stick condition.
     ! Stick condition in 1st and 2nd tangent condition are always the same. 
     !----------------------------------------------------------------------------------
     SUBROUTINE TangentContactSet() 
       
       INTEGER :: Removed0, Removed, Added
       REAL(KIND=dp) :: NodeLoad, TangentLoad, mustatic, mudynamic, stickcoeff, &
           Fstatic, Fdynamic, Ftangent, du(3), Slip
       INTEGER :: i,j,k,l,ind,IndN, IndT1, IndT2
       LOGICAL :: Found

       IF( FrictionContact .AND. &
           ListGetLogical( BC,'Stick Contact Global',Found ) ) THEN
        
         ! Sum up global normal and slide forces
         DO i = 1,Projector % NumberOfRows
           j = Projector % InvPerm( i ) 
           IF( j == 0 ) CYCLE
           k = FieldPerm( j ) 
           IF( k == 0 ) CYCLE
           k = NormalLoadVar % Perm(j)
                      
           ! If there is no contact there can be no stick either
           indN = Dofs * (i-1) + DofN
           IF( .NOT. MortarBC % Active(indN) ) CYCLE

           NodeLoad = NormalLoadVar % Values(k)
           TangentLoad = SlipLoadVar % Values(k)
         
           mustatic = ListGetRealAtNode( BC,'Static Friction Coefficient', j )
           mudynamic = ListGetRealAtNode( BC,'Dynamic Friction Coefficient', j )
           IF( mustatic <= mudynamic ) THEN
             CALL Warn('TangentContactSet','Static friction coefficient should be larger than dynamic!')
           END IF
           
           Fstatic = Fstatic + mustatic * ABS( NodeLoad ) 
           Fdynamic = Fdynamic + mudynamic * ABS( NodeLoad )
           Ftangent = Ftangent + ABS( TangentLoad ) 
           IF( Ftangent > Fstatic ) THEN
             SlipContact = .TRUE.
             FrictionContact = .FALSE.
           ELSE 
             GOTO 100
           END IF
         END DO
       END IF

       
       ! For stick and tie contact inherit the active flag from the normal component
       IF( SlipContact ) THEN
         MortarBC % Active( DofT1 :: Dofs ) = .FALSE.
         IF( Dofs == 3 ) THEN
            MortarBC % Active( DofT2 :: Dofs ) = .FALSE.
          END IF
          GOTO 100 
       ELSE IF( StickContact .OR. TieContact ) THEN
         MortarBC % Active( DofT1 :: Dofs ) = MortarBC % Active( DofN :: Dofs )
         IF( Dofs == 3 ) THEN
           MortarBC % Active( DofT2 :: Dofs ) = MortarBC % Active( DofN :: Dofs ) 
         END IF
         GOTO 100
       END IF

       CALL Info('TangentContactSet','Setting the stick set tangent components',Level=10)

       Removed0 = 0
       Removed = 0
       Added = 0        

       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = NormalLoadVar % Perm(j)

         indN = Dofs * (i-1) + DofN
         indT1 = ind - DofN + DofT1
         IF(Dofs == 3 ) indT2 = ind - DofN + DofT2

         ! If there is no contact there can be no stick either
         IF( .NOT. MortarBC % Active(indN) ) THEN
           IF( MortarBC % Active(indT1) ) THEN
             removed0 = removed0 + 1
             MortarBC % Active(indT1) = .FALSE.
             IF( Dofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           END IF
           CYCLE
         END IF

         ! Ok, we have normal contact what about stick
         ! Enforce stick condition
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Stick Active Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           IF( .NOT. MortarBC % Active(indT1) ) added = added + 1
           MortarBC % Active(indT1) = .TRUE.
           IF( Dofs == 3 ) MortarBC % Active(indT2) = .TRUE.
           CYCLE
         END IF

         ! Enforce no-stick condition (=slip)
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Stick Passive Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           IF( MortarBC % Active(IndT1) ) removed = removed + 1
           MortarBC % Active(indT1) = .FALSE.
           IF( Dofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           CYCLE
         END IF

         ! Remove nodes with too large tangent force
         !--------------------------------------------------------------------------       

         NodeLoad = NormalLoadVar % Values(k)
         TangentLoad = SlipLoadVar % Values(k)
         
         mustatic = ListGetRealAtNode( BC,'Static Friction Coefficient', j )
         mudynamic = ListGetRealAtNode( BC,'Dynamic Friction Coefficient', j )

         IF( mustatic <= mudynamic ) THEN
           CALL Warn('TangentContactSet','Static friction coefficient should be larger than dynamic!')
         END IF

         IF( MortarBC % Active(IndT1) ) THEN
           IF( TangentLoad > mustatic * ABS( NodeLoad ) ) THEN
             removed = removed + 1
             MortarBC % Active(indT1) = .FALSE.
             IF( Dofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           END IF
         ELSE              
           stickcoeff = ListGetRealAtNode( BC,'Stick Contact Coefficient', j, Found )
           IF( Found ) THEN
             DO l=1,Dofs             
               du(l) = VeloVar % Values( Dofs*(k-1)+l ) 
             END DO
             IF( Dofs == 3 ) THEN
               Slip = SQRT(du(dofT1)**2 + du(DofT2)**2)
             ELSE
               Slip = ABS( du(dofT1) )
             END IF
             IF( stickcoeff * slip  < mudynamic * ABS( NodeLoad ) ) THEN
               added = added + 1
               MortarBC % Active(indT1) = .TRUE.
               IF( Dofs == 3 ) MortarBC % Active(indT2) = .TRUE.
             END IF
           END IF
         END IF
       END DO

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' nodes to the stick set'
         CALL Info('DetermineContactSet',Message,Level=5)
       END IF
       
       IF(removed0 > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed0,' non-contact nodes from the stick set'
         CALL Info('DetermineContactSet',Message,Level=5)
       END IF

       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' sliding nodes from the stick set'
         CALL Info('DetermineContactSet',Message,Level=5)
       END IF


100    CALL Info('DetermineContactSet','Creating fields out of normal and stick contact sets',Level=10)

       DO i = 1, Projector % NumberOfRows
         j = Projector % InvPerm(i)
         IF( j == 0 ) CYCLE
         k = NormalActiveVar % Perm(j)

         IF( MortarBC % Active(Dofs*(i-1)+DofN) ) THEN
           NormalActiveVar % Values(k) = 1.0_dp
         ELSE
           NormalActiveVar % Values(k) = -1.0_dp
         END IF

         IF( MortarBC % Active(Dofs*(i-1)+DofT1) ) THEN
           StickActiveVar % Values(k) = 1.0_dp
         ELSE
           StickActiveVar % Values(k) = -1.0_dp
         END IF
       END DO

     END SUBROUTINE TangentContactSet



     ! Sets the diagonal entry for slip in the tangent direction(s).
     ! This coefficient may be used to relax the stick condition, and also to
     ! revert back nodes from slip to stick set. 
     !----------------------------------------------------------------------------------
     SUBROUTINE StickCoefficientSet() 
       
       REAL(KIND=dp) :: NodeLoad, TangentLoad
       INTEGER :: i,j,k,ind,IndN, IndT1, IndT2
       LOGICAL :: Found

       CALL Info('StickCoefficienttSet','Setting the stick coefficient entry for tangent components at stick',Level=10)

       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = NormalLoadVar % Perm(j)

         indN = Dofs * (i-1) + DofN
         indT1 = Dofs * (i-1) + DofT1
         IF(Dofs == 3 ) indT2 = Dofs * (i-1) + DofT2

         IF( .NOT. MortarBC % Active(indN) ) THEN
           ! If there is no contact there can be no stick either
           coeff = 0.0_dp            
         ELSE IF( .NOT. MortarBC % Active(indT1) ) THEN
           ! If there is no stick there can be no stick coefficient either
           coeff = 0.0_dp
         ELSE
           ! Get the stick coefficient
           coeff = ListGetRealAtNode( BC,'Stick Contact Coefficient', j )
         END IF

         MortarBC % Diag(indT1) = coeff
         IF( Dofs == 3 ) MortarBC % Diag(indT2) = coeff
       END DO
       
     END SUBROUTINE StickCoefficientSet



     ! Here we eliminate the middle nodes from the higher order elements if they 
     ! are different than both nodes of which they are associated with.
     ! There is no way geometric information could be accurate enough to allow
     ! such contacts to exist.
     !---------------------------------------------------------------------------
     SUBROUTINE QuadraticContactSet()

       LOGICAL :: ElemActive(8)
       INTEGER :: i,j,k,n,added, removed, elem, elemcode, ElemInds(8)
       INTEGER, POINTER :: Indexes(:)

       added = 0
       removed = 0

       DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( elem )         

         IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE

         Indexes => Element % NodeIndexes         
         n = Element % TYPE % NumberOfNodes
         elemcode = Element % Type % ElementCode

         DO i=1,n
           ElemActive(i) = MortarBC % Active( ElemInds(i) ) 
           IF(j>0) THEN
             ElemInds(i) = Dofs * ( j - 1) + DofN
             ElemActive(i) = MortarBC % Active( ElemInds(i) ) 
           ELSE
             ElemActive(i) = .FALSE.
           END IF
         END DO

         SELECT CASE ( elemcode ) 

         CASE( 202, 303, 404 ) 
           CONTINUE

         CASE( 203 )
           IF( ( ElemActive(1) .EQV. ElemActive(2) ) &
               .AND. ( ElemActive(1) .NEQV. ElemActive(3) ) ) THEN
             MortarBC % Active( ElemInds(3) ) = ElemActive(1) 
             IF( ElemActive(1) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

         CASE( 306 ) 
           IF( ( ElemActive(1) .EQV. ElemActive(2) ) &
               .AND. ( ElemActive(1) .NEQV. ElemActive(4) ) ) THEN
             MortarBC % Active( ElemInds(4) ) = ElemActive(1) 
             IF( ElemActive(1) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(2) .EQV. ElemActive(3) ) &
               .AND. ( ElemActive(2) .NEQV. ElemActive(5) ) ) THEN
             MortarBC % Active( ElemInds(5) ) = ElemActive(2) 
             IF( ElemActive(2) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(3) .EQV. ElemActive(1) ) &
               .AND. ( ElemActive(3) .NEQV. ElemActive(6) ) ) THEN
             MortarBC % Active( ElemInds(6) ) = ElemActive(3) 
             IF( ElemActive(3) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

         CASE( 408 ) 
           IF( ( ElemActive(1) .EQV. ElemActive(2) ) &
               .AND. ( ElemActive(1) .NEQV. ElemActive(5) ) ) THEN
             MortarBC % Active( ElemInds(5) ) = ElemActive(1) 
             IF( ElemActive(1) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(2) .EQV. ElemActive(3) ) &
               .AND. ( ElemActive(2) .NEQV. ElemActive(6) ) ) THEN
             MortarBC % Active( ElemInds(6) ) = ElemActive(2) 
             IF( ElemActive(2) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(3) .EQV. ElemActive(4) ) &
               .AND. ( ElemActive(3) .NEQV. ElemActive(7) ) ) THEN
             MortarBC % Active( ElemInds(7) ) = ElemActive(3) 
             IF( ElemActive(3) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(4) .EQV. ElemActive(1) ) &
               .AND. ( ElemActive(4) .NEQV. ElemActive(8) ) ) THEN
             MortarBC % Active( ElemInds(8) ) = ElemActive(4) 
             IF( ElemActive(4) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

         CASE DEFAULT
           CALL Fatal('NormalContactSet','Cannot deal with element: '//TRIM(I2S(elemcode)))

         END SELECT
       END DO

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' quadratic nodes to contact set'
         CALL Info('DetermineContactSet',Message,Level=5)
       END IF

       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' quadratic nodes from contact set'
         CALL Info('DetermineContactSet',Message,Level=5)
       END IF
         
     END SUBROUTINE QuadraticContactSet


     ! Project contact fields from slave to master
     !----------------------------------------------------------------------------------------
     SUBROUTINE ProjectFromSlaveToMaster()

       REAL(KIND=dp) :: Disp(3), Coord(3), PrevDisp(3), Velo(3), ContactDist(3), ContactVelo(3), &
           LocalNormal0(3), SlipCoord(3)
       REAL(KIND=dp), POINTER :: DispVals(:), PrevDispVals(:) 
       REAL(KIND=dp) :: MinDist, MaxDist, CoeffEps
       LOGICAL, ALLOCATABLE :: SlaveNode(:), NodeDone(:)
       REAL(KIND=dp), ALLOCATABLE :: CoeffTable(:), RealActive(:)
       INTEGER :: i,j,k,l,l2

       CALL Info('DetermineContact','Mapping entities from slave to master',Level=10)

       ALLOCATE( SlaveNode( Mesh % NumberOfNodes ) ) 
       SlaveNode = .FALSE.

       DO i=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( i )                  
         IF( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) THEN
           SlaveNode( Element % NodeIndexes ) = .TRUE.
         END IF
       END DO
 
       n = SIZE( DistVar % Values )
       ALLOCATE( CoeffTable( n ), NodeDone( n ) )
           

       CoeffTable = 0.0_dp
       NodeDone = .FALSE.
       

       DO i = 1,Projector % NumberOfRows             
         
         IF( Projector % InvPerm(i) == 0 ) CYCLE          
         l = DistVar % Perm( Projector % InvPerm(i) )
         
         DO j = Projector % Rows(i),Projector % Rows(i+1)-1
           k = Projector % Cols(j)
           
           IF( FieldPerm( k ) == 0 ) CYCLE               
           IF( SlaveNode( k ) ) CYCLE
           
           coeff = Projector % Values(j)
           
           l2 = DistVar % Perm( k )                
           
           IF(.NOT. NodeDone( l2 ) ) THEN
             DistVar % Values( l2 ) = 0.0_dp
             GapVar % Values( l2 ) = 0.0_dp
             NormalActiveVar % Values( l2 ) = 0.0_dp
             StickActiveVar % Values( l2 ) = 0.0_dp
             NormalLoadVar % Values( l2 ) = 0.0_dp
             SlipLoadVar % Values( l2 ) = 0.0_dp
             IF( CalculateVelocity ) THEN
               DO k=1,Dofs             
                 VeloVar % Values( Dofs*(l2-1)+k ) = 0.0_dp
               END DO
             END IF
             NodeDone( l2 ) = .TRUE.
           END IF

           CoeffTable( l2 ) = CoeffTable( l2 ) + coeff           
           DistVar % Values( l2 ) = DistVar % Values( l2 ) + coeff * DistVar % Values( l ) 
           GapVar % Values( l2 ) = GapVar % Values( l2 ) + coeff * GapVar % Values( l )
           NormalActiveVar % Values( l2 ) = NormalActiveVar % Values( l2 ) + coeff * NormalActiveVar % Values( l ) 
           StickActiveVar % Values( l2 ) = StickActiveVar % Values( l2 ) + coeff * StickActiveVar % Values( l ) 
           NormalLoadVar % Values( l2 ) = NormalLoadVar % Values( l2 ) + coeff * NormalLoadVar % Values( l ) 
           SlipLoadVar % Values( l2 ) = SlipLoadVar % Values( l2 ) + coeff * SlipLoadVar % Values( l ) 
           IF( CalculateVelocity ) THEN
             DO k=1,Dofs             
               VeloVar % Values( Dofs*(l2-1)+k ) = VeloVar % Values( Dofs*(l2-1)+k ) + &
                   coeff * VeloVar % Values( Dofs*(l-1)+k)
             END DO
           END IF
         END DO
       END DO
       
       CoeffEps = 1.0e-8 * MAXVAL( ABS( CoeffTable ) )
       DO i=1,SIZE( CoeffTable )            
         IF( NodeDone( i ) .AND. ( ABS( CoeffTable(i) ) > CoeffEps ) ) THEN
           DistVar % Values( i ) = DistVar % Values( i ) / CoeffTable( i ) 
           GapVar % Values( i ) = GapVar % Values( i ) / CoeffTable( i ) 
           NormalActiveVar % Values( i ) = NormalActiveVar % Values( i ) / CoeffTable( i ) 
           StickActiveVar % Values( i ) = StickActiveVar % Values( i ) / CoeffTable( i ) 

           IF( NormalActiveVar % Values( i ) >= 0.0_dp ) THEN
             NormalLoadVar % Values( i ) = NormalLoadVar % Values( i ) / CoeffTable( i ) 
             SlipLoadVar % Values( i ) = SlipLoadVar % Values( i ) / CoeffTable( i ) 
             IF( CalculateVelocity ) THEN
               DO k=1,Dofs
                 VeloVar % Values( Dofs*(i-1)+k ) = VeloVar % Values( Dofs*(i-1)+k ) / CoeffTable( i ) 
               END DO
             END IF
           ELSE
             NormalLoadVar % Values( i ) = 0.0_dp
             SlipLoadVar % Values( i ) = 0.0_dp
             IF( CalculateVelocity ) THEN
               DO k=1,Dofs
                 VeloVar % Values( Dofs*(i-1)+k ) = 0.0_dp
               END DO
             END IF             
           END IF

         END IF
       END DO

       DO i = 1, Projector % NumberOfRows
         j = Projector % InvPerm(i)
         IF( j == 0 ) CYCLE
         k = NormalActiveVar % Perm(j)
         
         IF( NormalActiveVar % Values( k ) < 0.0_dp ) THEN
           IF( CalculateVelocity ) THEN
             DO l=1,Dofs
               VeloVar % Values( Dofs*(k-1)+l ) = 0.0_dp
             END DO
           END IF
         END IF
       END DO


     END SUBROUTINE ProjectFromSlaveToMaster
   


     ! Set the friction in an implicit manner by copying matrix rows of the normal component
     ! to matrix rows of the tangential component multiplied by friction coefficient and 
     ! direction vector. 
     !---------------------------------------------------------------------------------------
     SUBROUTINE SetSlideFriction()

       REAL(KIND=dp), POINTER :: Values(:)
       LOGICAL, ALLOCATABLE :: NodeDone(:)
       REAL(KIND=dp) :: Coeff, ActiveLimit
       TYPE(Element_t), POINTER :: Element
       INTEGER, POINTER :: NodeIndexes(:)
       INTEGER :: i,j,k,k2,k3,l,l2,l3,n,t
       TYPE(Matrix_t), POINTER :: A
       LOGICAL :: Slave, Master, GivenDirection
       REAL(KIND=dp), POINTER :: VeloDir(:,:)
       REAL(KIND=dp) :: VeloCoeff(3),AbsVeloCoeff
       INTEGER :: VeloSign = 1


       IF(.NOT. ListCheckPresent( BC, 'Dynamic Friction Coefficient') ) RETURN
      
       CALL Info('DetermineContact','Setting contact friction for boundary',Level=10)

       GivenDirection = ListCheckPresent( BC,'Contact Velocity')
       IF(.NOT. GivenDirection ) THEN
         IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
           CALL Fatal('DetermineContact','Contact velocity must be given in some way')
         END IF
       END IF

       ActiveLimit = 0.0_dp
      
       Values => Solver % Matrix % values              
       ALLOCATE( NodeDone( SIZE( FieldPerm ) ) )
       A => Solver % Matrix        

       NodeDone = .FALSE.
       Coeff = 0.0_dp

       DO t = Mesh % NumberOfBulkElements+1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         Element => Mesh % Elements(t)
         
         Model % CurrentElement => Element
         
         Slave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag )
         Master = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag )
         
         IF( .NOT. ( Slave .OR. Master ) ) CYCLE

         NodeIndexes => Element % NodeIndexes
         n = Element % TYPE % NumberOfNodes
         
         DO i = 1, n
           j = Nodeindexes(i) 

           IF( NodeDone( j ) ) CYCLE
           IF( FieldPerm( j ) == 0 ) CYCLE

           ! Skipping the nodes not in the boundary
           k = NormalActiveVar % Perm( j )
           IF( k == 0 ) CYCLE
           
           ! Skipping the nodes not in contact. 
           IF( NormalActiveVar % Values( k ) <= -ActiveLimit ) CYCLE
           
           ! skipping the nodes in tangent stick
           IF( StickActiveVar % Values( k ) >= ActiveLimit ) CYCLE

           NodeDone( j ) = .TRUE.

           IF( Slave ) THEN
             Coeff = ListGetRealAtNode( BC,& 
                 'Dynamic Friction Coefficient', j, Found )
           ELSE
             Coeff = ListGetRealAtNode( MasterBC,& 
                 'Dynamic Friction Coefficient', j, Found )
             ! If friction not found in master then use the friction coefficient of the slave 
             ! Ideally they should be the same. 
             IF(.NOT. Found ) THEN
               Coeff = ListGetRealAtNode( BC,& 
                   'Dynamic Friction Coefficient', j, Found )               
             END IF
           END IF

           ! There is no point of setting too small friction coefficient
           IF(ABS(Coeff) < 1.0d-10) CYCLE

           IF( ThisRotatedContact ) THEN
             Rotated = GetSolutionRotation(NTT, j )
             LocalNormal = NTT(:,1)
             LocalT1 = NTT(:,2)
             IF( Dofs == 3 ) LocalT2 = NTT(:,3)
           ELSE
             Rotated = .FALSE.
             LocalNormal = ContactNormal
             LocalT1 = ContactT1
             IF( Dofs == 3 ) LocalT2 = ContactT2 
           END IF
           
           VeloCoeff = 0.0_dp           
           VeloSign = 1

           IF( GivenDirection ) THEN
             IF( Slave ) THEN
               VeloDir => ListGetConstRealArray( BC, &
                   'Contact Velocity', Found)
             ELSE
               VeloDir => ListGetConstRealArray( MasterBC, &
                   'Contact Velocity', Found)
               IF(.NOT. Found ) THEN
                 ! If velocity direction not found in master then use the opposite velocity of the slave
                 VeloDir => ListGetConstRealArray( BC, &
                     'Contact Velocity', Found)
                 VeloSign = -1
               END IF
             END IF
             VeloCoeff(DofT1) = SUM( VeloDir(1:3,1) * LocalT1 )
             IF( Dofs == 3 ) THEN
               VeloCoeff(DofT2) = SUM( VeloDir(1:3,1) * LocalT2 )
             END IF
           ELSE
             VeloCoeff(DofT1) = VeloVar % Values(Dofs*(k-1)+DofT1) 
             IF(Dofs==3) VeloCoeff(DofT2) = VeloVar % Values(Dofs*(k-1)+DofT2) 
             IF( .NOT. Slave .AND. .NOT. Rotated ) THEN
               VeloSign = -1
             END IF
           END IF

           ! Normalize coefficient to unity so that it only represents the direction of the force

           AbsVeloCoeff = SQRT( SUM( VeloCoeff**2 ) )
           IF( AbsVeloCoeff > TINY(AbsVeloCoeff) ) THEN
             VeloCoeff = VeloSign * VeloCoeff / AbsVeloCoeff
           ELSE
             CYCLE
           END IF

           ! Add the friction coefficient 
           VeloCoeff = Coeff * VeloCoeff 

           j = FieldPerm( j ) 
           k = DOFs * (j-1) + DofN 

           k2 = DOFs * (j-1) + DofT1 
           A % Rhs(k2) = A % Rhs(k2) - VeloCoeff(DofT1) * A % Rhs(k)

           IF( Dofs == 3 ) THEN
             k3 = DOFs * (j-1) + DofT2
             A % Rhs(k3) = A % Rhs(k3) - VeloCoeff(DofT2) * A % Rhs(k)             
           END IF

           DO l = A % Rows(k),A % Rows(k+1)-1
             DO l2 = A % Rows(k2), A % Rows(k2+1)-1
               IF( A % Cols(l2) == A % Cols(l) ) EXIT
             END DO

             A % Values(l2) = A % Values(l2) - VeloCoeff(DofT1) * A % Values(l)
             
             IF( Dofs == 3 ) THEN
               DO l3 = A % Rows(k3), A % Rows(k3+1)-1
                 IF( A % Cols(l3) == A % Cols(l) ) EXIT
               END DO
               A % Values(l3) = A % Values(l3) - VeloCoeff(DofT2) * A % Values(l)
             END IF
           END DO
         END DO
       END DO
       
       n = COUNT( NodeDone ) 
       CALL Info('SetSlideFriction','Number of friction nodes: '//TRIM(I2S(n)),Level=10)
       
       DEALLOCATE( NodeDone )

     END SUBROUTINE SetSlideFriction
     

   END SUBROUTINE DetermineContact




!------------------------------------------------------------------------------
!> Set dirichlet boundary condition for given dof. The conditions are
!> set based on the given name and applied directly to the matrix structure
!> so that a row is zeroed except for the diagonal which is set to one. 
!> Then the r.h.s. value determines the value of the field variable 
!> in the solution of the linear system.
!------------------------------------------------------------------------------
   SUBROUTINE SetDirichletBoundaries( Model, A, b, Name, DOF, NDOFs, Perm, &
       PermOffSet, OffDiagonalMatrix )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model          !< The current model structure
    TYPE(Matrix_t), POINTER :: A    !< The global matrix
    REAL(KIND=dp) :: b(:)           !< The global RHS vector
    CHARACTER(LEN=*) :: Name        !< Name of the dof to be set
    INTEGER :: DOF                  !< The order number of the dof
    INTEGER :: NDOFs                !< The total number of DOFs for this equation
    INTEGER :: Perm(:)              !< The node reordering info, this has been generated at the beginning of the 
                                    !< simulation for bandwidth optimization
    INTEGER, OPTIONAL :: PermOffSet  !< If the matrix and permutation vectors are not in sync the offset may used as a remedy. 
                                     !< Needed in fully coupled systems.
    LOGICAL, OPTIONAL :: OffDiagonalMatrix  !< For block systems the only the diagonal matrix should be given non-zero 
                                            !< entries for matrix and r.h.s., for off-diagonal matrices just set the row to zero.
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:), IndNodes(:), BCOrder(:)
    INTEGER, ALLOCATABLE :: Indexes(:), PassPerm(:)
    INTEGER :: BC,i,j,j2,k,l,m,n,t,k1,k2,OffSet
    LOGICAL :: GotIt, periodic, OrderByBCNumbering, ReorderBCs
    REAL(KIND=dp), POINTER :: MinDist(:)
    REAL(KIND=dp), POINTER :: WorkA(:,:,:) => NULL()
    REAL(KIND=dp) ::  s

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver

    LOGICAL :: Conditional
    LOGICAL, ALLOCATABLE :: DonePeriodic(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: CondName, DirName, PassName, PassCondName

    INTEGER :: NoNodes,NoDims,bf_id,nlen, NOFNodesFound, dim, &
        bndry_start, bndry_end, Upper
    REAL(KIND=dp), POINTER :: CoordNodes(:,:), Condition(:), Work(:),DiagScaling(:)
    REAL(KIND=dp) :: GlobalMinDist,Dist, Eps
    LOGICAL, ALLOCATABLE :: ActivePart(:), ActiveCond(:), ActivePartAll(:)
    TYPE(ValueList_t), POINTER :: ValueList, Params
    LOGICAL :: NodesFound, Passive, OffDiagonal, ApplyLimiter
    LOGICAL, POINTER :: LimitActive(:)
    TYPE(Variable_t), POINTER :: Var

    INTEGER :: ind, ElemFirst, ElemLast, bf, BCstr, BCend, BCinc
    REAL(KIND=dp) :: SingleVal
    LOGICAL :: AnySingleBC, AnySingleBF
    LOGICAL, ALLOCATABLE :: LumpedNodeSet(:)
    LOGICAL :: NeedListMatrix
    INTEGER, ALLOCATABLE :: Rows0(:), Cols0(:)
    REAL(KIND=dp), POINTER :: BulkValues0(:)

!------------------------------------------------------------------------------
! These logical vectors are used to minimize extra effort in setting up different BCs
!------------------------------------------------------------------------------

    DiagScaling => A % DiagScaling
    IF (.NOT.ASSOCIATED(DiagScaling)) THEN
      ALLOCATE(DiagScaling(A % NumberOFRows))
      DiagScaling = 1._dp
    END IF

    nlen = LEN_TRIM(Name)
    n = MAX( Model % NumberOfBodyForces,Model % NumberOfBCs)
    ALLOCATE( ActivePart(n), ActivePartAll(n), ActiveCond(n))
    CondName = Name(1:nlen) // ' Condition'
    PassName = Name(1:nlen) // ' Passive'
    PassCondName = Name(1:nlen) // ' Condition' // ' Passive'

    OffSet = 0
    OffDiagonal = .FALSE.
    IF( PRESENT( PermOffSet) ) OffSet = PermOffSet
    IF( PRESENT( OffDiagonalMatrix ) ) OffDiagonal = OffDiagonalMatrix

    ALLOCATE( Indexes(Model % Mesh % MaxElementDOFs) )
!------------------------------------------------------------------------------
! Go through the periodic BCs and set the linear dependence
!------------------------------------------------------------------------------

   ActivePart = .FALSE.
   DO BC=1,Model % NumberOfBCs
     IF ( ListGetLogical( Model % BCs(BC) % Values, &
         'Periodic BC ' // Name(1:nlen), GotIt ) ) ActivePart(BC) = .TRUE.
     IF ( ListGetLogical( Model % BCs(BC) % Values, &
         'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) ActivePart(BC) = .TRUE.
     IF ( ListCheckPresent( Model % BCs(BC) % Values, &
         'Periodic BC Scale ' // Name(1:nlen) ) ) ActivePart(BC) = .TRUE.
   END DO

   IF( ANY(ActivePart) ) THEN    
     IF( Offset > 0 ) THEN
       CALL Fatal('SetDirichletBoundaries','Periodicity not considered with offset')
     END IF

     ALLOCATE( DonePeriodic( Model % Mesh % NumberOFNodes ) )
     DonePeriodic = .FALSE.
     DO BC=1,Model % NumberOfBCs
       IF( ActivePart(BC) ) THEN
         CALL SetPeriodicBoundariesPass1( Model, A, b, Name, DOF, &
             NDOFs, Perm, BC, DonePeriodic )
       END IF
     END DO
     
     DonePeriodic = .FALSE.
     DO BC=1,Model % NumberOfBCs
       IF(ActivePart(BC)) THEN       
         CALL SetPeriodicBoundariesPass2( Model, A, b, Name, DOF, &
             NDOFs, Perm, BC, DonePeriodic )
       END IF
     END DO
     DEALLOCATE( DonePeriodic ) 

   END IF
   

! Add the possible friction coefficient
!----------------------------------------------------------
   IF ( ListCheckPresentAnyBC( Model,'Friction BC ' // Name(1:nlen) ) ) THEN
     CALL SetFrictionBoundaries( Model, A, b, Name, NDOFs, Perm )
   END IF


! Add the possible nodal jump in case of mortar projectors
!---------------------------------------------------------------
   IF( ListGetLogical( Model % Solver % Values,'Apply Mortar BCs',GotIt ) ) THEN
     CALL SetWeightedProjectorJump( Model, A, b, &
                      Name, DOF, NDOFs, Perm )
   END IF


!------------------------------------------------------------------------------
! Go through the normal Dirichlet BCs applied on the boundaries
!------------------------------------------------------------------------------

    ActivePart = .FALSE.
    ActiveCond = .FALSE.
    ActivePartAll = .FALSE.
    DO BC=1,Model % NumberOfBCs
      ActivePartAll(BC) = ListCheckPresent( &
            Model % BCs(bc) % Values, Name(1:nlen) // ' DOFs' )
      ActivePart(BC) = ListCheckPresent( Model % BCs(bc) % Values, Name ) 
      ActiveCond(BC) = ListCheckPresent( Model % BCs(bc) % Values, CondName )      
    END DO

    OrderByBCNumbering = ListGetLogical( Model % Simulation, &
       'Set Dirichlet BCs by BC Numbering', gotIt)

    BCOrder => ListGetIntegerArray( Model % Solver % Values, &
         'Dirichlet BC Order', ReorderBCs)
    IF(ReorderBCs) THEN
       IF(.NOT. OrderByBCNumbering) THEN
          CALL Warn('SetDirichletBoundaries',"Requested 'Dirichlet BC Order' but &
               &not 'Set Dirichlet BCs by BC Numbering', ignoring...")
       ELSE IF(SIZE(BCOrder) /= Model % NumberOfBCs) THEN
          CALL Fatal('SetDirichletBoundaries',"'Dirichlet BC Order' is the wrong length!")
       END IF
    END IF

    bndry_start = Model % NumberOfBulkElements+1
    bndry_end   = bndry_start+Model % NumberOfBoundaryElements-1

    ! check and set some flags for nodes belonging to n-t boundaries
    ! getting set by other bcs:
    ! --------------------------------------------------------------
    IF ( NormalTangentialNOFNodes>0 ) THEN
      IF ( OrderByBCNumbering ) THEN
        DO i=1,Model % NumberOfBCs
          BC = i
          IF(ReorderBCs) BC = BCOrder(BC)
          IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
          Conditional = ActiveCond(BC)

          DO t = bndry_start, bndry_end
            Element => Model % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                   Model % BCs(BC) % Tag ) CYCLE

            ValueList => Model % BCs(BC) % Values
            Model % CurrentElement => Element

            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              Indexes(1:n) = Element % NodeIndexes
            ELSE
              n = SgetElementDOFs( Indexes )
            END IF
            CALL CheckNTelement(n,t)
          END DO
        END DO
      ELSE
        DO t = bndry_start, bndry_end
          DO BC=1,Model % NumberOfBCs
            IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
            Conditional = ActiveCond(BC)
          
            Element => Model % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                 Model % BCs(BC) % Tag ) CYCLE
          
            ValueList => Model % BCs(BC) % Values
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              Indexes(1:n) = Element % NodeIndexes
            ELSE
              n = SgetElementDOFs( Indexes )
            END IF
            CALL CheckNTelement(n,t)
          END DO
        END DO
      END IF

      IF ( DOF<= 0 ) THEN
        DO t=bndry_start,bndry_end
          Element => Model % Elements(t)
          n = Element % TYPE % NumberOfNodes
          DO j=1,n
            k = BoundaryReorder(Element % NodeIndexes(j))
            IF (k>0) THEN
              NTelement(k,:)=0
              NTzeroing_done(k,:) = .FALSE.
            END IF
          END DO
        END DO
      END IF
    END IF

    ! Set the Dirichlet BCs from active boundary elements, if any...:
    !----------------------------------------------------------------
    IF( ANY(ActivePart) .OR. ANY(ActivePartAll) ) THEN    
      IF ( OrderByBCNumbering ) THEN
        DO i=1,Model % NUmberOfBCs
          BC = i
          IF(ReorderBCs) BC = BCOrder(BC)
          IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
          Conditional = ActiveCond(BC)

          DO t = bndry_start, bndry_end
            Element => Model % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                Model % BCs(BC) % Tag ) CYCLE
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              Indexes(1:n) = Element % NodeIndexes
            ELSE
              n = SgetElementDOFs( Indexes )
            END IF
            ValueList => Model % BCs(BC) % Values
            CALL SetElementValues(n,t)
          END DO
        END DO
      ELSE
        DO t = bndry_start, bndry_end
          DO BC=1,Model % NumberOfBCs
            IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
            Conditional = ActiveCond(BC)
            
            Element => Model % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                Model % BCs(BC) % Tag ) CYCLE
            
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              Indexes(1:n) = Element % NodeIndexes
            ELSE
              n = SgetElementDOFs( Indexes )
            END IF
            ValueList => Model % BCs(BC) % Values
            CALL SetElementValues(n,t)
          END DO
        END DO
      END IF
    END IF


!------------------------------------------------------------------------------
! Go through the Dirichlet conditions in the body force lists
!------------------------------------------------------------------------------
    
    ActivePart = .FALSE.
    ActiveCond = .FALSE.
    ActivePartAll = .FALSE.
    Passive = .FALSE.
    DO bf_id=1,Model % NumberOFBodyForces
      ValueList => Model % BodyForces(bf_id) % Values

      ActivePartAll(bf_id) = ListCheckPresent(ValueList, Name(1:nlen) // ' DOFs' ) 
      ActiveCond(bf_id) = ListCheckPresent( ValueList,CondName )      
      ActivePart(bf_id) = ListCheckPresent(ValueList, Name(1:nlen) ) 

      Passive = Passive .OR. ListCheckPresent(ValueList, PassName)
    END DO

    IF ( ANY(ActivePart) .OR. ANY(ActivePartAll) ) THEN

      Solver => Model % Solver
      Mesh   => Solver % Mesh
      ALLOCATE(PassPerm(Mesh % NumberOfNodes),NodeIndexes(1));PassPerm=0
      DO i=0,Mesh % PassBCCnt-1
        j=Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements-i
        PassPerm(Mesh % Elements(j) % NodeIndexes)=1
      END DO

      DO t = 1, Mesh % NumberOfBulkElements 
        Element => Mesh % Elements(t)
        bf_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force', GotIt)
        
        IF(.NOT. GotIt) CYCLE
        IF(.NOT. ActivePart(bf_id) .AND. .NOT. ActivePartAll(bf_id)) CYCLE
        Conditional = ActiveCond(bf_id)

        Model % CurrentElement => Element
        
        IF ( ActivePart(bf_id) ) THEN
          n = Element % TYPE % NumberOfNodes
          Indexes(1:n) = Element % NodeIndexes
        ELSE
          n = SgetElementDOFs( Indexes )
        END IF
        ValueList => Model % BodyForces(bf_id) % Values

        IF (ListGetLogical(ValueList,PassCondName,GotIt)) THEN
          IF (.NOT.CheckPassiveElement(Element)) CYCLE
          DO j=1,n
            NodeIndexes(1) = Indexes(j)
            IF(PassPerm(NodeIndexes(1))==0) CALL SetPointValues(1)
          END DO
        ELSE
          CALL SetElementValues( n,t )
        END IF
      END DO
      DEALLOCATE(NodeIndexes,PassPerm)
    END IF
    DEALLOCATE(ActivePart, ActiveCond)


!------------------------------------------------------------------------------
! Go through the pointwise Dirichlet BCs that are created on-the-fly
! Note that it is best that the coordinates are transformed to nodes using 
! the right variable. Otherwise it could point to nodes that are not active.
!------------------------------------------------------------------------------
     
    DO BC=1,Model % NumberOfBCs
      
      ValueList => Model % BCs(BC) % Values
      IF( .NOT. ListCheckPresent( ValueList,Name )) CYCLE
      NodesFound = ListCheckPresent( ValueList,'Target Nodes' )
      
      ! The coordinates are only requested for a body that has no list of nodes. 
      ! At the first calling the list of coorinates is transformed to list of nodes. 
      IF(.NOT. NodesFound) THEN
        CoordNodes => ListGetConstRealArray(ValueList,'Target Coordinates',GotIt)
        IF(GotIt) THEN
          Eps = ListGetConstReal( ValueList, 'Target Coordinates Eps', Gotit )
          IF ( .NOT. GotIt ) THEN
            Eps = HUGE(Eps)
          ELSE
            ! We are looking at square of distance
            Eps = Eps**2
          END IF

          NoNodes = SIZE(CoordNodes,1)
          NoDims = SIZE(CoordNodes,2)
          
          IF(NoNodes > 0) THEN               
            ALLOCATE( IndNodes(NoNodes), MinDist(NoNodes) )
            IndNodes = -1
            MinDist = HUGE( Dist )
            DO j=1,NoNodes
              DO i=1,Model % NumberOfNodes
                IF( Perm(i) == 0) CYCLE
                
                Dist = (Model % Mesh % Nodes % x(i) - CoordNodes(j,1))**2 
                IF(NoDims >= 2) Dist = Dist + (Model % Mesh % Nodes % y(i) - CoordNodes(j,2))**2
                IF(NoDims == 3) Dist = Dist + (Model % Mesh % Nodes % z(i) - CoordNodes(j,3))**2
                Dist = SQRT(Dist)
                
                IF(Dist < MinDist(j) .AND. Dist <= Eps ) THEN
                  MinDist(j) = Dist
                  IndNodes(j) = i
                END IF
              END DO
            END DO

            ! In parallel case eliminate all except the nearest node. 
            ! This relies on the fact that for each node partition the 
            ! distance to nearest node is computed accurately. 
            DO j=1,NoNodes
              GlobalMinDist = ParallelReduction( MinDist(j), 1 )
              IF( ABS( GlobalMinDist - MinDist(j) ) > TINY(Dist) ) THEN
                IndNodes(j) = 0
              END IF
            END DO

            NOFNodesFound = 0
            DO j=1,NoNodes
               IF ( IndNodes(j)>0 ) THEN
                 NOFNodesFound = NOFNodesFound+1
                 IndNodes(NOFNodesFound) = IndNodes(j)
               END IF
            END DO
            
            ! In the first time add the found nodes to the list structure
            IF ( NOFNodesFound > 0 ) THEN
              CALL ListAddIntegerArray( ValueList,'Target Nodes', &
                  NOFNodesFound, IndNodes) 
              NodesFound = .TRUE.            
            ELSE
              ! If no nodes found, add still an empty list and make sure the 
              ! zero is not treated later on. Otherwise this search would be 
              ! retreated each time. 
              CALL ListAddIntegerArray( ValueList,'Target Nodes', &
                  1, IndNodes) 
            END IF

            ! Finally deallocate the temporal vectors
            DEALLOCATE( IndNodes, MinDist ) 
          END IF
        END IF
      END IF
      
      ! If the target coordinates has already been assigned to an empty list 
      ! cycle over it by testing the 1st node. 
      IF( NodesFound ) THEN
        NodeIndexes => ListGetIntegerArray( ValueList,'Target Nodes')
        IF( NodeIndexes(1) == 0 ) NodesFound = .FALSE. 
      END IF

      IF(NodesFound) THEN           
        Conditional = ListCheckPresent( ValueList, CondName )      
        n = SIZE(NodeIndexes)
        CALL SetPointValues(n)
      END IF
    END DO


!------------------------------------------------------------------------------
!   Go through soft upper and lower limits
!------------------------------------------------------------------------------
    Params => Model % Solver % Values
    ApplyLimiter = ListGetLogical( Params,'Apply Limiter',GotIt) 

    IF( Dof/=0 .AND. ApplyLimiter ) THEN
      CALL Info('SetDirichletBoundaries','Applying limiters',Level=10)

      DO Upper=0,1
        
        ! The limiters have been implemented only componentwise
        !-------------------------------------------------------
        
        NULLIFY( LimitActive ) 
        Var => Model % Solver % Variable
        IF( Upper == 0 ) THEN
          IF( ASSOCIATED( Var % LowerLimitActive ) ) &
              LimitActive => Var % LowerLimitActive
        ELSE
          IF( ASSOCIATED( Var % UpperLimitActive ) ) &
              LimitActive => Var % UpperLimitActive
        END IF
        
        IF( .NOT. ASSOCIATED( LimitActive ) ) CYCLE
        
        IF( Upper == 0 ) THEN
          CondName = TRIM(name)//' Lower Limit' 
        ELSE
          CondName = TRIM(name)//' Upper Limit' 
        END IF
        
        ! check and set some flags for nodes belonging to n-t boundaries
        ! getting set by other bcs:
        ! --------------------------------------------------------------
        DO t = 1, Model % NumberOfBulkElements+Model % NumberOfBoundaryElements
          Element => Model % Elements(t)
          Model % CurrentElement => Element
          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          
          IF( t > Model % NumberOfBulkElements ) THEN
            DO bc = 1,Model % NumberOfBCs
              IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
              ValueList => Model % BCs(BC) % Values
              CALL SetLimiterValues(n)
            END DO
          ELSE             
            bf_id = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                'Body Force', GotIt)
            IF(.NOT. GotIt ) CYCLE
            ValueList => Model % Bodyforces(bf_id) % Values
            CALL SetLimiterValues(n)
          END IF
        END DO
      END DO
    END IF
    

    ! Check the boundaries and body forces for possible single nodes BCs that are used to fixed 
    ! the domain for undetermined equations. The loop is slower than optimal in the case that there is 
    ! a large amount of different boundaries that have a node to set. 
    !--------------------------------------------------------------------------------------------
    DirName = TRIM(Name)//' Single Node'
    AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
    AnySingleBF = ListCheckPresentAnyBodyForce( Model, DirName )

    IF( AnySingleBC .OR. AnySingleBF ) THEN
      DO bc = 1,Model % NumberOfBCs  + Model % NumberOfBodyForces    

        ! Make a distinction between BCs and BFs. 
        ! These are treated in the same loop because most of the logic is still the same. 
        IF( bc <= Model % NumberOfBCs ) THEN
          IF(.NOT. AnySingleBC ) CYCLE
          ValueList => Model % BCs(BC) % Values
          ElemFirst =  Model % NumberOfBulkElements + 1 
          ElemLast =  Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
        ELSE
          IF( .NOT. AnySingleBF ) CYCLE
          ValueList => Model % BodyForces(bc - Model % NumberOfBCs) % Values
          ElemFirst =  1
          ElemLast =  Model % NumberOfBulkElements
        END IF

        SingleVal = ListGetCReal( ValueList,DirName, GotIt) 
        IF( .NOT. GotIt ) CYCLE
        ind = ListGetInteger( ValueList,TRIM(Name)//' Single Node Index',GotIt )     
        
        ! On the first time find a one single uniquely defined node for setting 
        ! the value. In parallel it will be an unshared node with the highest possible 
        ! node number 
        IF( .NOT. GotIt ) THEN        

          ind = 0
          DO t = ElemFirst, ElemLast
            Element => Model % Elements(t)
            n = Element % TYPE % NumberOfNodes
            NodeIndexes => Element % NodeIndexes

            IF( bc <= Model % NumberOfBCs ) THEN
              IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
            ELSE
              bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values,'Body Force',GotIt)
              IF( bc - Model % NumberOfBCs /= bf ) CYCLE
            END IF

            DO i=1,n
              j = NodeIndexes(i)
              IF( Perm(j) == 0) CYCLE
              IF( ParEnv % PEs > 1 ) THEN
                IF( SIZE( Mesh % ParallelInfo % NeighbourList(j) % Neighbours) > 1 ) CYCLE               
                IF( Mesh % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPe ) CYCLE               
              END IF
              ind = j 
              EXIT
            END DO
            IF( ind > 0 ) EXIT
          END DO

          ! Find the maximum partition that owns the node. 
          ! It could be minimum also, just some convection is needed. 
          IF( ParEnv % PEs > 1 ) THEN
            k = -1
            IF( ind > 0 ) k = ParEnv % MyPe          
            k = NINT( ParallelReduction( 1.0_dp * k, 2 ) ) 
            IF( k == -1 ) THEN
              CALL Warn('SetDirichletBoundaries','Could not find node to set: '//TRIM(DirName))
            END IF
            IF( k /= ParEnv % MyPe ) ind = 0
          ELSE
            IF( ind == 0 ) THEN
              CALL Warn('SetDirichletBoundaries','Could not find node to set: '//TRIM(DirName))
            END IF
          END IF

          CALL ListAddInteger( ValueList,TRIM(Name)//' Single Node Index', ind )
        END IF

        ! Ok, if this is the partition where the single node to eliminate the floating should 
        ! be eliminated then set it here.         
        IF( ind > 0 ) THEN
          CALL SetSinglePoint(ind,DOF,SingleVal,.TRUE.)
        END IF
      END DO
    END IF


!------------------------------------------------------------------------------
!   Take care of the matrix entries of passive elements
!------------------------------------------------------------------------------

    IF ( Passive ) THEN
      Solver => Model % Solver
      Mesh => Solver % Mesh
      DO i=1,Solver % NumberOfActiveElements
        Element => Mesh % Elements(Solver % ActiveElements(i))
        IF (CheckPassiveElement(Element)) THEN
          n = sGetElementDOFs(Indexes,UElement=Element)
          DO j=1,n
            k=Indexes(j)
            IF (k<=0) CYCLE

            k=Perm(k)
            IF (k<=0) CYCLE

            s=0._dp
            DO l=1,NDOFs
              m=NDOFs*(k-1)+l
              s=s+ABS(A % Values(A % Diag(m)))
            END DO
            IF (s>EPSILON(s)) CYCLE
 
            DO l=1,NDOFs
              m = NDOFs*(k-1)+l
              CALL SetSinglePoint(k,l,Solver % Variable % Values(m),.FALSE.)
            END DO
          END DO
        END IF
      END DO
    END IF


    ! Check the boundaries and body forces for possible single nodes BCs that must have a constant
    ! value on that boundary / body force.
    !--------------------------------------------------------------------------------------------
    DirName = TRIM(Name)//' Constant'
    AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
    AnySingleBF = ListCheckPresentAnyBodyForce( Model, DirName )

    IF( AnySingleBC .OR. AnySingleBF ) THEN
      ALLOCATE( LumpedNodeSet( SIZE( Perm ) ) )

      IF( AnySingleBC ) CALL Info('SetDirichletBoundaries','Found BC constraint: '//TRIM(DirName))
      IF( AnySingleBF ) CALL Info('SetDirichletBoundaries','Found BodyForce constraint: '//TRIM(DirName))

      ! Improve the logic in future
      ! Now we assume that if the "supernode" has been found then also the matrix has the correct topology. 
      IF( AnySingleBC ) THEN
        NeedListMatrix = .NOT. ListCheckPresentAnyBC( Model, TRIM(Name)//' Constant Node Index')
      ELSE 
        NeedListMatrix = .NOT. ListCheckPresentAnyBodyForce( Model, TRIM(Name)//' Constant Node Index')
      END IF
      
      ! Move the list matrix because of its flexibility
      IF( NeedListMatrix ) THEN
        CALL Info('SetDirichletBoundaries','Using List maxtrix to set constant constraints',Level=8)
        CALL Info('SetDircihletBoundaries','Original matrix non-zeros: '&
            //TRIM(I2S(SIZE( A % Cols ))),Level=8)
        IF( ASSOCIATED( A % BulkValues ) ) THEN
          ALLOCATE( Cols0( SIZE( A % Cols ) ), Rows0( SIZE( A % Rows ) ) )
          Cols0 = A % Cols
          Rows0 = A % Rows
        END IF
        CALL List_toListMatrix(A)
      END IF

      DO bc = 1,Model % NumberOfBCs + Model % NumberOfBodyForces

        ! Make a distinction between BCs and BFs. 
        ! These are treated in the same loop because most of the logic is still the same. 
        IF( bc <= Model % NumberOfBCs ) THEN
          IF(.NOT. AnySingleBC ) CYCLE
          ValueList => Model % BCs(BC) % Values
          ElemFirst =  Model % NumberOfBulkElements + 1 
          ElemLast =  Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
        ELSE
          IF(.NOT. AnySingleBF ) CYCLE
          ValueList => Model % BodyForces(bc - Model % NumberOfBCs) % Values
          ElemFirst =  1
          ElemLast =  Model % NumberOfBulkElements
        END IF

        IF( .NOT. ListGetLogical( ValueList,DirName, GotIt) ) CYCLE

        ind = ListGetInteger( ValueList,TRIM(Name)//' Constant Node Index',GotIt )     
        

        ! On the first time find a one single uniquely defined node for setting 
        ! the value. In parallel it will be an unshared node with the highest possible 
        ! node number 
        IF( .NOT. GotIt ) THEN        

          ind = 0
          DO t = ElemFirst, ElemLast
            Element => Model % Elements(t)

            IF( bc <= Model % NumberOfBCs ) THEN
              IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
            ELSE
              bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values,'Body Force',GotIt)
              IF( bc - Model % NumberOfBCs /= bf ) CYCLE
            END IF            

            n = Element % TYPE % NumberOfNodes
            NodeIndexes => Element % NodeIndexes

            DO i=1,n
              j = NodeIndexes(i)
              IF( Perm(j) == 0) CYCLE
              IF( ParEnv % PEs > 1 ) THEN
                IF( SIZE( Mesh % ParallelInfo % NeighbourList(j) % Neighbours) > 1 ) CYCLE               
                IF( Mesh % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPe ) CYCLE               
               END IF
              ind = j 
              EXIT
            END DO
            IF( ind > 0 ) EXIT
          END DO

          ! Find the maximum partition that owns the node. 
          ! It could be minimum also, just some convection is needed. 
          IF( ParEnv % PEs > 1 ) THEN
            k = -1
            IF( ind > 0 ) k = ParEnv % MyPe          
            k = NINT( ParallelReduction( 1.0_dp * k, 2 ) ) 
            IF( k == -1 ) THEN
              CALL Warn('SetDirichletBoundaries','Could not find node to set: '//TRIM(DirName))
            END IF
            IF( k /= ParEnv % MyPe ) ind = 0
          ELSE
            IF( ind == 0 ) THEN
              CALL Warn('SetDirichletBoundaries','Could not find node to set: '//TRIM(DirName))
            END IF
          END IF

          CALL ListAddInteger( ValueList,TRIM(Name)//' Constant Node Index', ind )
          NeedListMatrix = .TRUE.
        END IF

        IF( ParEnv % PEs > 1 ) CALL Warn('SetDirichletBoundaries','Node index not set properly in parallel')
        IF( ind == 0 ) CYCLE

        ! Ok, now sum up the rows to the corresponding nodal index
        LumpedNodeSet = .FALSE.

        ! Don't lump the "supernode" and therefore mark it set already
        LumpedNodeSet(ind) = .TRUE.

        DO t = ElemFirst, ElemLast
          Element => Model % Elements(t)

          IF( bc <= Model % NumberOfBCs ) THEN
            IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
          ELSE
            bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values,'Body Force',GotIt)
            IF( bc - Model % NumberOfBCs /= bf ) CYCLE
          END IF

          n = Element % TYPE % NumberOfNodes
          Indexes(1:n) = Element % NodeIndexes

          CALL SetLumpedRows(ind,n)
        END DO

        n = COUNT( LumpedNodeSet ) 
        CALL Info('SetDirichletBoundaries','Number of lumped nodes set: '//TRIM(I2S(n)),Level=10)
      END DO

      IF( NeedListMatrix ) THEN
        DEALLOCATE( LumpedNodeSet )

        ! Revert back to CRS matrix
        CALL List_ToCRSMatrix(A)

        ! This is needed in order to copy the old BulkValues to a vector that 
        ! has the same size as the new matrix. Otherwise the matrix vector multiplication
        ! with the new Rows and Cols will fail. 
        IF( ASSOCIATED( A % BulkValues ) ) THEN
          BulkValues0 => A % BulkValues
          NULLIFY( A % BulkValues ) 
          ALLOCATE( A % BulkValues( SIZE( A % Values ) ) )
          A % BulkValues = 0.0_dp

          DO i=1,A % NumberOfRows
            DO j = Rows0(i), Rows0(i+1)-1
              k = Cols0(j) 
              DO j2 = A % Rows(i), A % Rows(i+1)-1
                k2 = A % Cols(j2)
                IF( k == k2 ) THEN
                  A % BulkValues(j2) = BulkValues0(j)
                  EXIT
                END IF
              END DO
            END DO
          END DO

          DEALLOCATE( Cols0, Rows0, BulkValues0 ) 
        END IF
        
        CALL Info('SetDircihletBoundaries','Modified matrix non-zeros: '&
            //TRIM(I2S(SIZE( A % Cols ))),Level=8)
      END IF
    END IF

    IF(.NOT.ASSOCIATED(A % DiagScaling,DiagScaling)) DEALLOCATE(DiagScaling)
    
!------------------------------------------------------------------------------

  CONTAINS

     ! Check n-t node setting element
     !-------------------------------
    SUBROUTINE CheckNTElement(n,elno)
      INTEGER :: n,elno
      INTEGER :: i,j,k,l,m,dim,kmax
      LOGICAL :: found
      REAL(KIND=dp) :: Condition(n), RotVec(3)
      
      dim = CoordinateSystemDimension()

      IF ( DOF <= 0 ) RETURN
      IF ( ALL(BoundaryReorder(Indexes(1:n))<1) ) RETURN
      IF ( .NOT. ListCheckPresent(ValueList, Name) ) RETURN
      IF ( ListGetLogical(ValueList,NormalTangentialName,Found) ) RETURN

      IF ( Conditional ) THEN
        Condition(1:n) = ListGetReal( ValueList, CondName, n, Indexes, gotIt )
        Conditional = Conditional .AND. GotIt
      END IF

      !
      ! Check for nodes belonging to n-t boundary getting set by other bcs.
      ! -------------------------------------------------------------------

      DO j=1,n
        IF ( Conditional .AND. Condition(j)<0.0_dp ) CYCLE
        k = Perm(Indexes(j))
        IF ( k > 0 ) THEN          
          k = k + OffSet
          m = BoundaryReorder(Indexes(j))
          IF ( m>0 ) THEN
            RotVec = 0._dp
            RotVec(DOF) = 1._dp
            CALL RotateNTSystem( RotVec, Indexes(j) )
            kmax = 1
            DO k=1,dim
              IF ( ABS(RotVec(k)) > ABS(RotVec(kmax)) ) kmax = k
            END DO
            NTelement(m,kmax)=elno
          END IF
        END IF
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE CheckNTElement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to a specific boundary or bulk element.
!------------------------------------------------------------------------------
    SUBROUTINE SetElementValues(n,elno)
!------------------------------------------------------------------------------
      INTEGER :: n,elno
      INTEGER :: i,j,k,l,m,dim,kmax,lmax
      LOGICAL :: CheckNT,found
      REAL(KIND=dp) :: Condition(n), Work(n), RotVec(3)
      
      dim = CoordinateSystemDimension()

      IF ( DOF > 0 ) THEN
        Work(1:n)  = ListGetReal( ValueList, Name, n, Indexes, gotIt )
        IF ( .NOT. GotIt ) THEN
          Work(1:n)  = ListGetReal( ValueList, Name(1:nlen) // ' DOFs', n, Indexes, gotIt )
        END IF
      ELSE
        CALL ListGetRealArray( ValueList, Name, WorkA, n, Indexes, gotIt )
      END IF
      
      IF ( gotIt ) THEN
        IF ( Conditional ) THEN
          Condition(1:n) = ListGetReal( ValueList, CondName, n, Indexes, gotIt )
          Conditional = Conditional .AND. GotIt
        END IF

       !
       ! Check for nodes belonging to n-t boundary getting set by other bcs.
       ! -------------------------------------------------------------------
        CheckNT = .FALSE.
        IF ( NormalTangentialNOFNodes>0 .AND. DOF>0 ) THEN
          CheckNT = .TRUE.
          IF ( ALL(BoundaryReorder(Indexes(1:n))<1) ) CheckNT = .FALSE.
          IF ( ListGetLogical(ValueList,NormalTangentialName,Found)) CheckNT=.FALSE.
        END IF
        
        DO j=1,n
          IF ( Conditional .AND. Condition(j) < 0.0d0 ) CYCLE

          k = Perm(Indexes(j))
          IF ( k > 0 ) THEN
            
            IF ( DOF>0 ) THEN
              m = 0
              IF ( NormalTangentialNOFNodes>0 ) m=BoundaryReorder(Indexes(j))
              IF ( m>0 .AND. CheckNT ) THEN
                RotVec = 0._dp
                RotVec(DOF) = 1._dp
                CALL RotateNTSystem( RotVec, Indexes(j) )

                ! When cartesian component "DOF" is defined set the N-T component
                ! closest to its direction. 
                kmax = 1 
                DO k=2,dim
                  IF ( ABS(RotVec(k)) > ABS(RotVec(kmax)) ) THEN
                    kmax = k
                  END IF
                END DO

                lmax = NDOFs * (Perm(Indexes(j))-1) + kmax
                IF ( .NOT. NTZeroing_done(m,kmax) ) THEN
                  NTZeroing_done(m,kmax) = .TRUE.
                  b(lmax) = 0._dp

                  IF( .NOT. OffDiagonal ) THEN
                    b(lmax) = b(lmax) + Work(j)/DiagScaling(lmax)
                  END IF

                  ! Consider all components of the cartesian vector mapped to the 
                  ! N-T coordinate system. Should this perhaps have scaling included?
                  CALL ZeroRow( A,lmax )
                  IF( .NOT. OffDiagonal) THEN
                    DO k=1,dim
                      l = NDOFs * (Perm(Indexes(j))-1) + k
                      CALL SetMatrixElement( A,lmax,l,RotVec(k) )
                    END DO
                  END IF
                  NTZeroing_done(m,kmax)   = .TRUE.
                  A % ConstrainedDOF(lmax) = .FALSE.
                END IF
              ELSE
                CALL SetSinglePoint(k,DOF,Work(j),.FALSE.)
              END IF
            ELSE
              DO l=1,MIN( NDOFs, SIZE(WorkA,1) )
                CALL SetSinglePoint(k,l,WorkA(l,1,j),.FALSE.)
              END DO
            END IF
          END IF
        END DO
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SetElementValues
!------------------------------------------------------------------------------
  


!------------------------------------------------------------------------------
!> Set values related to a specific boundary or bulk element.
!> If scaling has been applied the rows need to be scaled when
!> they are moved.
!------------------------------------------------------------------------------
    SUBROUTINE SetLumpedRows(ind0,n)
!------------------------------------------------------------------------------
      INTEGER :: ind0,n
      INTEGER :: ind,i,j,k,k0
      REAL(KIND=dp) :: Coeff
      ! -------------------------------------------------------------------        

      
      DO j=1,n
        ind = Indexes(j)

        IF( LumpedNodeSet(ind) ) CYCLE
        LumpedNodeSet(ind) = .TRUE.

        IF ( DOF > 0 ) THEN
          k0 = Offset + NDOFs * (Perm(ind0)-1) + DOF
          k = OffSet + NDOFs * (Perm(ind)-1) + DOF

          Coeff = DiagScaling(k0) / DiagScaling(k)

          CALL MoveRow( A, k, k0, Coeff )
          b(k0) = b(k0) + Coeff * b(k)

          CALL AddToMatrixElement( A, k, k, 1.0_dp )
          CALL AddToMatrixElement( A, k, k0, -Coeff )
          b(k) = 0.0_dp
        ELSE
          DO l = 1, NDOFs
            k0 = Offset + NDOFs + (Perm(ind0)-1) * DOF
            k = OffSet + NDOFs * (Perm(ind)-1) + l

            Coeff = DiagScaling(k0) / DiagScaling(k)

            CALL MoveRow( A, k, k0, Coeff )
            b(k0) = b(k0) + Coeff * b(k)
          
            CALL AddToMatrixElement( A, k, k, 1.0_dp )
            CALL AddToMatrixElement( A, k, k0, -Coeff )
          END DO
        END IF
      END DO

!------------------------------------------------------------------------------
    END SUBROUTINE SetLumpedRows
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Set values related to individual points.
!------------------------------------------------------------------------------
    SUBROUTINE SetPointValues(n)
!------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp) :: Work(n), Condition(n)        

      INTEGER :: i,j,k,k1,l

      IF ( DOF > 0 ) THEN
        Work(1:n) = ListGetReal( ValueList, Name, n, NodeIndexes, gotIt )
      ELSE
        CALL ListGetRealArray( ValueList, Name, WorkA, n, NodeIndexes, gotIt )
      END IF

      IF ( gotIt ) THEN

        Condition(1:n) = 1.0d0
        IF ( Conditional ) THEN
          Condition(1:n) = ListGetReal( ValueList, CondName, n, NodeIndexes, gotIt )
          Conditional = Conditional .AND. GotIt
        END IF

        DO j=1,n
          k = Perm(NodeIndexes(j))
          IF( k == 0 ) CYCLE

          IF ( Conditional .AND. Condition(j) < 0.0d0 ) CYCLE
          IF ( NodeIndexes(j) > SIZE(Perm) .OR. NodeIndexes(j) < 1 ) THEN
            CALL Warn('SetDirichletBoundaries','Invalid Node Number')
            CYCLE
          END IF

          IF ( DOF>0 ) THEN
            CALL SetSinglePoint(k,DOF,Work(j),.FALSE.)
          ELSE
            DO l=1,MIN( NDOFs, SIZE(Worka,1) )
              CALL SetSinglePoint(k,l,WorkA(l,1,j),.FALSE.)
            END DO
          END IF

        END DO
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SetPointValues
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to one single point
!------------------------------------------------------------------------------
    SUBROUTINE SetSinglePoint(ind,DOF,val,ApplyPerm)
!------------------------------------------------------------------------------
      LOGICAL :: ApplyPerm
      INTEGER :: ind, DOF
      REAL(KIND=dp) :: val

      REAL(KIND=dp) :: s
      INTEGER :: i,j,k,k1,l

      k = ind
      IF (ApplyPerm) k = Perm(ind)
      IF( k == 0 ) RETURN
      
      k = OffSet + NDOFs * (k-1) + DOF

      ! Do not add non-zero entries to pure halo nodes which are not associated with the partition.
      ! These are nodes could be created by the -halobc flag in ElmerGrid.
      IF( ParEnv % PEs > 1 ) THEN
        IF( .NOT. ANY( A % ParallelInfo % NeighbourList(k) % Neighbours == ParEnv % MyPe ) ) THEN
           RETURN
        END IF
      END IF

      IF ( A % FORMAT == MATRIX_SBAND ) THEN
        CALL SBand_SetDirichlet( A,b,k,s*val )
      ELSE
        IF (.NOT.A % NoDirichlet ) CALL ZeroRow( A,k )

        IF( .NOT. OffDiagonal ) THEN
          IF(ALLOCATED(A % Dvalues)) A % Dvalues(k) =  val / DiagScaling(k)
          IF( .NOT.A % NoDirichlet ) THEN
            CALL SetMatrixElement( A,k,k,1._dp )
            IF(.NOT.ALLOCATED(A % Dvalues)) b(k) = val / DiagScaling(k)
          END IF
        END IF
        IF(ALLOCATED(A % ConstrainedDOF)) A % ConstrainedDOF(k) = .TRUE.
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SetSinglePoint
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to upper and lower limiters.
!------------------------------------------------------------------------------
    SUBROUTINE SetLimiterValues(n)
!------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp) :: Work(n)

      Work(1:n)  = ListGetReal( ValueList, CondName, n, NodeIndexes, gotIt )

      IF ( gotIt ) THEN
        DO j=1,n
          k = Perm(NodeIndexes(j))
          IF( k == 0 ) CYCLE

          IF( .NOT. LimitActive(nDofs*(k-1)+dof)) CYCLE
          CALL SetSinglePoint(k,DOF,Work(j),.FALSE.)
        END DO
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SetLimiterValues
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> At first pass sum together the rows related to the periodic dofs.
!------------------------------------------------------------------------------
   SUBROUTINE SetPeriodicBoundariesPass1( Model, A, b, &
                      Name, DOF, NDOFs, Perm, This, Done )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), POINTER :: A  !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    LOGICAL :: Done(:)            !< Has the node already been done. 
    INTEGER :: This               !< Number of the current boundary.
    INTEGER :: DOF                !< The order number of the dof
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
                                  !< beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    INTEGER :: p,q,i,j,k,l,m,n,nn,ii,nlen,jmp,size0
    INTEGER, POINTER :: PPerm(:)
    LOGICAL :: GotIt, Found, Jump
    REAL(KIND=dp) :: Scale, weight, coeff
    TYPE(Matrix_t), POINTER :: F, G, Projector, Projector1
    TYPE(Variable_t), POINTER :: Var, WeightVar
    TYPE(ValueList_t), POINTER :: BC
    TYPE(MortarBC_t), POINTER :: MortarBC
!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)
    BC => Model % BCs(This) % Values

    IF ( ListGetLogical( BC,& 
        'Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      Scale = -1.0_dp
    ELSE IF ( ListGetLogical( BC, &
        'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      Scale = 1.0d0
    ELSE 
      Scale = ListGetConstReal( BC, &
          'Periodic BC Scale ' // Name(1:nlen), GotIt) 
      IF(.NOT. GotIt ) RETURN      
    END IF
    
    Projector => Model % BCs(This) % PMatrix
    IF ( .NOT. ASSOCIATED(Projector) ) RETURN
    
!   For explicit conditions just create the dependency almost like a normal Dirichlet BC, 
!   For implicit one (otherwise) do the assembly of the projector:
!   ---------------------------------
    IF ( ListGetLogical( BC, &
        'Periodic BC Explicit', Found ) ) THEN
      
      Var => VariableGet( Model % Variables,Name(1:nlen) ) 
      
      DO i=1,Projector % NumberOfRows
        ii = Projector % InvPerm(i)
        IF( ii == 0 ) CYCLE
        k = Perm(ii)
        IF ( .NOT. Done(ii) .AND. k>0 ) THEN
          k = NDOFs * (k-1) + DOF
          IF( .NOT.A % NoDirichlet ) THEN
            CALL ZeroRow( A,k )
            CALL AddToMatrixElement( A, k, k, 1.0_dp )
          ELSE
          END IF
          IF(ALLOCATED(A % Dvalues)) A % Dvalues(k) = 0._dp
          IF(ALLOCATED(A % ConstrainedDOF)) A % ConstrainedDOF(k) = .TRUE.
          
          DO l = Projector % Rows(i), Projector % Rows(i+1)-1
            IF ( Projector % Cols(l) <= 0 ) CYCLE
            m = Perm( Projector % Cols(l) )
            IF ( m > 0 ) THEN
              m = NDOFs * (m-1) + DOF
              IF(ALLOCATED(A % Dvalues)) THEN
                A % Dvalues(k) = A % Dvalues(k) - Scale * Projector % Values(l) * &
                    Var % Values(m)/DiagScaling(k)
              ELSE
                b(k) = b(k) - Scale * Projector % Values(l) * &
                    Var % Values(m)/DiagScaling(k)
              END IF
            END IF
          END DO
        END IF
      END DO
      
    ELSE IF ( ListGetLogical( BC, &
        'Periodic BC Use Lagrange Coefficient', Found ) ) THEN

      Jump = ListCheckPresent( BC, &
          'Periodic BC Coefficient '//Name(1:nlen))
      
      IF( .NOT. ASSOCIATED( Model % Solver % MortarBCs ) ) THEN
        CALL Info('SetPeriodicBoundariesPass1',&
            'Allocating mortar BCs for solver',Level=5)
        ALLOCATE( Model % Solver % MortarBCs( Model % NumberOfBCs ) )
        DO i=1, Model % NumberOfBCs
          Model % Solver % MortarBCs(i) % Projector => NULL()
        END DO
      END IF
      
      IF( ASSOCIATED( Projector, &
          Model % Solver % MortarBCs(This) % Projector) ) THEN
        CALL Info('SetPeridociBoundariesPass1','Using existing projector: '&
            //TRIM(I2S(This)),Level=8)
        RETURN
      END IF
      
      Model % Solver % MortarBCs(This) % Projector => Projector
      CALL Info('SetPeridociBoundariesPass1','Using projector as mortar constraint: '&
          //TRIM(I2S(This)),Level=5)

      MortarBC => Model % Solver % MortarBCs(This)      
      IF( Jump ) THEN
        IF( ASSOCIATED( MortarBC % Diag ) ) THEN
          IF( SIZE( MortarBC % Diag ) < NDofs * Projector % NumberOfRows ) THEN
            DEALLOCATE( MortarBC % Diag ) 
          END IF
        END IF
        IF( .NOT. ASSOCIATED( MortarBC % Diag ) ) THEN
          CALL Info('SetWeightedPeridocBCsPass1','Allocating projector mortar diag',Level=10)
          ALLOCATE( MortarBC % Diag( NDofs * Projector % NumberOfRows ) )
          MortarBC % Diag = 0.0_dp
        ELSE
          MortarBC % Diag( DOF::NDOFs ) = 0.0_dp
        END IF

        IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
          IF( SIZE( MortarBC % Rhs ) < NDofs * Projector % NumberOfRows ) THEN
            DEALLOCATE( MortarBC % Rhs ) 
          END IF
        END IF
        IF( .NOT. ASSOCIATED( MortarBC % Rhs ) ) THEN
          CALL Info('SetWeightedProjectorPass1','Allocating projector mortar rhs',Level=10)
          ALLOCATE( MortarBC % Rhs( NDofs * Projector % NumberOfRows ) )
          MortarBC % Rhs = 0.0_dp
        ELSE
          MortarBC % Rhs( DOF::NDOFs ) = 0.0_dp
        END IF
      END IF

      ! Create the permutation that is later need in putting the diag and rhs to correct position
      IF( ASSOCIATED( MortarBC % Perm ) ) THEN
        IF( SIZE( MortarBC % Perm ) < SIZE( Perm ) ) THEN
          DEALLOCATE( MortarBC % Perm ) 
        END IF
      END IF
      IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN
        CALL Info('SetWeightedProjectorPass1','Allocating projector mortar perm',Level=10)
        ALLOCATE( MortarBC % Perm( SIZE( Perm ) ) )
      END IF
      
      MortarBC % Perm = 0
      DO i=1,SIZE( Projector % InvPerm )
        j = Projector % InvPerm(i) 
        IF( j > 0 .AND. j <= SIZE( Perm ) ) THEN
          MortarBC % Perm( j ) = i
        END IF
      END DO
      
      ! We can use directly the nodal projector
      MortarBC % Projector => Projector
      MortarBC % SlaveScale = -Scale
      MortarBC % MasterScale = -1.0_dp
 
      IF( Jump ) THEN
        PPerm => Perm
        CALL CalculateNodalWeights(Model % Solver,.TRUE.,&
            PPerm,'Periodic BC Coefficient '//Name(1:nlen),WeightVar )
        IF(.NOT. ASSOCIATED( WeightVar ) ) THEN
          CALL Fatal('SetPeriodicBoundariesPass1','Nodal weights needed for setting jumps!')
        END IF
        
        DO i=1,Projector % NumberOfRows
          k = Projector % InvPerm(i)
          IF ( k<=0 ) CYCLE
          
          ! Add the diagonal unity projector (scaled)
          weight = WeightVar % Values( PPerm( k ) )
          coeff = ListGetRealAtNode( BC,'Periodic BC Coefficient '&
              //Name(1:nlen), k, Found )

          ! For Nodal projector the entry is 1/(weight*coeff)
          ! For Galerkin projector the is weight/coeff 
          IF( Found ) THEN
            MortarBC % Diag( NDOFS* (i-1) + DOF ) = 1.0_dp / ( weight * coeff ) 
          END IF
        END DO
      END IF

      Model % Solver % MortarBCsChanged = .TRUE.
      
    ELSE

      ALLOCATE(F)
      F = A
      F % RHS => F % BulkRHS
      F % Values => F % BulkValues

      DO i=1,Projector % NumberOfRows
         ii = Projector % InvPerm(i)
         IF( ii == 0 ) CYCLE
         k = Perm(ii)
         IF ( .NOT. Done(ii) .AND. k>0 ) THEN
            k = NDOFs * (k-1) + DOF
            DO l=Projector % Rows(i),Projector % Rows(i+1)-1
              IF ( Projector % Cols(l) <= 0 .OR. Projector % Values(l)==0.0d0 ) CYCLE

              m = Perm(Projector % Cols(l))
              IF ( m > 0 ) THEN
                m = NDOFs*(m-1) + DOF
                DO nn=A % Rows(k),A % Rows(k+1)-1
                   CALL AddToMatrixElement( A, m, A % Cols(nn), &
                          -scale*Projector % Values(l) * A % Values(nn) / &
                          DiagScaling(k) * DiagScaling(m))
                   IF (ASSOCIATED(F % Values)) THEN
                     CALL AddToMatrixElement( F, m, F % Cols(nn), &
                          -scale*Projector % Values(l) * F % Values(nn) )
                   END IF
                END DO
                b(m)=b(m) - scale*Projector % Values(l)*b(k)*&
                        DiagScaling(m) / DiagScaling(k)
                IF (ASSOCIATED(F % RHS)) THEN
                  F % RHS(m) = F % RHS(m) - scale*Projector % Values(l)*F % RHS(k)
                END IF
              END IF
            END DO
         END IF
         Done(ii) = .TRUE.
      END DO
      DEALLOCATE(F)
    END IF

!------------------------------------------------------------------------------
   END SUBROUTINE SetPeriodicBoundariesPass1
!------------------------------------------------------------------------------


!> At second pass add the [...1 .. -1 ...] row structure that results to the 
!> equality of the periodic dofs. 
!------------------------------------------------------------------------------
   SUBROUTINE SetPeriodicBoundariesPass2( Model, A, b, &
                      Name, DOF, NDOFs, Perm, This, Done )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), POINTER :: A  !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    LOGICAL :: Done(:)            !< Has the node already been done. 
    INTEGER :: This               !< Number of the current boundary.
    INTEGER :: DOF                !< The order number of the dof
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,m,n,nn,ii,nlen
    LOGICAL :: GotIt, Jump, Found
    REAL(KIND=dp) :: Scale,s,ValueOffset,val,coeff,weight
    TYPE(Matrix_t), POINTER :: Projector
    INTEGER, POINTER :: PPerm(:)
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Variable_t), POINTER :: WeightVar
!------------------------------------------------------------------------------


    BC =>  Model % BCs(This) % Values
    IF ( ListGetLogical( BC, &
         'Periodic BC Use Lagrange Coefficient', GotIt ) ) RETURN

    IF ( ListGetLogical( BC, &
         'Periodic BC Explicit', GotIt ) ) RETURN

    nlen = LEN_TRIM(Name)

    IF ( ListGetLogical( BC, &
       'Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      Scale = -1.0_dp
    ELSE IF ( ListGetLogical( BC, &
        'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      Scale = 1.0d0
    ELSE 
      Scale = ListGetCReal( BC, &
          'Periodic BC Scale ' // Name(1:nlen), GotIt) 
      IF(.NOT. GotIt ) RETURN      
    END IF

    ValueOffset = ListGetCReal( BC, &
          'Periodic BC Offset ' // Name(1:nlen), GotIt) 

    Jump = ListCheckPresent( BC, &
        'Periodic BC Coefficient '//Name(1:nlen))
    IF( Jump ) THEN
      PPerm => Perm
      CALL CalculateNodalWeights(Model % Solver,.TRUE.,&
          PPerm,'Periodic BC Coefficient '//Name(1:nlen),WeightVar )
      IF(.NOT. ASSOCIATED( WeightVar ) ) THEN
        CALL Fatal('SetPeriodicBoundariesPass1','Nodal weights needed for setting jumps!')
      END IF
    END IF

    Projector => Model % BCs(This) % PMatrix
    IF ( .NOT. ASSOCIATED(Projector) ) RETURN


!   Do the assembly of the projector:
!   ---------------------------------
    DO i=1,Projector % NumberOfRows
       ii = Projector % InvPerm(i)
       IF( ii == 0 ) CYCLE

       k = Perm( ii )
       IF ( .NOT. Done(ii) .AND. k > 0 ) THEN

         IF( Jump ) THEN
           weight = WeightVar % Values( k )
           coeff = ListGetRealAtNode( BC,'Periodic BC Coefficient '&
               //Name(1:nlen),ii, Found )
           val = -weight * coeff * DiagScaling(k)**2
           scale = -1.0
         ELSE         
           val = 1.0_dp
         END IF

          k = NDOFs * (k-1) + DOF
          IF(.NOT. Jump) THEN
            CALL ZeroRow( A,k )
            b(k) = 0.0_dp
          END IF

          DO l = Projector % Rows(i), Projector % Rows(i+1)-1
             IF ( Projector % Cols(l) <= 0 ) CYCLE
             m = Perm( Projector % Cols(l) )
             IF ( m > 0 ) THEN
               m = NDOFs * (m-1) + DOF
               CALL AddToMatrixElement( A,k,m,val * Projector % Values(l) * &
                   ( DiagScaling(m) / DiagScaling(k) ) )
             END IF
          END DO

          b(k) = b(k) - ValueOffset / DiagScaling(k)
          CALL AddToMatrixElement( A,k,k,scale*val )
       END IF
       Done(ii) = .TRUE.
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE SetPeriodicBoundariesPass2
!------------------------------------------------------------------------------




!> At first pass sum together the rows related to the periodic dofs.
!------------------------------------------------------------------------------
   SUBROUTINE SetFrictionBoundaries( Model, A, b, &
                      Name, NDOFs, Perm )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), POINTER :: A  !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
                                  !< beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    INTEGER :: t,u,j,k,k2,l,l2,n,bc_id,nlen,NormalInd
    LOGICAL :: Found
    REAL(KIND=dp),ALLOCATABLE :: Coeff(:)
    LOGICAL, ALLOCATABLE :: NodeDone(:)
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
!------------------------------------------------------------------------------


    CALL Info('SetFrictionBoundaries','Setting friction boundaries for variable: '//TRIM(Name),&
        Level=8)

    IF( NDOFs /= 2 ) THEN
      CALL Warn('SetFrictionBoundaries','Assumes friction only in 2D system')
    END IF

    nlen = LEN_TRIM(Name)
    Mesh => Model % Mesh

    ALLOCATE( NodeDone( SIZE( Perm ) ) )
    ALLOCATE( Coeff( Mesh % MaxElementNodes ) )

    NodeDone = .FALSE.
    Coeff = 0.0_dp
    
    DO t = Mesh % NumberOfBulkElements+1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t)
      
      Model % CurrentElement => Element
            
      DO bc_id = 1,Model % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) EXIT
      END DO
      IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE     
      BC => Model % BCs(bc_id) % Values

      IF( .NOT. ListGetLogical( BC,& 
          'Friction BC ' // Name(1:nlen), Found ) ) CYCLE

      NodeIndexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes
      
      Coeff(1:n) = ListGetReal( BC,& 
          'Friction Coefficient ' // Name(1:nlen), n, NodeIndexes )
      IF( ListGetLogical( BC,& 
          'Normal-Tangential ' // Name(1:nlen), Found ) ) THEN
        NormalInd = 1 
      ELSE
        NormalInd = ListGetInteger( BC,& 
            'Friction Normal Component ' // Name(1:nlen) )
      END IF

      DO i = 1, n
        j = Perm( Nodeindexes(i) )
        IF( NodeDone( j ) ) CYCLE

        k = NDOFs * (j-1) + NormalInd
        k2 = NDOFs * (j-1) + ( 3 - NormalInd ) 

        DO l = A % Rows(k),A % Rows(k+1)-1
          DO l2 = A % Rows(k2), A % Rows(k2+1)-1
            IF( A % Cols(l2) == A % Cols(l) ) EXIT
          END DO
          A % Values(l2) = A % Values(l2) - Coeff(i) * A % Values(l)
        END DO
        A % Rhs(k2) = A % Rhs(k2) - Coeff(i) * A % Rhs(k)
        NodeDone( j ) = .TRUE.
      END DO
    END DO

    n = COUNT( NodeDone ) 
    CALL Info('SetFrictionBoundaries','Number of friction nodes: '//TRIM(I2S(n)),Level=10)

    DEALLOCATE( NodeDone, Coeff )

!------------------------------------------------------------------------------
  END SUBROUTINE SetFrictionBoundaries
!------------------------------------------------------------------------------


!> Set the diagonal entry related to mortar BCs.
!> This implements the implicit jump condition. 
!------------------------------------------------------------------------------
   SUBROUTINE SetWeightedProjectorJump( Model, A, b, &
       Name, DOF, NDOFs, Perm )
     !------------------------------------------------------------------------------
     TYPE(Model_t) :: Model        !< The current model structure
     TYPE(Matrix_t), POINTER :: A  !< The global matrix
     REAL(KIND=dp) :: b(:)         !< The global RHS vector
     CHARACTER(LEN=*) :: Name      !< name of the dof to be set
     INTEGER :: DOF                !< The order number of the dof
     INTEGER :: NDOFs              !< the total number of DOFs for this equation
     INTEGER, TARGET :: Perm(:)    !< The node reordering info
     !------------------------------------------------------------------------------
     INTEGER :: i,j,k,i2,j2,k2,n,u,v,node,totsize,nodesize,bc_ind,&
         nnodes,nlen,TargetBC
     INTEGER, POINTER :: PPerm(:)
     INTEGER, ALLOCATABLE :: InvInvPerm(:)
     LOGICAL :: Found, AddRhs, AddCoeff, AddRes
     LOGICAL, ALLOCATABLE :: NodeDone(:)
     REAL(KIND=dp) :: coeff, weight, voff, res
     TYPE(Matrix_t), POINTER :: Projector
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Element_t), POINTER :: Element, Left, Right
     LOGICAL :: SomethingDone
     TYPE(MortarBC_t), POINTER :: MortarBC
     !------------------------------------------------------------------------------

     ! If there is no mortar projector then nothing to do
     SomethingDone = .FALSE.

     ! Go through the projectors and check for jumps
     ! If there is a jump add an entry to the diagonal-to-be
     DO bc_ind=1,Model % NumberOFBCs

       MortarBC => Model % Solver % MortarBCs(bc_ind) 

       Projector => MortarBC % Projector
       IF( .NOT. ASSOCIATED( Projector ) ) CYCLE

       ! For this boundary there should also be a coefficient 
       ! otherwise nothing needs to be done. 
       nlen = LEN_TRIM(Name)
       BC => Model % BCs(bc_ind) % Values

       AddCoeff = ListCheckPresent( BC,'Mortar BC Coefficient '//Name(1:nlen))
       AddRes = ListCheckPresent( BC,'Mortar BC Resistivity '//Name(1:nlen))
       AddRhs = ListCheckPresent( BC,'Mortar BC Offset '//Name(1:nlen))

       IF( .NOT. (AddCoeff .OR. AddRes .OR. AddRhs) ) CYCLE

       Model % Solver % MortarBCsChanged = .TRUE.

       IF( .NOT. ASSOCIATED( Projector % InvPerm ) ) THEN
         CALL Fatal('SetWeightedProjectorJump','The > Projector % InvPerm < is really needed here!')
       END IF

       totsize = MAXVAL( Perm )
       nodesize = MAXVAL( Perm(1:Model % Mesh % NumberOfNodes ) )

       IF( AddCoeff .OR. AddRes ) THEN
         IF( ASSOCIATED( MortarBC % Diag ) ) THEN
           IF( SIZE( MortarBC % Diag ) < NDofs * Projector % NumberOfRows ) THEN
             DEALLOCATE( MortarBC % Diag ) 
           END IF
         END IF
         IF( .NOT. ASSOCIATED( MortarBC % Diag ) ) THEN
           CALL Info('SetWeightedProjectorJump','Allocating projector mortar diag',Level=10)
           ALLOCATE( MortarBC % Diag( NDofs * Projector % NumberOfRows ) )
           MortarBC % Diag = 0.0_dp
         ELSE
           MortarBC % Diag(DOF::NDOFs) = 0.0_dp
         END IF
       END IF

       IF( AddRhs ) THEN
         IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
           IF( SIZE( MortarBC % Rhs ) < NDofs * Projector % NumberOfRows ) THEN
             DEALLOCATE( MortarBC % Rhs ) 
           END IF
         END IF
         IF( .NOT. ASSOCIATED( MortarBC % Rhs ) ) THEN
           CALL Info('SetWeightedProjectorJump','Allocating projector mortar rhs',Level=10)
           ALLOCATE( MortarBC % Rhs( NDofs * Projector % NumberOfRows ) )
           MortarBC % Rhs = 0.0_dp
         ELSE
           MortarBC % Rhs(DOF::NDOFs) = 0.0_dp
         END IF
       END IF

       ! Create the permutation that is later need in putting the diag and rhs to correct position
       IF( ASSOCIATED( MortarBC % Perm ) ) THEN
         IF( SIZE( MortarBC % Perm ) < SIZE( Perm ) ) THEN
           DEALLOCATE( MortarBC % Perm ) 
         END IF
       END IF
       IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN
         CALL Info('SetWeightedProjectorJump','Allocating projector mortar perm',Level=10)
         ALLOCATE( MortarBC % Perm( SIZE( Perm ) ) )
       END IF

       MortarBC % Perm = 0
       DO i=1,SIZE( Projector % InvPerm )
         j = Projector % InvPerm(i) 
         IF( j > 0 .AND. j <= nodesize ) THEN
           MortarBC % Perm( j ) = i
         END IF
       END DO


       TargetBC = ListGetInteger( BC,'Mortar BC',Found ) 

       CALL Info('SetWeightedProjectorJump','Setting jump to mortar projector in BC '&
           //TRIM(I2S(bc_ind)),Level=7)
    
       ! Create a table that shows how the additional degrees of freedom map
       ! to their corresponding regular dof. This is needed when creating the jump.
       ALLOCATE( NodeDone( Projector % NumberOfRows ) )
       NodeDone = .FALSE.
       
       ! Looping through elements rather than looping through projector rows directly
       ! is done in order to be able to refer to boundary properties associated 
       ! with the element. 
       DO t=1,Model % Mesh % NumberOfBoundaryElements
         Element => Model % Mesh % Elements( t + Model % Mesh % NumberOfBulkElements )

         IF( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE

         ! Outside code this tells the active element
         Model % CurrentElement => Element

         Left => Element % BoundaryInfo % Left
         Right => Element % BoundaryInfo % Right 
        
         IF( TargetBC > 0 ) THEN
           IF( ASSOCIATED( Left ) ) THEN
             IF( Left % PartIndex /= ParEnv % myPE ) CYCLE
           ELSE IF ( ASSOCIATED( Right ) ) THEN
             IF( Left % PartIndex /= ParEnv % myPE ) CYCLE
           ELSE
             CYCLE
           END IF
         ELSE
           ! This case is for the case when TargetBC = 0 i.e. for Discontinuous BC
           ! These are conditions that resulted to creation of zero 
           ! constraint matrix entries in this partition so no need to do them.
           IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) THEN
             CYCLE
           END IF
           
           ! For this we have a zero mass matrix entry so don't bother to add zero
!          IF( Left % PartIndex /= ParEnv % myPE .AND. &
!              Right % PartIndex /= ParEnv % myPe ) THEN
!            CYCLE
!          END IF
         END IF

         nnodes = Element % TYPE % NumberOfNodes
         DO u=1, nnodes
           node = Element % NodeIndexes(u)

           IF( Perm( node ) == 0 ) CYCLE

           i = MortarBC % Perm( node ) 
           IF( i == 0 ) CYCLE

           IF( NodeDone( i ) ) CYCLE
           NodeDone( i ) = .TRUE. 

           Found = .FALSE.

           IF( AddCoeff ) THEN
             coeff = ListGetRealAtNode( BC,'Mortar BC Coefficient '&
                 //Name(1:nlen),node, Found )        
             res = 1.0_dp / coeff
           END IF

           IF( AddRes ) THEN
             res = ListGetRealAtNode( BC,'Mortar BC Resistivity '&
                 //Name(1:nlen),node, Found )
           END IF

           ! For Nodal projector the entry is 1/(weight*coeff)
           ! For Galerkin projector the is weight/coeff 
           IF( Found ) THEN 
             IF( AddCoeff .OR. Addres ) THEN
               MortarBC % Diag(NDOFs*(i-1)+DOF) = res
             END IF
           END IF

           IF( AddRhs ) THEN
             voff = ListGetRealAtNode( BC,'Mortar BC Offset '&
                 //Name(1:nlen),node, Found )        
             IF( Found ) THEN
               MortarBC % Rhs(NDofs*(i-1)+DOF) = voff
             END IF
           END IF

         END DO
       END DO
       
       SomethingDone = .TRUE.

       DEALLOCATE( NodeDone )
     END DO

     IF( SomethingDone ) THEN
       CALL Info('setWeightedProjectorJump','Created a jump for weighted projector',Level=7)
     END IF
 
!------------------------------------------------------------------------------
   END SUBROUTINE SetWeightedProjectorJump
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichletBoundaries
!------------------------------------------------------------------------------



!> Prepare to set Dirichlet conditions for attachment DOFs in the case of
!> component mode synthesis
!------------------------------------------------------------------------------
  SUBROUTINE SetConstraintModesBoundaries( Model, A, b, &
      Name, NDOFs, Perm )
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), POINTER :: A  !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
                                  !< beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    INTEGER :: i,t,u,j,k,k2,l,l2,n,bc_id,nlen,NormalInd
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: Var
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER, ALLOCATABLE :: BCPerm(:)

!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)
    Mesh => Model % Mesh
    Solver => Model % Solver
    Var => Solver % Variable 
    
    ! This needs to be allocated only once, hence return if already set
    IF( Var % NumberOfConstraintModes > 0 ) RETURN

    CALL Info('SetConstraintModesBoundaries','Setting constraint modes boundaries for variable: '&
        //TRIM(Name),Level=7)

    IF( .NOT. ListGetLogicalAnyBC( Model,'Constraint Modes ' // Name(1:nlen)) ) THEN
      CALL Warn('SetConstraintModesBoundaries',&
          'Constraint Modes Analysis requested but no constrained BCs given!')
      RETURN
    END IF

    ! Allocate the indeces for the constraint modes
    ALLOCATE( Var % ConstraintModesIndeces( A % NumberOfRows ) )
    Var % ConstraintModesIndeces = 0
    
    ALLOCATE( BCPerm( Model % NumberOfBCs ) ) 
    BCPerm = 0
    
    j = 0 
    DO bc_id = 1,Model % NumberOfBCs
      BC => Model % BCs(bc_id) % Values        
      IF( ListGetLogical( BC,& 
          'Constraint Modes ' // Name(1:nlen), Found ) ) THEN
        j = j + 1
        BCPerm(bc_id) = j
      END IF
    END DO
    CALL Info('SetConstraintModesBoundaries','Number of active constraint modes boundaries: '&
        //TRIM(I2S(j)),Level=7)


    Var % NumberOfConstraintModes = NDOFS * j 
    
    
    DO t = Mesh % NumberOfBulkElements+1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t)
      
      DO bc_id = 1,Model % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) EXIT
      END DO
      IF( bc_id > Model % NumberOfBCs ) CYCLE      
      IF( BCPerm(bc_id) == 0 ) CYCLE
      
      NodeIndexes => Element % NodeIndexes

      ! For vector valued problems treat each component as separate dof
      DO k=1,NDOFs
        Var % ConstraintModesIndeces( NDOFS*(Perm(NodeIndexes)-1)+k) = NDOFS*(BCPerm(bc_id)-1)+k
      END DO
    END DO
    
    ! The constraint modes can be either lumped or not.
    ! If they are not lumped then mark each individually
    IF( .NOT. ListGetLogical(Solver % Values,'Constraint Modes Lumped',Found ) ) THEN
      j = 0
      DO i=1,A % NumberOfRows
        IF( Var % ConstraintModesIndeces(i) > 0 ) THEN
          j = j + 1
          Var % ConstraintModesIndeces(i) = j
        END IF
      END DO
      CALL Info('SetConstraintModesBoundaries','Number of active constraint modes: '&
          //TRIM(I2S(j)),Level=7)
      Var % NumberOfConstraintModes = j 
    END IF
    

    ! Manipulate the boundaries such that we need to modify only the r.h.s. in the actual linear solver
    DO k=1,A % NumberOfRows       
      IF( Var % ConstraintModesIndeces(k) == 0 ) CYCLE
      b(k) = 0.0_dp
      CALL ZeroRow(A,k)
      CALL SetMatrixElement( A,k,k,1._dp )
    END DO

    
    ALLOCATE( Var % ConstraintModes( Var % NumberOfConstraintModes, A % NumberOfRows ) )
    Var % ConstraintModes = 0.0_dp

    DEALLOCATE( BCPerm ) 
    
    CALL Info('SetConstraintModesBoundarues','All done',Level=10)

!------------------------------------------------------------------------------
  END SUBROUTINE SetConstraintModesBoundaries
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Sets just one Dirichlet point in contrast to setting the whole field.
!> This is a lower order routine that the previous one. 
!------------------------------------------------------------------------------
  SUBROUTINE SetDirichletPoint( A, b,DOF, NDOFs, Perm, NodeIndex, NodeValue) 
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: b(:)
    REAL(KIND=dp) :: NodeValue
    INTEGER :: DOF, NDOFs, Perm(:), NodeIndex
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: s
    INTEGER :: PermIndex

!------------------------------------------------------------------------------
    PermIndex = Perm(NodeIndex)
    IF ( PermIndex > 0 ) THEN
      PermIndex = NDOFs * (PermIndex-1) + DOF

      IF ( A % FORMAT == MATRIX_SBAND ) THEN
        CALL SBand_SetDirichlet( A,b,PermIndex,NodeValue )
      ELSE IF ( A % FORMAT == MATRIX_CRS .AND. &
          A % Symmetric.AND..NOT. A % NoDirichlet ) THEN

         CALL CRS_SetSymmDirichlet(A,b,PermIndex,NodeValue)
      ELSE                  
        s = A % Values(A % Diag(PermIndex))
        b(PermIndex) = NodeValue * s
        IF(.NOT. A % NoDirichlet) THEN
           CALL ZeroRow( A,PermIndex )
           CALL SetMatrixElement( A,PermIndex,PermIndex,s )
        END IF
      END IF
      IF(ALLOCATED(A % ConstrainedDOF)) A % ConstrainedDOF(PermIndex) = .TRUE.
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichletPoint
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Sets nodal loads directly to the matrix structure. 
!> The intended use for this, is for example, in multiphysics coupling where
!> the nodal loads may have been computed by another solver. 
!------------------------------------------------------------------------------
   SUBROUTINE SetNodalLoads( Model, A, b, Name, DOF, NDOFs, Perm )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model         !< The current model structure
    TYPE(Matrix_t), POINTER :: A   !< The global matrix
    REAL(KIND=dp) :: b(:)          !< The global RHS vector
    CHARACTER(LEN=*) :: Name       !< Name of the dof to be set
    INTEGER :: DOF                 !< The order number of the dof
    INTEGER :: NDOFs               !< The total number of DOFs for this equation
    INTEGER :: Perm(:)             !< The node reordering info, this has been generated at the
                                   !< beginning of the simulation for bandwidth optimization.
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: Element
    INTEGER, ALLOCATABLE :: Indexes(:)
    INTEGER, POINTER :: NodeIndexes(:), Neigh(:)
    INTEGER :: BC,i,j,k,l,n,t,k1,k2
    LOGICAL :: GotIt
    REAL(KIND=dp), POINTER :: WorkA(:,:,:) => NULL()
    REAL(KIND=dp) ::  s

    LOGICAL :: Conditional
    CHARACTER(LEN=MAX_NAME_LEN) :: LoadName

    INTEGER, POINTER :: IndNodes(:)
    INTEGER :: NoNodes,NoDims,bf_id,nlen,NOFNodesFound
    REAL(KIND=dp), POINTER :: CoordNodes(:,:), DiagScaling(:),MinDist(:)
    REAL(KIND=dp) :: GlobalMinDist,Dist,Eps
    LOGICAL, ALLOCATABLE :: ActivePart(:), ActivePartAll(:), DoneLoad(:)
    LOGICAL :: NodesFound
    TYPE(ValueList_t), POINTER :: ValueList

    LoadName = TRIM(Name) // ' Load'
    nlen = LEN_TRIM(LoadName)
    
    CALL Info('SetNodalLoads','Checking for nodal loads for variable: '//TRIM(Name),Level=10)

    n = MAX(Model % NumberOfBCs, Model % NumberOFBodyForces) 
    ALLOCATE( ActivePart(n), ActivePartAll(n) )

    ALLOCATE( Indexes(Model % Solver % Mesh % MaxElementDOFs) )
!------------------------------------------------------------------------------
! Go through the boundaries
!------------------------------------------------------------------------------

    DiagScaling => A % DiagScaling
    IF (.NOT.ASSOCIATED(DiagScaling)) THEN
      ALLOCATE(DiagScaling(A % NumberOFRows))
      DiagScaling=1._dp
    END IF

    ActivePart = .FALSE.
    ActivePartAll = .FALSE.
    DO BC=1,Model % NumberOfBCs
      IF(.NOT. ListCheckPresent( Model % BCs(BC) % Values,'Target Boundaries')) CYCLE
      ActivePart(BC) = ListCheckPresent( Model % BCs(BC) % Values, LoadName )
      ActivePartAll(BC) = ListCheckPresent( &
          Model % BCs(BC) % Values, LoadName(1:nlen) // ' DOFs' )
    END DO

    IF ( ANY(ActivePart) .OR. ANY(ActivePartAll) ) THEN
      CALL Info('SetNodalLoads','Setting nodal loads on boundaries: '//TRIM(LoadName),Level=9)
      ALLOCATE(DoneLoad( SIZE(b)/NDOFs) )
      DoneLoad = .FALSE.

      DO BC=1,Model % NumberOfBCs
        IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE

        DO t = Model % NumberOfBulkElements + 1, &
          Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

          Element => Model % Elements(t)
          IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BC) % Tag ) CYCLE
          
          Model % CurrentElement => Element
          IF ( ActivePart(BC) ) THEN
            n = Element % TYPE % NumberOfNodes
            Indexes(1:n) = Element % NodeIndexes
          ELSE
            n = SgetElementDOFs( Indexes )
          END IF
          ValueList => Model % BCs(BC) % Values

          CALL SetElementLoads( n )
        END DO
      END DO
    END IF

!------------------------------------------------------------------------------
! Go through the nodal load conditions for the body force list
!------------------------------------------------------------------------------

    ActivePart = .FALSE.
    ActivePartAll = .FALSE.
    DO bf_id=1,Model % NumberOFBodyForces
      ActivePart(bf_id) = ListCheckPresent( Model % BodyForces(bf_id) % Values, LoadName ) 
      ActivePartAll(bf_id) = ListCheckPresent( &
            Model % BodyForces(bf_id) % Values, LoadName(1:nlen) // ' DOFs' ) 
    END DO

    IF ( ANY( ActivePart ) .OR. ANY(ActivePartAll) ) THEN
      CALL Info('SetNodalLoads','Setting nodal loads on body force: '//TRIM(LoadName),Level=9)
      IF(.NOT. ALLOCATED(DoneLoad)) ALLOCATE(DoneLoad( SIZE(b)/NDOFs) )      
      DoneLoad = .FALSE.

      DO t = 1, Model % NumberOfBulkElements 
        Element => Model % Elements(t)
        bf_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force', GotIt)
        
        IF(.NOT. GotIt) CYCLE
        IF(.NOT. ActivePart(bf_id) .AND. .NOT. ActivePartAll(bf_id) ) CYCLE

        Model % CurrentElement => Element
        IF ( ActivePart(bf_id) ) THEN
          n = Element % TYPE % NumberOfNodes
          Indexes(1:n) = Element % NodeIndexes
        ELSE
          n = SgetElementDOFs( Indexes )
        END IF
        ValueList => Model % BodyForces(bf_id) % Values

        CALL SetElementLoads( n )
      END DO
    END IF
   
    DEALLOCATE(ActivePart)
    IF(ALLOCATED(DoneLoad)) DEALLOCATE(DoneLoad)


!------------------------------------------------------------------------------
! Go through the point loads that are created on-the-fly
!------------------------------------------------------------------------------

    DO BC=1,Model % NumberOfBCs
      ValueList => Model % BCs(BC) % Values
      IF( .NOT. ListCheckPresent( ValueList,LoadName )) CYCLE
      NodesFound = ListCheckPresent(ValueList,'Target Nodes')
      
      ! At the first calling the list of coorinates is transformed to list of nodes. 
      IF(.NOT. NodesFound) THEN
        CoordNodes => ListGetConstRealArray(ValueList, 'Target Coordinates',GotIt)
        IF(GotIt) THEN
          Eps = ListGetConstReal( ValueList, 'Target Coordinates Eps', Gotit )
          IF ( .NOT. GotIt ) THEN
            Eps = HUGE(Eps)
          ELSE
            ! We are looking at square of distance
            Eps = Eps**2
          END IF

          NoNodes = SIZE(CoordNodes,1)
          NoDims = SIZE(CoordNodes,2)
          
          IF(NoNodes > 0) THEN               
            ALLOCATE( IndNodes(NoNodes), MinDist(NoNodes) )
            IndNodes = -1
            MinDist = HUGE( Dist )
            DO j=1,NoNodes
              DO i=1,Model % NumberOfNodes
                IF( Perm(i) == 0) CYCLE
                
                Dist = (Model % Mesh % Nodes % x(i) - CoordNodes(j,1))**2 
                IF(NoDims >= 2) Dist = Dist + (Model % Mesh % Nodes % y(i) - CoordNodes(j,2))**2
                IF(NoDims == 3) Dist = Dist + (Model % Mesh % Nodes % z(i) - CoordNodes(j,3))**2
                Dist = SQRT(Dist)
                
                IF(Dist < MinDist(j) .AND. Dist <= Eps ) THEN
                  MinDist(j) = Dist
                  IndNodes(j) = i
                END IF
              END DO
            END DO

            ! In parallel case eliminate all except the nearest node. 
            ! This relies on the fact that for each node partition the 
            ! distance to nearest node is computed accurately. 
            DO j=1,NoNodes
              GlobalMinDist = ParallelReduction( MinDist(j), 1 )
              IF(ABS(GlobalMinDist - MinDist(j) )>TINY(Dist)) THEN
                IndNodes(j) = 0
              ELSE IF(ParEnv % PEs>1) THEN
                ! In parallel apply load only on the owner partition:
                ! ---------------------------------------------------
                neigh=>Model % Mesh % ParallelInfo % NeighbourList(IndNodes(j)) % Neighbours
                DO i=1,SIZE(Neigh)
                  IF(ParEnv % Active(neigh(i))) EXIT
                END DO
                IF(neigh(i)/=ParEnv % MyPE) IndNodes(j) = 0
              END IF
            END DO

            NOFNodesFound = 0
            DO j=1,NoNodes
               IF ( IndNodes(j)>0 ) THEN
                 NOFNodesFound = NOFNodesFound+1
                 IndNodes(NOFNodesFound) = IndNodes(j)
               END IF
            END DO
            
            ! In the first time add the found nodes to the list structure
            IF ( NOFNodesFound > 0 ) THEN
              CALL ListAddIntegerArray( ValueList,'Target Nodes', &
                  NOFNodesFound, IndNodes) 
              NodesFound = .TRUE.            
            ELSE
              ! If no nodes found, add still an empty list and make sure the 
              ! zero is not treated later on. Otherwise this search would be 
              ! retreated each time. 
              CALL ListAddIntegerArray( ValueList,'Target Nodes', 0, IndNodes) 
            END IF

            ! Finally deallocate the temporal vectors
            DEALLOCATE( IndNodes, MinDist ) 
          END IF
        END IF
      END IF
      
      IF(NodesFound) THEN           
        CALL Info('SetNodalLoads','Setting nodal loads on target nodes: '//TRIM(Name),Level=9)
        NodeIndexes => ListGetIntegerArray( ValueList,'Target Nodes')
        n = SIZE(NodeIndexes)
        CALL SetPointLoads(n)
      END IF

    END DO

    DEALLOCATE( Indexes )
    IF(.NOT.ASSOCIATED(A % DiagScaling,DiagScaling)) DEALLOCATE(DiagScaling)

    CALL Info('SetNodalLoads','Finished checking for nodal loads',Level=12)


CONTAINS

     SUBROUTINE SetElementLoads(n)
       INTEGER :: n
       REAL(KIND=dp) :: Work(n)
       
       NodeIndexes => Element % NodeIndexes(1:n)
       
       IF ( DOF > 0 ) THEN
         Work(1:n) = ListGetReal( ValueList, LoadName, n, Indexes, gotIt )
         IF ( .NOT. Gotit ) THEN
           Work(1:n) = ListGetReal( ValueList, LoadName(1:nlen) // ' DOFs', n, Indexes, gotIt )
         END IF
       ELSE
         CALL ListGetRealArray( ValueList, LoadName, WorkA, n, Indexes, gotIt )
       END IF

       IF ( gotIt ) THEN

         DO j=1,n
           k = Perm(Indexes(j))
           
           IF ( k > 0 ) THEN
             IF ( DoneLoad(k) ) CYCLE
             DoneLoad(k) = .TRUE.

             IF ( DOF>0 ) THEN
               k = NDOFs * (k-1) + DOF
               IF( ParEnv % Pes > 1 ) THEN
                  IF(  A % ParallelInfo % NeighbourList(k) % Neighbours(1) /= ParEnv % MyPe ) CYCLE
               END IF
               b(k) = b(k) + Work(j) * DiagScaling(k)
             ELSE
               DO l=1,MIN( NDOFs, SIZE(Worka,1) )
                 k1 = NDOFs * (k-1) + l
                 b(k1) = b(k1) + WorkA(l,1,j) * DiagScaling(k1)
               END DO
             END IF
           END IF
         END DO
       END IF
       
     END SUBROUTINE SetElementLoads
     
     
     SUBROUTINE SetPointLoads(n)
       INTEGER :: n
       REAL(KIND=dp) :: Work(n)

       IF(n<=0) RETURN
       
       IF ( DOF > 0 ) THEN
         Work(1:n) = ListGetReal( ValueList, LoadName, n, NodeIndexes, gotIt )
       ELSE
         CALL ListGetRealArray( ValueList, LoadName, WorkA, n, NodeIndexes, gotIt )
       END IF
       
       IF ( GotIt ) THEN
         DO j=1,n
           IF ( NodeIndexes(j) > SIZE(Perm) .OR. NodeIndexes(j) < 1 ) THEN
             CALL Warn('SetNodalLoads','Invalid Node Number')
             CYCLE
           END IF
         
           k = Perm(NodeIndexes(j))
           IF ( k > 0 ) THEN
             IF ( DOF>0 ) THEN
               k = NDOFs * (k-1) + DOF
               b(k) = b(k) + Work(j) * DiagScaling(k)
             ELSE
               DO l=1,MIN( NDOFs, SIZE(WorkA,1) )
                 k1 = NDOFs * (k-1) + l
                 b(k1) = b(k1) + WorkA(l,1,j) * DiagScaling(k1)
               END DO
             END IF
           END IF
         END DO
       END IF

     END SUBROUTINE SetPointLoads
     
!------------------------------------------------------------------------------
   END SUBROUTINE SetNodalLoads
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION sGetElementDOFs( Indexes, UElement, USolver )  RESULT(NB)
!------------------------------------------------------------------------------
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     INTEGER :: Indexes(:)

     TYPE(Solver_t),  POINTER :: Solver
     TYPE(Element_t), POINTER :: Element, Parent

     LOGICAL :: Found, GB
     INTEGER :: nb,i,j,EDOFs, FDOFs, BDOFs,FaceDOFs, EdgeDOFs, BubbleDOFs

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     NB = 0

     IF ( ListGetLogical( Solver % Values, 'Discontinuous Galerkin', Found ) ) THEN
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

     DO i=1,Element % NDOFs
        NB = NB + 1
        Indexes(NB) = Element % NodeIndexes(i)
     END DO

     FaceDOFs   = Solver % Mesh % MaxFaceDOFs
     EdgeDOFs   = Solver % Mesh % MaxEdgeDOFs
     BubbleDOFs = Solver % Mesh % MaxBDOFs

     IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
        DO j=1,Element % TYPE % NumberOFEdges
          EDOFs = Solver % Mesh % Edges( Element % EdgeIndexes(j) ) % BDOFs
          DO i=1,EDOFs
             NB = NB + 1
             Indexes(NB) = EdgeDOFs*(Element % EdgeIndexes(j)-1) + &
                      i + Solver % Mesh % NumberOfNodes
          END DO
        END DO
     END IF

     IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
        DO j=1,Element % TYPE % NumberOFFaces
           FDOFs = Solver % Mesh % Faces( Element % FaceIndexes(j) ) % BDOFs
           DO i=1,FDOFs
              NB = NB + 1
              Indexes(NB) = FaceDOFs*(Element % FaceIndexes(j)-1) + i + &
                 Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges
           END DO
        END DO
     END IF

     GB = ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
     IF (.NOT. Found) GB = .TRUE.

     IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
       IF (.NOT. isActivePElement(Element) ) RETURN

       Parent => Element % BoundaryInfo % Left
       IF (.NOT.ASSOCIATED(Parent) ) &
         Parent => Element % BoundaryInfo % Right
       IF (.NOT.ASSOCIATED(Parent) ) RETURN

       IF ( ASSOCIATED( Parent % EdgeIndexes ) ) THEN
         EDOFs = Element % BDOFs
         DO i=1,EDOFs
           NB = NB + 1
           Indexes(NB) = EdgeDOFs*(Parent % EdgeIndexes(Element % PDefs % LocalNumber)-1) + &
                    i + Solver % Mesh % NumberOfNodes
         END DO
       END IF

       IF ( ASSOCIATED( Parent % FaceIndexes ) ) THEN
         FDOFs = Element % BDOFs
         DO i=1,FDOFs
           NB = NB + 1
           Indexes(NB) = FaceDOFs*(Parent % FaceIndexes(Element % PDefs % LocalNumber)-1) + i + &
              Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges
         END DO
       END IF
     ELSE IF ( GB ) THEN
        IF ( ASSOCIATED( Element % BubbleIndexes ) ) THEN
           DO i=1,Element % BDOFs
              NB = NB + 1
              Indexes(NB) = FaceDOFs*Solver % Mesh % NumberOfFaces + &
                 Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges + &
                   Element % BubbleIndexes(i)
           END DO
        END IF
     END IF
!------------------------------------------------------------------------------
  END FUNCTION SgetElementDOFs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if Normal / Tangential vector boundary conditions present and
!> allocate space for normals, and if in 3D for two tangent direction
!> vectors.
!------------------------------------------------------------------------------
   SUBROUTINE CheckNormalTangentialBoundary( Model, VariableName, &
     NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals,     &
        BoundaryTangent1, BoundaryTangent2, dim )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model

    CHARACTER(LEN=*) :: VariableName

    INTEGER, POINTER :: BoundaryReorder(:)
    INTEGER :: NumberOfBoundaryNodes,dim

    REAL(KIND=dp), POINTER :: BoundaryNormals(:,:),BoundaryTangent1(:,:), &
                       BoundaryTangent2(:,:)
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER :: i,j,k,n,t,ierr,iter, proc
    LOGICAL :: GotIt, Found, Conditional
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp), ALLOCATABLE :: Condition(:)

    TYPE buff_t
      INTEGER, ALLOCATABLE :: buff(:)
    END TYPE buff_t
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER, POINTER :: nlist(:)
    TYPE(Buff_t), ALLOCATABLE, TARGET :: n_index(:)
    INTEGER, ALLOCATABLE :: n_count(:), gbuff(:)
!------------------------------------------------------------------------------

    ! need an early initialization to avarage normals across partitions:
    !-------------------------------------------------------------------
    IF ( Parenv  % PEs >1 ) THEN
      IF (.NOT. ASSOCIATED(Model % Solver % Matrix % ParMatrix) ) &
         CALL ParallelInitMatrix( Model % Solver, Model % Solver % Matrix )
    END IF

    NumberOfBoundaryNodes = 0

    Found = .FALSE.
    DO i=1,Model % NumberOfBCs
      IF ( ListGetLogical(Model % BCs(i) % Values, VariableName, Gotit) ) THEN
        Found = ListGetLogical( Model % BCs(i) % Values, &
           TRIM(VariableName) // ' Rotate',Gotit )
        IF (.NOT. Gotit ) Found = .TRUE.
        IF ( Found ) EXIT
      END IF
    END DO
    IF ( .NOT. Found ) RETURN

    Mesh => Model % Mesh
    n = Mesh % NumberOFNodes

    IF ( .NOT. ASSOCIATED( BoundaryReorder ) ) THEN
      ALLOCATE( BoundaryReorder(n) )
    ELSE IF ( SIZE(BoundaryReorder)<n ) THEN
      DEALLOCATE( BoundaryReorder )
      ALLOCATE( BoundaryReorder(n) )
    END IF
    BoundaryReorder = 0

!------------------------------------------------------------------------------
    DO t=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
                  Mesh % NumberOfBoundaryElements

      CurrentElement => Model % Elements(t)
      IF ( CurrentElement % TYPE % ElementCode == 101 )  CYCLE

      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      ALLOCATE( Condition(n)  )
      DO i=1,Model % NumberOfBCs
        IF ( CurrentElement % BoundaryInfo % Constraint == &
                  Model % BCs(i) % Tag ) THEN
          IF ( ListGetLogical( Model % BCs(i) % Values,VariableName, gotIt) ) THEN
            Found = ListGetLogical( Model % BCs(i) % Values, &
                 TRIM(VariableName) // ' Rotate',gotIt)
            IF ( Found .OR. .NOT. GotIt ) THEN
              Condition(1:n) = ListGetReal( Model % BCs(i) % Values, &
                 TRIM(VariableName) // ' Condition', n, NodeIndexes, Conditional )

              DO j=1,n
                IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE

                k = NodeIndexes(j)
                IF ( BoundaryReorder(k)==0 ) THEN
                  NumberOfBoundaryNodes = NumberOfBoundaryNodes + 1
                  BoundaryReorder(k) = NumberOfBoundaryNodes
                END IF
              END DO
            END IF
          END IF
        END IF
      END DO
      DEALLOCATE( Condition )
    END DO

    IF (ParEnv % PEs>1 )  THEN
!------------------------------------------------------------------------------
!   If parallel execution, check for parallel matrix initializations
!------------------------------------------------------------------------------
      ALLOCATE( n_count(ParEnv% PEs),n_index(ParEnv % PEs) )
      n_count = 0
      IF ( NumberOfBoundaryNodes>0 ) THEN
        DO i=1,Mesh % NumberOfNodes
          IF (BoundaryReorder(i)<=0 ) CYCLE
          IF (.NOT.Mesh % ParallelInfo % INTERFACE(i) ) CYCLE

          nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
          DO j=1,SIZE(nlist)
            k = nlist(j)+1
            IF ( k-1 == ParEnv % myPE ) CYCLE
            n_count(k) = n_count(k)+1
          END DO
        END DO
        DO i=1,ParEnv % PEs
          IF ( n_count(i)>0 ) &
            ALLOCATE( n_index(i) % buff(n_count(i)) )
        END DO
        n_count = 0
        DO i=1,Mesh % NumberOfNodes
          IF (BoundaryReorder(i)<=0 ) CYCLE
          IF (.NOT.Mesh % ParallelInfo % INTERFACE(i) ) CYCLE

          nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
          DO j=1,SIZE(nlist)
            k = nlist(j)+1
            IF ( k == ParEnv % myPE+1 ) CYCLE
            n_count(k) = n_count(k)+1
            n_index(k) % buff(n_count(k)) = Mesh % Parallelinfo % &
                 GlobalDOFs(i)
          END DO
        END DO
      END IF

      DO i=1,ParEnv % PEs
        IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
           CALL MPI_BSEND( n_count(i), 1, MPI_INTEGER, i-1, &
                800, MPI_COMM_WORLD, ierr )
           IF ( n_count(i)>0 ) &
             CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                 801, MPI_COMM_WORLD, ierr )
        END IF
      END DO

      DO i=1,ParEnv % PEs
        IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff)

        IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
           CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                800, MPI_COMM_WORLD, status, ierr )
           IF ( n>0 ) THEN
             ALLOCATE( gbuff(n) )
             proc = status(MPI_SOURCE)
             CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                 801, MPI_COMM_WORLD, status, ierr )

             DO j=1,n
               k = SearchNodeL( Mesh % ParallelInfo, gbuff(j), Mesh % NumberOfNodes )
               IF ( k>0 ) THEN
                 IF ( BoundaryReorder(k)<= 0 ) THEN
                   NumberOfBoundaryNodes = NumberOfBoundaryNodes + 1
                   BoundaryReorder(k) = NumberOfBoundaryNodes
                 END IF
               END IF
             END DO
             DEALLOCATE(gbuff)
           END IF
        END IF
      END DO
      DEALLOCATE( n_index, n_count )
    END IF

!------------------------------------------------------------------------------

    IF ( NumberOfBoundaryNodes == 0 ) THEN
!     DEALLOCATE( BoundaryReorder )
!     NULLIFY( BoundaryReorder, BoundaryNormals,BoundaryTangent1, &
!                        BoundaryTangent2)
    ELSE
      IF ( ASSOCIATED(BoundaryNormals) ) THEN
        DEALLOCATE( BoundaryNormals, BoundaryTangent1, &
                    BoundaryTangent2, NTelement, NTzeroing_done)
      END IF

      ALLOCATE( NTelement(NumberOfBoundaryNodes,3) )
      ALLOCATE( NTzeroing_done(NumberOfBoundaryNodes,3) )
      ALLOCATE( BoundaryNormals(NumberOfBoundaryNodes,3)  )
      ALLOCATE( BoundaryTangent1(NumberOfBoundaryNodes,3) )
      ALLOCATE( BoundaryTangent2(NumberOfBoundaryNodes,3) )

      BoundaryNormals  = 0.0d0
      BoundaryTangent1 = 0.0d0
      BoundaryTangent2 = 0.0d0
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE CheckNormalTangentialBoundary
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Average boundary normals for nodes. The average boundary normals
!> may be beneficial as they provide more continuous definition of normal
!> over curved boundaries. 
!------------------------------------------------------------------------------
   SUBROUTINE AverageBoundaryNormals( Model, VariableName,    &
     NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals, &
       BoundaryTangent1, BoundaryTangent2, dim )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model

    INTEGER, POINTER :: BoundaryReorder(:)
    INTEGER :: NumberOfBoundaryNodes,DIM

    REAL(KIND=dp), POINTER :: BoundaryNormals(:,:),BoundaryTangent1(:,:), &
                       BoundaryTangent2(:,:)

    CHARACTER(LEN=*) :: VariableName
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: i,j,k,l,m,n,t, iBC, ierr, proc
    LOGICAL :: GotIt, Found, PeriodicNormals, Conditional
    REAL(KIND=dp) :: s,Bu,Bv,Nrm(3),Basis(32),DetJ
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Matrix_t), POINTER :: Projector
    REAL(KIND=dp), ALLOCATABLE :: Condition(:)

    TYPE(Variable_t), POINTER :: NrmVar, Tan1Var, Tan2Var

    LOGICAL, ALLOCATABLE :: Done(:), NtMasterBC(:), NtSlaveBC(:)
  
    REAL(KIND=dp), POINTER :: SetNormal(:,:), Rot(:,:)

    REAL(KIND=dp), TARGET :: x(Model % MaxElementNodes)
    REAL(KIND=dp), TARGET :: y(Model % MaxElementNodes)
    REAL(KIND=dp), TARGET :: z(Model % MaxElementNodes)

    TYPE buff_t
      INTEGER, ALLOCATABLE :: buff(:)
      REAL(KIND=dp), ALLOCATABLE :: normals(:)
    END TYPE buff_t
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER, POINTER :: nlist(:)
    TYPE(Buff_t), ALLOCATABLE :: n_index(:)
    REAL(KIND=dp), ALLOCATABLE :: nbuff(:)
    INTEGER, ALLOCATABLE :: n_count(:), gbuff(:), n_comp(:)

    LOGICAL :: MassConsistent, LhsSystem
    LOGICAL, ALLOCATABLE :: LhsTangent(:),RhsTangent(:)
    INTEGER :: LhsConflicts

    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------

    ElementNodes % x => x
    ElementNodes % y => y
    ElementNodes % z => z

    Mesh => Model % Mesh
    NrmVar => VariableGet( Mesh % Variables, 'Normals' )

    !dim = CoordinateSystemDimension() 


    IF ( ASSOCIATED(NrmVar) ) THEN

      IF ( NumberOfBoundaryNodes >0 ) THEN
        BoundaryNormals = 0._dp
        DO i=1,Model % NumberOfNodes
           k = BoundaryReorder(i)
           IF (k>0 ) THEN
             DO l=1,NrmVar % DOFs
                BoundaryNormals(k,l) = NrmVar % Values( NrmVar % DOFs* &
                             (NrmVar % Perm(i)-1)+l)
             END DO
           END IF
         END DO
      END IF

    ELSE

!------------------------------------------------------------------------------
!   Compute sum of elementwise normals for nodes on boundaries
!------------------------------------------------------------------------------
      ALLOCATE( n_comp(Model % NumberOfNodes) )
      n_comp = 0

      IF ( NumberOfBoundaryNodes>0 ) THEN
        BoundaryNormals = 0._dp

        DO t=Model % NumberOfBulkElements + 1, Model % NumberOfBulkElements + &
                      Model % NumberOfBoundaryElements
          Element => Model % Elements(t)
          IF ( Element % TYPE  % ElementCode < 200 ) CYCLE

          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes

          ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes)

          ALLOCATE(Condition(n))

          DO i=1,Model % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
              IF ( ListGetLogical( Model % BCs(i) % Values,VariableName, gotIt) ) THEN
                Found = ListGetLogical( Model % BCs(i) % Values, &
                    TRIM(VariableName) // ' Rotate',gotIt)
                IF ( Found .OR. .NOT. Gotit ) THEN
                  MassConsistent=ListGetLogical(Model % BCs(i) % Values, &
                          'Mass Consistent Normals',gotIt)

                  Condition(1:n) = ListGetReal( Model % BCs(i) % Values, &
                       TRIM(VariableName) // ' Condition', n, NodeIndexes, Conditional )

                  DO j=1,n
                    IF ( Conditional .AND. Condition(j) < 0._dp ) CYCLE

                    k = BoundaryReorder( NodeIndexes(j) )
                    IF (k>0) THEN
                      nrm = 0._dp
                      IF (MassConsistent) THEN
                        CALL IntegMassConsistent(j,n,nrm)
                      ELSE
                        Bu = Element % TYPE % NodeU(j)
                        Bv = Element % TYPE % NodeV(j)
                        nrm = NormalVector(Element,ElementNodes,Bu,Bv,.TRUE.)
                      END IF
                      n_comp(NodeIndexes(j)) = 1
                      BoundaryNormals(k,:) = BoundaryNormals(k,:) + nrm
                    END IF
                  END DO
                END IF
              END IF
            END IF
          END DO
          DEALLOCATE(Condition)
        END DO

        DO iBC=1,Model % NumberOfBCs
          Projector => Model % BCs(iBC) % PMatrix
          IF ( .NOT. ASSOCIATED( Projector ) ) CYCLE

          !
          ! TODO: consistent normals, if rotations given:
          ! ---------------------------------------------
          Rot => ListGetConstRealArray(Model % BCs(iBC) % Values, &
                  'Periodic BC Rotate', Found )
          IF ( Found .AND. ASSOCIATED(Rot) ) THEN
            IF ( ANY(Rot/=0) ) THEN
              ALLOCATE( Done(SIZE(BoundaryNormals,1)) )
              Done=.FALSE.
              DO i=1,Projector % NumberOfRows
                 k = BoundaryReorder(Projector % InvPerm(i))
                 IF ( k <= 0 ) CYCLE
                 DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                   IF ( Projector % Cols(l) <= 0 ) CYCLE
                   m = BoundaryReorder(Projector % Cols(l))
                   IF ( m>0 ) THEN
                     IF ( .NOT.Done(m) ) THEN
                       Done(m) = .TRUE.
                       BoundaryNormals(m,:) = -BoundaryNormals(m,:)
                     END IF
                   END IF
                 END DO
              END DO
              DEALLOCATE(Done)
              CYCLE
            END IF
          END IF

          DO i=1,Projector % NumberOfRows
            k = BoundaryReorder(Projector % InvPerm(i))
            IF ( k <= 0 ) CYCLE
            DO l=Projector % Rows(i),Projector % Rows(i+1)-1
              IF ( Projector % Cols(l) <= 0 ) CYCLE
              m = BoundaryReorder(Projector % Cols(l))
              IF ( m>0 ) BoundaryNormals(m,:) = 0._dp
            END DO
          END DO
        END DO

        DO iBC=1,Model % NumberOfBCs
           Projector => Model % BCs(iBC) % PMatrix
           IF ( .NOT. ASSOCIATED( Projector ) ) CYCLE

           !
           ! TODO: consistent normals, if rotations given:
           ! ---------------------------------------------
           Rot => ListGetConstRealArray(Model % BCs(iBC) % Values, &
                   'Periodic BC Rotate', Found )
           IF ( Found .AND. ASSOCIATED(Rot) ) THEN
             IF ( ANY(Rot/=0) ) CYCLE
           END IF

           DO i=1,Projector % NumberOfRows
              k = BoundaryReorder(Projector % InvPerm(i))
              IF ( k <= 0 ) CYCLE
              DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                IF ( Projector % Cols(l) <= 0 ) CYCLE
                m = BoundaryReorder(Projector % Cols(l))
                IF ( m > 0 ) &
                   BoundaryNormals(m,:) = BoundaryNormals(m,:) + &
                     Projector % Values(l) * BoundaryNormals(k,:)
              END DO
           END DO
        END DO
      END IF

      IF (ParEnv % PEs>1 ) THEN
        ALLOCATE( n_count(ParEnv% PEs),n_index(ParEnv % PEs) )
        n_count = 0

        IF ( NumberOfBoundaryNodes>0 ) THEN
          DO i=1,Mesh % NumberOfNodes
            IF (BoundaryReorder(i)<=0 .OR. n_comp(i)<=0 ) CYCLE
            IF (.NOT.Mesh % ParallelInfo % INTERFACE(i) ) CYCLE
  
            nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
            DO j=1,SIZE(nlist)
              k = nlist(j)+1
              IF ( k-1 == ParEnv % myPE ) CYCLE
              n_count(k) = n_count(k)+1
            END DO
          END DO
          DO i=1,ParEnv % PEs
            IF ( n_count(i)>0 ) &
                ALLOCATE( n_index(i) % buff(n_count(i)), &
                        n_index(i) % normals(3*n_count(i)) )
          END DO

          n_count = 0
          DO i=1,Model % NumberOfNodes
            IF (BoundaryReorder(i)<=0 .OR. n_comp(i)<=0 ) CYCLE
            IF (.NOT.Mesh % ParallelInfo % INTERFACE(i) ) CYCLE

            nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
            DO j=1,SIZE(nlist)
              k = nlist(j)+1
              IF ( k-1 == ParEnv % myPE ) CYCLE
              n_count(k) = n_count(k)+1
              n_index(k) % buff(n_count(k)) = Mesh % Parallelinfo % &
                 GlobalDOFs(i)
              l = BoundaryReorder(i)
              n_index(k) % normals(3*n_count(k)-2)=BoundaryNormals(l,1)
              n_index(k) % normals(3*n_count(k)-1)=BoundaryNormals(l,2)
              n_index(k) % normals(3*n_count(k)-0)=BoundaryNormals(l,3)
            END DO
          END DO
        END IF

        DO i=1,ParEnv % PEs
          IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
            CALL MPI_BSEND( n_count(i), 1, MPI_INTEGER, i-1, &
                900, MPI_COMM_WORLD, ierr )
            IF ( n_count(i)>0 ) THEN
              CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                  901, MPI_COMM_WORLD, ierr )
              CALL MPI_BSEND( n_index(i) % normals, 3*n_count(i), MPI_DOUBLE_PRECISION, &
                    i-1,  902, MPI_COMM_WORLD, ierr )
            END IF
          END IF
        END DO
        DO i=1,ParEnv % PEs
          IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff, n_index(i) % Normals)

          IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
             CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                    900, MPI_COMM_WORLD, status, ierr )
             IF ( n>0 ) THEN
               proc = status(MPI_SOURCE)
               ALLOCATE( gbuff(n), nbuff(3*n) )
               CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                   901, MPI_COMM_WORLD, status, ierr )

               CALL MPI_RECV( nbuff, 3*n, MPI_DOUBLE_PRECISION, proc, &
                    902, MPI_COMM_WORLD, status, ierr )

               DO j=1,n
                 k = SearchNodeL( Mesh % ParallelInfo, gbuff(j), Mesh % NumberOfNodes )
                 IF ( k>0 ) THEN
                   n_comp(k) = n_comp(k)+1
                   l = BoundaryReorder(k)
                   IF ( l>0 ) THEN
                     BoundaryNormals(l,1)=BoundaryNormals(l,1)+nbuff(3*j-2)
                     BoundaryNormals(l,2)=BoundaryNormals(l,2)+nbuff(3*j-1)
                     BoundaryNormals(l,3)=BoundaryNormals(l,3)+nbuff(3*j-0)
                   END IF
                 END IF
               END DO
               DEALLOCATE(gbuff, nbuff)
             END IF
          END IF
        END DO
        DEALLOCATE( n_index, n_count )
      END IF

      DEALLOCATE(n_comp)
    END IF

!------------------------------------------------------------------------------
!   normalize 
!------------------------------------------------------------------------------
    IF ( NumberOfBoundaryNodes>0 ) THEN

      LhsSystem = ListGetLogical(Model % Simulation,'Use Lhs System',Found) 
      IF(.NOT. Found ) LhsSystem = ( dim == 3 )

      IF( LhsSystem ) THEN
        ALLOCATE( NtMasterBC( Model % NumberOfBCs ), NtSlaveBC( Model % NumberOfBCs ) )
        NtMasterBC = .FALSE.; NtSlaveBC = .FALSE.

        DO i = 1, Model % NumberOfBcs
          IF( .NOT. ListCheckPrefix( Model % BCs(i) % Values,'Normal-Tangential') ) CYCLE
          
          j = ListGetInteger( Model % BCs(i) % Values,'Mortar BC',Found )
          IF( .NOT. Found ) THEN
            j = ListGetInteger( Model % BCs(i) % Values,'Contact BC',Found )
          END IF
          IF( j == 0 .OR. j > Model % NumberOfBCs ) CYCLE

          NtSlaveBC( i ) = .TRUE.
          NtMasterBC( j ) = .TRUE.
        END DO
        LhsSystem = ANY( NtMasterBC )
      END IF

      IF( LhsSystem ) THEN
        DO i = 1, Model % NumberOfBcs
          IF( NtSlaveBC( i ) .AND. NtMasterBC( i ) ) THEN
            CALL Warn('AverageBoundaryNormals','BC '//TRIM(I2S(i))//' is both N-T master and slave!')
          END IF
        END DO

        ALLOCATE( LhsTangent( Model % NumberOfNodes ) )
        LhsTangent = .FALSE.

        ALLOCATE( RhsTangent( Model % NumberOfNodes ) )
        RhsTangent = .FALSE. 

        DO t=Model % NumberOfBulkElements + 1, Model % NumberOfBulkElements + &
            Model % NumberOfBoundaryElements
          Element => Model % Elements(t)
          IF ( Element % TYPE  % ElementCode < 200 ) CYCLE
          
          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          
          DO i=1,Model % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
              IF( NtMasterBC(i) ) LhsTangent( NodeIndexes ) = .TRUE.
              IF( NtSlaveBC(i) ) RhsTangent( NodeIndexes ) = .TRUE.
              EXIT
            END IF
          END DO
        END DO

        LhsConflicts = COUNT( LhsTangent .AND. RhsTangent )
        IF( LhsConflicts > 0 ) THEN
          CALL Warn('AverageBoundaryNormals',&
              'There are '//TRIM(I2S(LhsConflicts))//' nodes that could be both rhs and lhs!')
        END IF
      END IF


      DO i=1,Model % NumberOfNodes
        k = BoundaryReorder(i) 
        IF ( k > 0 ) THEN
          s = SQRT( SUM( BoundaryNormals(k,:)**2 ) )
          IF ( s /= 0.0d0 ) &
            BoundaryNormals(k,:) = BoundaryNormals(k,:) / s
          IF ( dim > 2 ) THEN
            CALL TangentDirections( BoundaryNormals(k,:),  &
                BoundaryTangent1(k,:), BoundaryTangent2(k,:) )
            IF( LhsSystem ) THEN
              IF( LhsTangent(i) ) THEN
                BoundaryTangent2(k,:) = -BoundaryTangent2(k,:)
              END IF
            END IF
          END IF
        END IF
      END DO
      
      IF( ListGetLogical( Model % Simulation,'Save Averaged Normals',Found ) ) THEN
        CALL Info('AverageBoundaryNormals','Saving averaged boundary normals to variable: Averaged Normals')
        NrmVar => VariableGet( Mesh % Variables, 'Averaged Normals' )
        
        IF(.NOT. ASSOCIATED( NrmVar ) ) THEN
          CALL VariableAddVector( Mesh % Variables, Mesh, Model % Solver,'Averaged Normals',3,&
              Perm = BoundaryReorder )
          NrmVar => VariableGet( Mesh % Variables, 'Averaged Normals' )
        END IF
            
        DO i=1,Model % NumberOfNodes
          k = BoundaryReorder(i)
          IF (k>0 ) THEN
            DO l=1,NrmVar % DOFs
              NrmVar % Values( NrmVar % DOFs* &
                  (NrmVar % Perm(i)-1)+l)  = BoundaryNormals(k,l)
            END DO
          END IF
        END DO

        IF( dim > 2 .AND. ListGetLogical( Model % Simulation,'Save Averaged Tangents',Found ) ) THEN
          Tan1Var => VariableGet( Mesh % Variables, 'Averaged First Tangent' )
          Tan2Var => VariableGet( Mesh % Variables, 'Averaged Second Tangent' )

          IF(.NOT. ASSOCIATED( Tan1Var ) ) THEN
            CALL VariableAddVector( Mesh % Variables, Mesh, Model % Solver,&
                'Averaged First Tangent',3, Perm = BoundaryReorder )
            Tan1Var => VariableGet( Mesh % Variables, 'Averaged First Tangent' )
            CALL VariableAddVector( Mesh % Variables, Mesh, Model % Solver,&
                'Averaged Second Tangent',3, Perm = BoundaryReorder )
            Tan2Var => VariableGet( Mesh % Variables, 'Averaged Second Tangent' )
          END IF
          
          DO i=1,Model % NumberOfNodes
            k = BoundaryReorder(i)
            IF (k>0 ) THEN
              DO l=1,Tan1Var % DOFs
                Tan1Var % Values( Tan1Var % DOFs* &
                    (Tan1Var % Perm(i)-1)+l)  = BoundaryTangent1(k,l)
                Tan2Var % Values( Tan2Var % DOFs* &
                    (Tan2Var % Perm(i)-1)+l)  = BoundaryTangent2(k,l)
              END DO
            END IF
          END DO
        END IF
      END IF
    END IF



 CONTAINS

    SUBROUTINE IntegMassConsistent(j,n,nrm)
      INTEGER :: t,j,n
      LOGICAL :: stat
      REAL(KIND=dp) :: detJ,Basis(n),nrm(:),lnrm(3)

      TYPE(GaussIntegrationPoints_t) :: IP

      !----------------------
      IP = GaussPoints(Element)
      DO t=1,IP % n
        stat = ElementInfo(Element, ElementNodes, IP % U(t), &
               IP % v(t), IP % W(t), detJ, Basis)

        lnrm = NormalVector(Element,ElementNodes, &
              IP % U(t),IP % v(t),.TRUE.)

        nrm = nrm + IP % s(t) * lnrm * detJ * Basis(j)
      END DO
    END SUBROUTINE IntegMassConsistent

!------------------------------------------------------------------------------
  END SUBROUTINE AverageBoundaryNormals
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Search an element QueriedNode from an ordered set Nodes and return
!> Index to Nodes structure. Return value -1 means QueriedNode was
!> not found.
!------------------------------------------------------------------------------
FUNCTION SearchNodeL( ParallelInfo, QueriedNode,n ) RESULT(Indx)

  USE Types
  IMPLICIT NONE

  TYPE (ParallelInfo_t) :: ParallelInfo
  INTEGER :: QueriedNode, Indx,n

  ! Local variables

  INTEGER :: Lower, Upper, Lou, i

!------------------------------------------------------------------------------

  Indx = -1
  Upper = n
  Lower = 1

  ! Handle the special case

  IF ( Upper == 0 ) RETURN

10 CONTINUE
  IF ( ParallelInfo % GlobalDOFs(Lower) == QueriedNode ) THEN
     Indx = Lower
     RETURN
  ELSE IF ( ParallelInfo % GlobalDOFs(Upper) == QueriedNode ) THEN
     Indx = Upper
     RETURN
  END IF

  IF ( (Upper - Lower) > 1 ) THEN
     Lou = ISHFT((Upper + Lower), -1)
     IF ( ParallelInfo % GlobalDOFs(Lou) < QueriedNode ) THEN
        Lower = Lou
        GOTO 10
     ELSE
        Upper = Lou
        GOTO 10
     END IF
  END IF

  RETURN
!------------------------------------------------------------------------------
END FUNCTION SearchNodeL
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Initialize solver for next timestep.
!------------------------------------------------------------------------------
  SUBROUTINE InitializeTimestep( Solver )
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver  !< Solver to be initialized.
!------------------------------------------------------------------------------
     CHARACTER(LEN=MAX_NAME_LEN) :: Method
     LOGICAL :: GotIt
     INTEGER :: i, Order,ndofs
     REAL(KIND=dp), POINTER :: Work(:), SaveValues(:)
     TYPE(Matrix_t), POINTER :: A

!------------------------------------------------------------------------------
     Solver % DoneTime = Solver % DoneTime + 1
!------------------------------------------------------------------------------

     IF ( .NOT. ASSOCIATED( Solver % Matrix ) .OR. &
          .NOT. ASSOCIATED( Solver % Variable % Values ) ) RETURN

     IF ( Solver % TimeOrder <= 0 ) RETURN
!------------------------------------------------------------------------------

     Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )
     IF ( Method == 'none' ) RETURN
    
     IF ( .NOT.GotIt ) THEN

        Solver % Beta = ListGetConstReal( Solver % Values, 'Newmark Beta', GotIt )
        IF ( .NOT. GotIt ) THEN
           Solver % Beta = ListGetConstReal( CurrentModel % Simulation, 'Newmark Beta', GotIt )
       END IF

       IF ( .NOT.GotIt ) THEN
         CALL Warn( 'InitializeTimestep', &
               'Timestepping method defaulted to IMPLICIT EULER' )

         Solver % Beta = 1.0D0
         Method = 'implicit euler'
       END IF

     ELSE

       Solver % Beta = 1._dp
       SELECT CASE( Method )
         CASE('implicit euler')
           Solver % Beta = 1.0d0

         CASE('explicit euler')
           Solver % Beta = 0.0d0

         CASE('runge-kutta')
           Solver % Beta = 0.0d0

         CASE('crank-nicolson')
           Solver % Beta = 0.5d0

         CASE('fs')
           Solver % Beta = 0.5d0

         CASE('newmark')
           Solver % Beta = ListGetConstReal( Solver % Values, 'Newmark Beta', GotIt )
           IF ( .NOT. GotIt ) THEN
              Solver % Beta = ListGetConstReal( CurrentModel % Simulation, &
                              'Newmark Beta', GotIt )
           END IF

           IF ( Solver % Beta<0 .OR. Solver % Beta>1 ) THEN
             WRITE( Message, * ) 'Invalid value of Beta ', Solver % Beta
             CALL Warn( 'InitializeTimestep', Message )
           END IF

         CASE('bdf')
           IF ( Solver % Order < 1 .OR. Solver % Order > 5  ) THEN
             WRITE( Message, * ) 'Invalid order BDF ',  Solver % Order
             CALL Fatal( 'InitializeTimestep', Message )
           END IF

         CASE DEFAULT 
           WRITE( Message, * ) 'Unknown timestepping method: ',Method
           CALL Fatal( 'InitializeTimestep', Message )
       END SELECT

     END IF

     ndofs = Solver % Matrix % NumberOfRows

     IF ( Method /= 'bdf' .OR. Solver % TimeOrder > 1 ) THEN
       IF ( Solver % DoneTime == 1 .AND. Solver % Beta /= 0.0d0 ) THEN
         Solver % Beta = 1.0d0
       END IF
 
       SELECT CASE( Solver % TimeOrder )
         CASE(1)
           Order = MIN(Solver % DoneTime, Solver % Order)
           DO i=Order, 2, -1
             Solver % Variable % PrevValues(:,i) = &
                   Solver % Variable % PrevValues(:,i-1)
           END DO
           Solver % Variable % PrevValues(:,1) = Solver % Variable % Values
           Solver % Matrix % Force(:,2) = Solver % Matrix % Force(:,1)

         CASE(2)
           SELECT CASE(Method)
           CASE DEFAULT
             Solver % Alpha = ListGetConstReal( Solver % Values, &
                        'Bossak Alpha', GotIt )
             IF ( .NOT. GotIt ) THEN
                 Solver % Alpha = ListGetConstReal( CurrentModel % Simulation, &
                            'Bossak Alpha', GotIt )
             END IF
             IF ( .NOT. GotIt ) Solver % Alpha = -0.05d0

             Solver % Variable % PrevValues(:,3) = &
                                 Solver % Variable % Values
             Solver % Variable % PrevValues(:,4) = &
                        Solver % Variable % PrevValues(:,1)
             Solver % Variable % PrevValues(:,5) = &
                        Solver % Variable % PrevValues(:,2)
           END SELECT
       END SELECT
     ELSE
       Order = MIN(Solver % DoneTime, Solver % Order)
       DO i=Order, 2, -1
         Solver % Variable % PrevValues(:,i) = &
               Solver % Variable % PrevValues(:,i-1)
       END DO
       Solver % Variable % PrevValues(:,1) = Solver % Variable % Values
     END IF


     IF( ListGetLogical( Solver % Values,'Nonlinear Timestepping', GotIt ) ) THEN
       IF( Solver % DoneTime > 1 ) THEN
         A => Solver % Matrix
         CALL Info('InitializeTimestep','Saving previous linear system for timestepping',Level=12)
         IF( .NOT. ASSOCIATED( A % BulkValues ) ) THEN
           CALL Fatal('InitializeTimestep','BulkValues should be associated!')
         END IF
         
         IF( .NOT. ASSOCIATED( A % BulkResidual ) ) THEN
           ALLOCATE( A % BulkResidual( SIZE( A % BulkRhs ) ) )
         END IF
         
         SaveValues => A % Values
         A % Values => A % BulkValues
         CALL MatrixVectorMultiply( A, &
             Solver % Variable % Values, A % BulkResidual )
         A % Values => SaveValues
         A % BulkResidual = A % BulkResidual - A % BulkRhs
       END IF
     END IF

!------------------------------------------------------------------------------
  END SUBROUTINE InitializeTimestep
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Update force vector AFTER ALL OTHER ASSEMBLY STEPS BUT BEFORE SETTING
!> DIRICHLET CONDITIONS. Required only for time dependent simulations..
!------------------------------------------------------------------------------
  SUBROUTINE FinishAssembly( Solver, ForceVector )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: ForceVector(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: Method, Simulation
    INTEGER :: Order
    LOGICAL :: Found
!------------------------------------------------------------------------------

    IF ( Solver % Matrix % FORMAT == MATRIX_LIST ) THEN
      CALL List_toCRSMatrix(Solver % Matrix)
    END IF

    Simulation = ListGetString( CurrentModel % Simulation, 'Simulation Type' )
    IF ( Simulation == 'transient' ) THEN
      Method = ListGetString( Solver % Values, 'Timestepping Method' )
      Order = MIN(Solver % DoneTime, Solver % Order)

      IF ( Order <= 0 .OR. Solver % TimeOrder /= 1 .OR. Method=='bdf' ) RETURN

      IF ( Solver % Beta /= 0.0d0 ) THEN
        ForceVector = ForceVector + ( Solver % Beta - 1 ) * &
            Solver % Matrix % Force(:,1) + &
                ( 1 - Solver % Beta ) * Solver % Matrix % Force(:,2)
      END IF
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE FinishAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE InvalidateVariable( TopMesh,PrimaryMesh,Name )
!------------------------------------------------------------------------------
    CHARACTER(LEN=*) :: Name
    TYPE(Mesh_t),  POINTER :: TopMesh,PrimaryMesh
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: tmpname
    INTEGER :: i
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var,Var1, PrimVar
!------------------------------------------------------------------------------
    Mesh => TopMesh

    PrimVar => VariableGet( PrimaryMesh % Variables, Name, ThisOnly=.TRUE.)
    IF ( .NOT.ASSOCIATED( PrimVar) ) RETURN

    DO WHILE( ASSOCIATED(Mesh) )
      IF ( .NOT.ASSOCIATED( PrimaryMesh, Mesh) ) THEN
        Var => VariableGet( Mesh % Variables, Name, ThisOnly=.TRUE.)
        IF ( ASSOCIATED( Var ) ) THEN
          Var % Valid = .FALSE.
          Var % PrimaryMesh => PrimaryMesh
        END IF

        IF ( PrimVar % DOFs > 1 ) THEN
          IF ( .FALSE. .AND. PrimVar % Name == 'flow solution' ) THEN
            Var1 => VariableGet( Mesh % Variables, 'Velocity 1', .TRUE.)
            IF ( ASSOCIATED( Var1 ) ) THEN
               Var1 % Valid = .FALSE.
               Var1 % PrimaryMesh => PrimaryMesh
            END IF
            Var1 => VariableGet( Mesh % Variables, 'Velocity 2', .TRUE.)
            IF ( ASSOCIATED( Var1 ) ) THEN
               Var1 % Valid = .FALSE.
               Var1 % PrimaryMesh => PrimaryMesh
            END IF
            Var1 => VariableGet( Mesh % Variables, 'Velocity 3', .TRUE.)
            IF ( ASSOCIATED( Var1 ) ) THEN
               Var1 % Valid = .FALSE.
               Var1 % PrimaryMesh => PrimaryMesh
            END IF
            Var1 => VariableGet( Mesh % Variables, 'Pressure', .TRUE.)
            IF ( ASSOCIATED( Var1 ) ) THEN
               Var1 % Valid = .FALSE.
               Var1 % PrimaryMesh => PrimaryMesh
            END IF
            Var1 => VariableGet( Mesh % Variables, 'Surface', .TRUE.)
            IF ( ASSOCIATED( Var1 ) ) THEN
               Var1 % Valid = .FALSE.
               Var1 % PrimaryMesh => PrimaryMesh
            END IF
          ELSE
            DO i=1,PrimVar % DOFs
              tmpname = ComponentName( Name, i )
              Var1 => VariableGet( Mesh % Variables, tmpname, .TRUE. )
              IF ( ASSOCIATED( Var1 ) ) THEN
                 Var1 % Valid = .FALSE.
                 Var1 % PrimaryMesh => PrimaryMesh
              END IF
            END DO
          END IF
        END IF
      END IF
      Mesh => Mesh % Next
    END DO 

    PrimVar % ValuesChanged = .TRUE.
    IF ( PrimVar % DOFs > 1 ) THEN
      IF ( .FALSE. .AND. PrimVar % Name == 'flow solution' ) THEN
        Var => VariableGet( PrimaryMesh % Variables, 'Surface', .TRUE.)
        IF ( ASSOCIATED(Var) ) Var % ValuesChanged = .TRUE.
        Var => VariableGet( PrimaryMesh % Variables, 'Pressure', .TRUE.)
        IF ( ASSOCIATED(Var) ) Var % ValuesChanged = .TRUE.
        Var => VariableGet( PrimaryMesh % Variables, 'Velocity 1', .TRUE.)
        IF ( ASSOCIATED(Var) ) Var % ValuesChanged = .TRUE.
        Var => VariableGet( PrimaryMesh % Variables, 'Velocity 2', .TRUE.)
        IF ( ASSOCIATED(Var) ) Var % ValuesChanged = .TRUE.
        Var => VariableGet( PrimaryMesh % Variables, 'Velocity 3', .TRUE.)
        IF ( ASSOCIATED(Var) ) Var % ValuesChanged = .TRUE.
      ELSE
        DO i=1,PrimVar % DOFs
          tmpname = ComponentName( Name, i )
          Var => VariableGet( PrimaryMesh % Variables, tmpname, .TRUE. )
          IF ( ASSOCIATED(Var) ) Var % ValuesChanged = .TRUE.
        END DO
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE InvalidateVariable
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Rotate a vector to normal-tangential coordinate system.
!------------------------------------------------------------------------------
  SUBROUTINE RotateNTSystem( Vec, NodeNumber )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Vec(:)
     INTEGER :: NodeNumber
!------------------------------------------------------------------------------
     INTEGER :: i,j,k, dim
     REAL(KIND=dp) :: Bu,Bv,Bw,RM(3,3)
!------------------------------------------------------------------------------

     IF ( NormalTangentialNOFNodes <= 0 ) RETURN

     dim = CoordinateSystemDimension()

     k = BoundaryReorder(NodeNumber)
     IF ( k <= 0 ) RETURN

     IF ( dim < 3 ) THEN
       Bu = Vec(1)
       Bv = Vec(2)
       Vec(1) =  BoundaryNormals(k,1)*Bu + BoundaryNormals(k,2)*Bv
       Vec(2) = -BoundaryNormals(k,2)*Bu + BoundaryNormals(k,1)*Bv
     ELSE
       Bu = Vec(1)
       Bv = Vec(2)
       Bw = Vec(3)

       RM(:,1) = BoundaryNormals(k,:)
       RM(:,2) = BoundaryTangent1(k,:)
       RM(:,3) = BoundaryTangent2(k,:)

       Vec(1) = RM(1,1)*Bu + RM(2,1)*Bv + RM(3,1)*Bw
       Vec(2) = RM(1,2)*Bu + RM(2,2)*Bv + RM(3,2)*Bw
       Vec(3) = RM(1,3)*Bu + RM(2,3)*Bv + RM(3,3)*Bw
     END IF
!------------------------------------------------------------------------------
  END SUBROUTINE RotateNTSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
!> Rotate all components of a solution vector to normal-tangential coordinate system
!------------------------------------------------------------------------------------
  SUBROUTINE RotateNTSystemAll( Solution, Perm, NDOFs )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Solution(:)
    INTEGER :: Perm(:), NDOFs
!------------------------------------------------------------------------------
    INTEGER :: i,j,k, dim
    REAL(KIND=dp) :: Bu,Bv,Bw,RM(3,3)
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()

    IF ( NormalTangentialNOFNodes<=0.OR.ndofs<dim ) RETURN

    DO i=1,SIZE(BoundaryReorder)
       k = BoundaryReorder(i)
       IF ( k <= 0 ) CYCLE
       j = Perm(i)
       IF ( j <= 0 ) CYCLE

       IF ( dim < 3 ) THEN
          Bu = Solution(NDOFs*(j-1)+1)
          Bv = Solution(NDOFs*(j-1)+2)

          Solution(NDOFs*(j-1)+1) = BoundaryNormals(k,1)*Bu + BoundaryNormals(k,2)*Bv
          Solution(NDOFs*(j-1)+2) = -BoundaryNormals(k,2)*Bu + BoundaryNormals(k,1)*Bv

       ELSE
          Bu = Solution(NDOFs*(j-1)+1)
          Bv = Solution(NDOFs*(j-1)+2)
          Bw = Solution(NDOFs*(j-1)+3)
 
          RM(:,1) = BoundaryNormals(k,:)
          RM(:,2) = BoundaryTangent1(k,:)
          RM(:,3) = BoundaryTangent2(k,:)

          Solution(NDOFs*(j-1)+1) = RM(1,1)*Bu + RM(2,1)*Bv + RM(3,1)*Bw
          Solution(NDOFs*(j-1)+2) = RM(1,2)*Bu + RM(2,2)*Bv + RM(3,2)*Bw
          Solution(NDOFs*(j-1)+3) = RM(1,3)*Bu + RM(2,3)*Bv + RM(3,3)*Bw
       END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE RotateNTSystemAll
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Backrotate a solution from normal-tangential coordinate system to cartesian one.
!------------------------------------------------------------------------------
  SUBROUTINE BackRotateNTSystem( Solution, Perm, NDOFs )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Solution(:)
     INTEGER :: Perm(:), NDOFs
!------------------------------------------------------------------------------
     INTEGER :: i,j,k, dim
     REAL(KIND=dp) :: Bu,Bv,Bw,RM(3,3)
!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()

     IF ( NormalTangentialNOFNodes<=0.OR.ndofs<dim ) RETURN

     DO i=1,SIZE(BoundaryReorder)
       k = BoundaryReorder(i)
       IF ( k <= 0 ) CYCLE
       j = Perm(i)
       IF ( j <= 0 ) CYCLE

       IF ( dim < 3 ) THEN
         Bu = Solution(NDOFs*(j-1)+1)
         Bv = Solution(NDOFs*(j-1)+2)

         Solution(NDOFs*(j-1)+1) = BoundaryNormals(k,1) * Bu - &
                         BoundaryNormals(k,2) * Bv

         Solution(NDOFs*(j-1)+2) = BoundaryNormals(k,2) * Bu + &
                         BoundaryNormals(k,1) * Bv
       ELSE
         Bu = Solution(NDOFs*(j-1)+1)
         Bv = Solution(NDOFs*(j-1)+2)
         Bw = Solution(NDOFs*(j-1)+3)

         RM(1,:) = BoundaryNormals(k,:)
         RM(2,:) = BoundaryTangent1(k,:)
         RM(3,:) = BoundaryTangent2(k,:)

         Solution(NDOFs*(j-1)+1) = RM(1,1)*Bu + RM(2,1)*Bv + RM(3,1)*Bw
         Solution(NDOFs*(j-1)+2) = RM(1,2)*Bu + RM(2,2)*Bv + RM(3,2)*Bw
         Solution(NDOFs*(j-1)+3) = RM(1,3)*Bu + RM(2,3)*Bv + RM(3,3)*Bw
       END IF
     END DO 
!------------------------------------------------------------------------------
  END SUBROUTINE BackRotateNTSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetSolutionRotation(A,n) RESULT(rotated)
!------------------------------------------------------------------------------
    INTEGER :: n
    LOGICAL :: rotated
    REAL(KIND=dp) :: A(3,3)
!------------------------------------------------------------------------------
    INTEGER :: k,dim
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()

    Rotated=.FALSE.

    A=0._dp
    A(1,1)=1._dp
    A(2,2)=1._dp
    A(3,3)=1._dp
    IF (NormalTangentialNOFNodes<=0) RETURN

    k = BoundaryReorder(n)
    IF (k>0) THEN
      Rotated = .TRUE.
      IF (dim==2) THEN
        A(1,1)= BoundaryNormals(k,1)
        A(1,2)=-BoundaryNormals(k,2)
        A(2,1)= BoundaryNormals(k,2)
        A(2,2)= BoundaryNormals(k,1)
      ELSE
        A(:,1)=BoundaryNormals(k,:)
        A(:,2)=BoundaryTangent1(k,:)
        A(:,3)=BoundaryTangent2(k,:)
      END IF
    END IF
!------------------------------------------------------------------------------
  END FUNCTION GetSolutionRotation
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Computes the norm related to a solution vector of the Solver.
!------------------------------------------------------------------------------
  FUNCTION ComputeNorm(Solver, nin, values) RESULT (Norm)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Solver_t), TARGET :: Solver
    INTEGER :: nin
    REAL(KIND=dp), TARGET, OPTIONAL :: values(:)
    
    INTEGER :: NormDim, NormDofs, Dofs,i,j,k,n,totn,PermStart
    INTEGER, POINTER :: NormComponents(:)
    INTEGER, ALLOCATABLE :: iPerm(:)
    REAL(KIND=dp) :: Norm, nscale, val
    LOGICAL :: Stat, ComponentsAllocated, ConsistentNorm
    REAL(KIND=dp), POINTER :: x(:)
    REAL(KIND=dp), ALLOCATABLE, TARGET :: y(:)

    CALL Info('ComputeNorm','Computing norm of solution',Level=10)

    IF(PRESENT(values)) THEN
      x => values
    ELSE
      x => Solver % Variable % Values
    END IF
    
    NormDim = ListGetInteger(Solver % Values,'Nonlinear System Norm Degree',Stat)
    IF(.NOT. Stat) NormDim = 2

    Dofs = Solver % Variable % Dofs

    ComponentsAllocated = .FALSE.
    NormComponents => ListGetIntegerArray(Solver % Values,&
        'Nonlinear System Norm Components',Stat)
    IF(Stat) THEN
      NormDofs = SIZE( NormComponents ) 
    ELSE
      NormDofs = ListGetInteger(Solver % Values,'Nonlinear System Norm Dofs',Stat)
      IF(Stat) THEN
        ALLOCATE(NormComponents(NormDofs))
        ComponentsAllocated = .TRUE.
        DO i=1,NormDofs
          NormComponents(i) = i
        END DO
      ELSE
        NormDofs = Dofs        
      END IF
    END IF
 
    n = nin

    IF( ParEnv % PEs > 1 ) THEN
      ConsistentNorm = ListGetLogical(Solver % Values,'Nonlinear System Consistent Norm',Stat)
      CALL Info('ComputeNorm','Using consistent norm in parallel',Level=10)
    ELSE
      ConsistentNorm = .FALSE.
    END IF


    PermStart = ListGetInteger(Solver % Values,'Norm Permutation',Stat)
    IF ( Stat ) THEN
      ALLOCATE(iPerm(SIZE(Solver % Variable % Perm))); iPerm=0
      n = 0
      DO i=PermStart,SIZE(iPerm)
        IF ( Solver % Variable % Perm(i)>0 ) THEN
          n = n + 1
          iPerm(n) = Solver % Variable % Perm(i)
        END IF
      END DO
      ALLOCATE(y(n))
      y = x(iPerm(1:n))
      x => y
      DEALLOCATE(iPerm)
    END IF


    IF( NormDofs < Dofs ) THEN
      IF( ConsistentNorm ) THEN
        CALL Warn('ComputeNorm','Consistent norm not implemented for selective norm')
      END IF

      totn = NINT( ParallelReduction(1._dp*n) )
      nscale = NormDOFs*totn/(1._dp*DOFs)
      Norm = 0.0_dp

      SELECT CASE(NormDim)
      CASE(0)
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = MAX(Norm, MAXVAL( ABS(x(j::Dofs))) )
        END DO
        Norm = ParallelReduction(Norm,2)
      CASE(1)
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = Norm + SUM( ABS(x(j::Dofs)) )
        END DO
        Norm = ParallelReduction(Norm)/nscale
      CASE(2)
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = Norm + SUM( x(j::Dofs)**2 )
        END DO
        Norm = SQRT(ParallelReduction(Norm)/nscale)
      CASE DEFAULT
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = Norm + SUM( x(j::Dofs)**NormDim )
        END DO
        Norm = (ParallelReduction(Norm)/nscale)**(1.0d0/NormDim)
      END SELECT
    ELSE IF( ConsistentNorm ) THEN
      ! In consistent norm we have to skip the dofs not owned by the partition in order
      ! to count each dof only once. 

      Norm = 0.0_dp
      totn = 0
      DO j=1,n
        IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
            == ParEnv % MyPE ) totn = totn + 1
      END DO        

      totn = NINT( ParallelReduction(1._dp*totn) )
      nscale = 1.0_dp * totn

      SELECT CASE(NormDim)
        
      CASE(0) 
        DO j=1,n
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = MAX( Norm, ABS( val ) )
        END DO
        
      CASE(1)
        DO j=1,n
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = Norm + ABS(val)
        END DO
        
      CASE(2)          
        DO j=1,n
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          
          Norm = Norm + val**2 
        END DO
        
      CASE DEFAULT
        DO j=1,n
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = Norm + val**NormDim 
        END DO
      END SELECT

      SELECT CASE(NormDim)
      CASE(0)
        Norm = ParallelReduction(Norm,2)
      CASE(1)
        Norm = ParallelReduction(Norm) / nscale
      CASE(2)
        Norm = SQRT(ParallelReduction(Norm)/nscale)
      CASE DEFAULT
        Norm = (ParallelReduction(Norm)/nscale)**(1.0d0/NormDim)
      END SELECT
      
    ELSE
      totn = NINT( ParallelReduction(1._dp*n) )
      nscale = 1.0_dp * totn

      SELECT CASE(NormDim)
      CASE(0)
        Norm = ParallelReduction(MAXVAL(ABS(x(1:n))),2)
      CASE(1)
        Norm = ParallelReduction(SUM(ABS(x(1:n))))/nscale
      CASE(2)
        Norm = SQRT(ParallelReduction(SUM(x(1:n)**2))/nscale)
      CASE DEFAULT
        Norm = (ParallelReduction(SUM(x(1:n)**NormDim))/nscale)**(1.0d0/NormDim)
      END SELECT
    END IF

    IF( ComponentsAllocated ) THEN
      DEALLOCATE( NormComponents ) 
    END IF
!------------------------------------------------------------------------------
  END FUNCTION ComputeNorm
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
!> When a new field has been computed compare it to the previous one.
!> Different convergence measures may be used. 
!> Also performs relaxation if a non-unity relaxation factor is given.
!------------------------------------------------------------------------------
  SUBROUTINE ComputeChange(Solver,SteadyState,nsize,values,values0,Matrix,RHS)
!------------------------------------------------------------------------------
    TYPE(Solver_t), TARGET :: Solver
    LOGICAL :: SteadyState
    TYPE(Matrix_t), OPTIONAL, TARGET :: Matrix
    INTEGER, OPTIONAL :: nsize
    REAL(KIND=dp), OPTIONAL, TARGET :: values(:), values0(:), RHS(:)
!------------------------------------------------------------------------------
    INTEGER :: i, n, nn, RelaxAfter, IterNo
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp), POINTER :: b(:), x(:), r(:)
    REAL(KIND=dp), POINTER :: x0(:)
    REAL(KIND=dp) :: Norm, PrevNorm, rNorm, bNorm, Change, Relaxation, tmp(1),dt, &
        Tolerance, MaxNorm, eps, Ctarget
    CHARACTER(LEN=MAX_NAME_LEN) :: ConvergenceType
    INTEGER, TARGET  ::  Dnodes(1)
    INTEGER, POINTER :: Indexes(:)
    TYPE(Variable_t), POINTER :: iterV, VeloVar, TimestepVar, WeightVar
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, str
    LOGICAL :: Stat, ConvergenceAbsolute, Relax, RelaxBefore, DoIt, Skip, &
        SkipConstraints, ResidualMode 

    TYPE(Matrix_t), POINTER :: MMatrix
    REAL(KIND=dp), POINTER CONTIG :: Mx(:), Mb(:), Mr(:)
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: TmpXVec, TmpRVec, TmpRHSVec
    INTEGER :: ipar(1)
    TYPE(ValueList_t), POINTER :: SolverParams

    SolverParams => Solver % Values
    
  
    IF(SteadyState) THEN	
      Skip = ListGetLogical( SolverParams,'Skip Compute Steady State Change',Stat)
      IF( Skip ) RETURN

      ! No residual mode for steady state analysis
      ResidualMode = .FALSE.

      ConvergenceType = ListGetString(SolverParams,&
          'Steady State Convergence Measure',Stat)
      IF(.NOT. Stat) ConvergenceType = 'norm' 

      ConvergenceAbsolute = &
          ListGetLogical(SolverParams,'Steady State Convergence Absolute',Stat)
      IF(.NOT. Stat) ConvergenceAbsolute = &
          ListGetLogical(SolverParams,'Use Absolute Norm for Convergence',Stat)

      Relaxation = ListGetCReal( SolverParams, &
          'Steady State Relaxation Factor', Relax )
      Relax = Relax .AND. ABS(Relaxation-1.0_dp) > TINY(Relaxation)

      iterV => VariableGet( Solver % Mesh % Variables, 'coupled iter' )
      IterNo = iterV % Values(1)
      IF( Relax ) THEN
        RelaxAfter = ListGetInteger(SolverParams,'Steady State Relaxation After',Stat)
        IF( Stat .AND. RelaxAfter >= IterNo ) Relax = .FALSE.
      END IF	

      RelaxBefore = .TRUE.
      IF(Relax) THEN
        RelaxBefore = ListGetLogical( SolverParams, &
            'Steady State Relaxation Before', Stat )      
        IF (.NOT. Stat ) RelaxBefore = .TRUE.
      END IF

      ! Steady state system has never any constraints
      SkipConstraints = .FALSE.

    ELSE
      iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
      IterNo = iterV % Values(1)
      Solver % Variable % NonlinIter = iterV % Values(1)
      iterV % Values(1) = iterV % Values(1) + 1 

      Skip = ListGetLogical( SolverParams,'Skip Compute Nonlinear Change',Stat)
      IF(Skip) RETURN

      ResidualMode = ListGetLogical( SolverParams,'Linear System Residual Mode',Stat)

      ConvergenceType = ListGetString(SolverParams,&
          'Nonlinear System Convergence Measure',Stat)
      IF(.NOT. stat) ConvergenceType = 'norm' 

      ConvergenceAbsolute = &
          ListGetLogical(SolverParams,'Nonlinear System Convergence Absolute',Stat)
      IF(.NOT. Stat) ConvergenceAbsolute = &
          ListGetLogical(SolverParams,'Use Absolute Norm for Convergence',Stat)

      Relaxation = ListGetCReal( SolverParams, &
          'Nonlinear System Relaxation Factor', Relax )
      Relax = Relax .AND. (Relaxation /= 1.0d0)
      IF( Relax ) THEN
        RelaxAfter = ListGetInteger(SolverParams,'Nonlinear System Relaxation After',Stat)
        IF( Stat .AND. RelaxAfter >= Solver % Variable % NonlinIter ) Relax = .FALSE.
      END IF	

      SkipConstraints = ListGetLogical(SolverParams,&
          'Nonlinear System Convergence Without Constraints',Stat) 

      RelaxBefore = .TRUE.
      IF(Relax) THEN
        RelaxBefore = ListGetLogical( SolverParams, &
            'Nonlinear System Relaxation Before', Stat )
        IF (.NOT. Stat ) RelaxBefore = .TRUE.
      END IF
    END IF


    IF(PRESENT(values)) THEN
      x => values
    ELSE
      x => Solver % Variable % Values      
    END IF

    IF ( .NOT. ASSOCIATED(x) ) THEN
      Solver % Variable % Norm = 0.0d0 
      IF(SteadyState) THEN
        Solver % Variable % SteadyChange = 0.0d0
      ELSE
        Solver % Variable % NonlinChange = 0.0d0
      END IF
      RETURN
    END IF

    IF(PRESENT(nsize)) THEN
      n = nsize 
    ELSE 
      n = SIZE( x )
    END IF

    IF( SkipConstraints ) n = MIN( n, Solver % Matrix % NumberOfRows )

    Stat = .FALSE.
    IF(PRESENT(values0)) THEN
      x0 => values0
      Stat = .TRUE.
    ELSE IF(SteadyState) THEN
      IF( ASSOCIATED(Solver % Variable % SteadyValues) ) THEN
        x0 => Solver % Variable % SteadyValues
        Stat = .TRUE.
      END IF
    ELSE 
      IF( ASSOCIATED(Solver % Variable % NonlinValues)) THEN
        x0 => Solver % Variable % NonlinValues
        Stat = .TRUE.
      END IF
    END IF
    IF(Stat .AND. .NOT. SkipConstraints ) THEN
      IF (SIZE(x0) /= SIZE(x)) CALL Warn('ComputeChange','Possible mismatch in length of vectors!')
    END IF

    IF( ResidualMode ) THEN
      IF(Relax .AND. RelaxBefore) THEN
        x(1:n) = x0(1:n) + Relaxation*x(1:n)
      ELSE
        x(1:n) = x0(1:n) + x(1:n)
      END IF
    ELSE 
      IF(Relax .AND. RelaxBefore) THEN
        x(1:n) = (1-Relaxation)*x0(1:n) + Relaxation*x(1:n)
      END IF
    END IF

    IF(SteadyState) THEN
      PrevNorm = Solver % Variable % PrevNorm
    ELSE
      PrevNorm = Solver % Variable % Norm
    END IF

    Norm = ComputeNorm(Solver, n, x)
    Solver % Variable % Norm = Norm

    !--------------------------------------------------------------------------
    ! The norm should be bounded in order to reach convergence
    !--------------------------------------------------------------------------
!   IF( ISNAN(Norm) ) THEN ! ISNAN not avaiable in all compilers
    IF( Norm /= Norm ) THEN
      CALL Fatal('ComputeChange','Norm of solution appears to be NaN')
    END IF

    MaxNorm = ListGetCReal( SolverParams, &
        'Nonlinear System Max Norm', Stat )
    IF( .NOT. Stat) MaxNorm = HUGE( Norm )

    IF(  Norm > MaxNorm ) THEN
      WRITE( Message, *) 'Computed Norm:',Norm
      CALL Info('ComputeChange',Message)
      CALL Fatal('ComputeChange','Norm of solution exceeded given bounds')
    END IF
      
    SELECT CASE( ConvergenceType )
        
    CASE('residual')
      !--------------------------------------------------------------------------
      ! x is solution of A(x0)x=b(x0), thus residual should be reall r=b(x)-A(x)x 
      ! Instead we use r=b(x0)-A(x0)x0 which unfortunately is one step behind.
      !--------------------------------------------------------------------------
      IF(PRESENT(Matrix)) THEN
        A => Matrix
      ELSE
        A => Solver % Matrix
      END IF

      IF(PRESENT(RHS)) THEN
        b => RHS
      ELSE
        b => Solver % Matrix % rhs
      END IF
      
      ALLOCATE(r(n)); r=0._dp

      IF (Parenv % Pes>1) THEN
        ALLOCATE( TmpRHSVec(n), TmpXVec(n) )

        nn = A % ParMatrix % SplittedMatrix % InsideMatrix % NumberOfRows

        TmpRhsVec = b
        CALL ParallelInitSolve( A, tmpXVec, TmpRhsVec, r)

        tmpXvec = x0(1:n)
        CALL ParallelVector(a,TmpXvec)
        CALL ParallelVector(A,tmpRhsvec)

        CALL ParallelMatrixVector(A, TmpXvec, r)
        DO i=1,nn
          r(i) = r(i) - tmprhsvec(i)
        END DO

        Change = ParallelNorm(nn,r)
        bNorm =  ParallelNorm(nn,tmpRhsVec)
      ELSE
        CALL MatrixVectorMultiply( A, x0, r)
        DO i=1,n
          r(i) = r(i) - b(i)
        END DO
        Change = ComputeNorm(Solver, n, r)
        bNorm  = ComputeNorm(Solver, n, b)
      END IF


      IF(.NOT. ConvergenceAbsolute) THEN
        IF(bNorm > 0.0) THEN
          Change = Change / bNorm
        END IF
      END IF
      DEALLOCATE(r)
      
    CASE('linear system residual')
      !--------------------------------------------------------------------------
      ! Here the true linear system residual r=b(x0)-A(x0)x is computed.
      ! This option is useful for certain special solvers.  
      !--------------------------------------------------------------------------
      A => Solver % Matrix
      b => Solver % Matrix % rhs
      
      IF (ParEnv % Pes > 1) THEN

        ALLOCATE( TmpRHSVec(n), TmpXVec(n), TmpRVec(n) )
        TmpRHSVec(1:n) = b(1:n)
        TmpXVec(1:n) = x(1:n)
        TmpRVec(1:n) = 0.0d0

        CALL ParallelVector(A, TmpRHSVec)
        CALL ParallelVector(A, TmpXVec)       
        CALL SParMatrixVector( TmpXVec, TmpRVec, ipar )
 
        nn = A % ParMatrix % SplittedMatrix % InsideMatrix % NumberOfRows

        DO i=1, nn
          TmpRVec(i) = TmpRHSVec(i) - TmpRVec(i)
        END DO

        Change = ParallelNorm( nn, TmpRVec )

        IF(.NOT. ConvergenceAbsolute) THEN
          bNorm = ParallelNorm( nn, TmpRHSVec )
          IF(bNorm > 0.0) THEN
            Change = Change / bNorm
          END IF
        END IF
        DEALLOCATE( TmpRHSVec, TmpXVec, TmpRVec )
      ELSE	
        ALLOCATE(r(n)) 
        CALL MatrixVectorMultiply( A, x, r)
        DO i=1,n
          r(i) = r(i) - b(i)
        END DO
        Change = SQRT( DOT_PRODUCT( r(1:n), r(1:n) ) )
        IF(.NOT. ConvergenceAbsolute) THEN
          bNorm = SQRT( DOT_PRODUCT( b(1:n), b(1:n) ) )
          IF(bNorm > 0.0) THEN
            Change = Change / bNorm
          END IF
        END IF
        DEALLOCATE(r)	
      END IF
      
    CASE('solution')      

      ALLOCATE(r(n))
      r = x(1:n)-x0(1:n)
      Change = ComputeNorm(Solver, n, r)
      IF( .NOT. ConvergenceAbsolute .AND. Norm + PrevNorm > 0.0) THEN
        Change = Change * 2.0/ (Norm+PrevNorm)
      END IF
      DEALLOCATE(r)      

    CASE('norm')
      Change = ABS( Norm-PrevNorm )
      IF( .NOT. ConvergenceAbsolute .AND. Norm + PrevNorm > 0.0) THEN
        Change = Change * 2.0/ (Norm+PrevNorm)
      END IF
      
    CASE DEFAULT
      CALL Warn('ComputeChange','Unknown convergence measure: '//TRIM(ConvergenceType))    
      
    END SELECT
    
    !--------------------------------------------------------------------------
    ! Check for convergence: 0/1
    !--------------------------------------------------------------------------
    IF(SteadyState) THEN
      Solver % Variable % SteadyChange = Change
      Tolerance = ListGetCReal( SolverParams,'Steady State Convergence Tolerance',Stat)
      IF( Stat ) THEN
        IF( Change <= Tolerance ) THEN
          Solver % Variable % SteadyConverged = 1
        ELSE
          Solver % Variable % SteadyConverged = 0
        END IF          
      END IF
    ELSE
      Solver % Variable % NonlinChange = Change
      Tolerance = ListGetCReal( SolverParams,'Nonlinear System Convergence Tolerance',Stat)
      IF( Stat ) THEN
        IF( Change <= Tolerance ) THEN
          Solver % Variable % NonlinConverged = 1
        ELSE
          Solver % Variable % NonlinConverged = 0
        END IF          
      END IF
      Tolerance = ListGetCReal( SolverParams,'Nonlinear System Exit Condition',Stat)
      IF( Stat .AND. Tolerance > 0.0 ) THEN
        CALL Info('ComputeChange','Nonlinear iteration condition enforced by exit condition')
        Solver % Variable % NonlinConverged = 1
      END IF
    END IF

    IF(Relax .AND. .NOT. RelaxBefore) THEN
      x(1:n) = (1-Relaxation)*x0(1:n) + Relaxation*x(1:n)
      Solver % Variable % Norm = ComputeNorm(Solver,n,x)
    END IF

    ! Steady state output is done in MainUtils
    SolverName = ListGetString( SolverParams, 'Equation',Stat)
    IF(.NOT. Stat) SolverName = Solver % Variable % Name

    IF(SteadyState) THEN        
      WRITE( Message, '(a,g15.8,g15.8,a)') &
         'SS (ITER='//TRIM(i2s(IterNo))//') (NRM,RELC): (',Norm, Change,&
          ' ) :: '// TRIM(SolverName)
    ELSE
      WRITE( Message, '(a,g15.8,g15.8,a)') &
         'NS (ITER='//TRIM(i2s(IterNo))//') (NRM,RELC): (',Norm, Change,&
          ' ) :: '// TRIM(SolverName)
    END IF
    CALL Info( 'ComputeChange', Message, Level=3 )


    ! The update of exported variables may be done internally to allow some nonlinear features	   
    ! or in steady state level to allow coupling to other solvers.
    !-----------------------------------------------------------------------------------------
    DoIt = .FALSE.
    IF( SteadyState ) THEN 
      DoIt = ListGetLogical( SolverParams,&
          'Update Exported Variables',Stat)
    ELSE 
      DoIt = ListGetLogical( SolverParams,&
          'Nonlinear Update Exported Variables',Stat)
    END IF
    IF( DoIt ) CALL UpdateExportedVariables( Solver )	


    ! Optional a posteriori scaling for the computed fields
    ! May be usefull for some floating systems where one want to impose some intergral 
    ! constraints without actually using them. Then first use just one Dirichlet point
    ! and then fix the level a posteriori using this condition. 
    !----------------------------------------------------------------------------------
    DoIt = .FALSE.
    IF( SteadyState ) THEN 
      DoIt = ListGetLogical( SolverParams,&
          'Nonlinear System Set Average Solution',Stat)
    ELSE 
      DoIt = ListGetLogical( SolverParams,&
          'Linear System Set Average Solution',Stat)
    END IF
    IF( DoIt ) THEN
      IF( ParEnv % PEs > 1 ) THEN
        CALL Fatal('ComputeChange','Setting average value not implemented in parallel!')
      END IF
      Ctarget = ListGetCReal( SolverParams,'Average Solution Value',Stat)      
      str = ListGetString( SolverParams,'Average Solution Weight Variable',Stat)
      IF( Stat ) THEN
        WeightVar => VariableGet( Solver % Mesh % Variables, str )
        IF( .NOT. ASSOCIATED( WeightVar ) ) THEN
          CALL Fatal('ComputeChange','> Average Solution Weight < missing: '//TRIM(str))
        END IF
        IF( SIZE(x) /= SIZE(WeightVar % Values ) ) THEN
          CALL Fatal('ComputeChange','Field and weight size mismatch: '//TRIM(str))          
        END IF
        Ctarget = Ctarget - SUM( WeightVar % Values * x ) / SUM( WeightVar % Values )
      ELSE
        Ctarget = Ctarget - SUM(x) / SIZE(x)
      END IF
      x = x + Ctarget
    END IF

    ! Compute derivative of solution with time i.e. velocity 
    ! For 2nd order schemes there is direct pointer to the velocity component
    ! Thus only 1st order schemes need to be computed.
    IF( Solver % TimeOrder == 1) THEN
      DoIt = .FALSE.
      IF( SteadyState ) THEN
        DoIt = ListGetLogical( SolverParams,'Calculate Velocity',Stat)
      ELSE
        DoIt = ListGetLogical( SolverParams,'Nonlinear Calculate Velocity',Stat)
      END IF
      IF( DoIt ) THEN
        TimestepVar => VariableGet( Solver % Mesh % Variables, 'timestep size' )
        dt = TimestepVar % Values(1) 
        str = TRIM( Solver % Variable % Name ) // ' Velocity'
        VeloVar => VariableGet( Solver % Mesh % Variables, str )        
        VeloVar % Values = (x - Solver % Variable % PrevValues(:,1)) / dt
      END IF
    END IF


    ! Calculate derivative a.k.a. sensitivity
    IF( SteadyState ) THEN

      IF( ListGetLogical( SolverParams,'Calculate Derivative',Stat) ) THEN

        IF( IterNo > 1 ) THEN
          TimestepVar => VariableGet( Solver % Mesh % Variables, 'derivative eps' )
          IF( ASSOCIATED( TimestepVar ) ) THEN
            eps = TimestepVar % Values(1)
            Stat = .TRUE.
          ELSE
            eps = ListGetCReal( SolverParams,'derivative eps',Stat)
          END IF
          IF(.NOT. Stat) THEN
            CALL Warn('ComputeChange','Derivative Eps not given, using one')
            Eps = 1.0_dp
          END IF

          str = TRIM( Solver % Variable % Name ) // ' Derivative'
          VeloVar => VariableGet( Solver % Mesh % Variables, str )
          IF( ASSOCIATED( VeloVar ) ) THEN
            CALL Info('ComputeChange','Computing variable:'//TRIM(str))
            VeloVar % Values = (x - x0) / eps
          ELSE
            CALL Warn('ComputeChange','Derivative variable not present')
          END IF
        END IF
      END IF
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE ComputeChange
!------------------------------------------------------------------------------
    


!------------------------------------------------------------------------------
!> Checks stepsize of a linear system so that the error has decreased.
!> Various indicatators and search algorithms have been implemented,
!------------------------------------------------------------------------------
  FUNCTION CheckStepSize(Solver,FirstIter,&
      nsize,values,values0) RESULT ( ReduceStep ) 
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    LOGICAL :: FirstIter
    INTEGER, OPTIONAL :: nsize
    REAL(KIND=dp), OPTIONAL, TARGET :: values(:), values0(:)
    LOGICAL :: ReduceStep
!------------------------------------------------------------------------------
    INTEGER :: MaxTests=0,tests,MaxNonlinIter,NonlinIter, Dofs
    REAL(KIND=dp) :: Residual0, Residual1, Residual
    INTEGER :: i,n,m,ForceDof, SearchMode, CostMode, iter = 0
    TYPE(Matrix_t), POINTER :: A, MP
    TYPE(Variable_t), POINTER :: IterVar, Var
    REAL(KIND=dp), POINTER :: b(:), x(:), x0(:), r(:), x1(:), x2(:), mr(:), mx(:), mb(:)
    REAL(KIND=dp) :: Norm, PrevNorm, rNorm, bNorm, Relaxation, Alpha, Myy, &
        NonlinTol, LineTol, Cost0(4), Cost1(4), Cost(4), OrthoCoeff, x0norm, x1norm, Change, &
        LinTol
    REAL(KIND=dp), ALLOCATABLE :: TempRHS(:)
    INTEGER, POINTER :: Indexes(:)
    LOGICAL :: Stat, Init, Newton, Ortho, Debug, SaveToFile
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, ConvergenceType, FileName
    TYPE(ValueList_t), POINTER :: SolverParams


    SAVE SolverParams, Alpha, Myy, Relaxation, MaxTests, tests, &
        Residual, NonlinTol, LinTol, x1, x0, LineTol, CostMode, SearchMode, &
        Cost0, Residual0, Cost1, n, Dofs, ForceDof, Ortho, Newton, &
        ConvergenceType, Norm, PrevNorm, iter, FileName, SaveToFile

    Debug = .FALSE.
    
    SolverParams => Solver % Values
    Var => Solver % Variable
    Dofs = Var % Dofs
 
    IF(PRESENT(values)) THEN
      x => values
    ELSE 
      x => Var % Values      
    END IF


    ! Assembly the vectors, if needed, and 
    ! also at first time get the line search parameters.
    !----------------------------------------------------
    IF( FirstIter ) THEN
      CALL Info('CheckStepSize','Initializing step-size search',Level=6)

      IF(PRESENT(nsize)) THEN
        n = nsize
      ELSE
        n = SIZE(x)
      END IF
      
      IF( ASSOCIATED( x0 ) ) THEN
        IF( SIZE(x0) /= n ) DEALLOCATE( x0 )
      END IF
      
      IF( PRESENT( values0 ) ) THEN
        x0 => values0 
      ELSE
        IF( .NOT. ASSOCIATED( x0 ) ) THEN
          ALLOCATE( x0(n) )
        END IF
      END IF
      
      IF( ASSOCIATED( x1 ) ) THEN
        IF( SIZE(x1) /= n ) DEALLOCATE( x1 )
      END IF
      IF( .NOT. ASSOCIATED( x1 ) ) THEN
        ALLOCATE( x1(n) )
      END IF

      Norm = 0.0_dp
      Var % NonlinConverged = 0
      Var % NonlinChange = 1.0_dp
      
      ! 1 - Residual norm : |Ax-b| 
      ! 2 - Quadratic functional : x^T(Ax-2b)/2
      ! 3 - Weighted residual : x^T(Ax-b)
      ! 4 - Lumped force : SUM(r_i)
      !------------------------------------------------------------
      CostMode = ListGetInteger( SolverParams,'Nonlinear System Linesearch Cost Mode',Stat)
      IF(.NOT. Stat) CostMode = 1
      
      ! 1 - Armijo-Goldstein criterion & successive relaxation 
      ! 2 - Minimize cost by bisection
      ! 3 - Find the zero cost by bisection
      !------------------------------------------------------------
      SearchMode = ListGetInteger( SolverParams,'Nonlinear System Linesearch Search Mode',Stat)
      IF(.NOT. Stat) SearchMode = 1

      ! Should the search direction be orthogonalized 
      !-----------------------------------------------------------
      Ortho = ListGetLogical( SolverParams,'Nonlinear System Linesearch Orthogonal',Stat)

      ! Is the outer ieration performed by Newton i.e. the search 
      ! should always be differential. 
      !-----------------------------------------------------------
      Newton = ListGetLogical( SolverParams,'Nonlinear System Linesearch Newton',Stat)

      NonlinTol = ListGetConstReal( SolverParams, &
          'Nonlinear System Convergence Tolerance', Stat )
      LinTol = ListGetConstReal( SolverParams, &
          'Linear System Convergence Tolerance', Stat )

      MaxNonlinIter = ListGetInteger( SolverParams,&
            'Nonlinear System Max Iterations',Stat)
      IF( MaxNonlinIter <= 2 ) THEN
        CALL Warn('CheckStepSize','For linesearch to work the nonlin iterations should be larger: '&
            //I2S(MaxNonlinIter))
      END IF

      ConvergenceType = ListGetString(SolverParams,&
          'Nonlinear System Convergence Measure',Stat)
      IF(.NOT. Stat) ConvergenceType = 'norm'

      ! Parameters related to line search algorithms
      !------------------------------------------------
      MaxTests = ListGetInteger( SolverParams,&
          'Nonlinear System Linesearch Iterations',Stat)
      IF( .NOT. Stat ) MaxTests = 10

      Myy = ListGetConstReal( SolverParams, &
          'Nonlinear System Linesearch Limit', Stat )
      IF(.NOT. Stat) Myy = 0.5_dp

      Relaxation = ListGetConstReal( SolverParams, &
          'Nonlinear System Linesearch Factor', Stat )
      IF(.NOT. Stat) Relaxation = 0.5_dp

      LineTol = ListGetConstReal( SolverParams, &
          'Nonlinear System Linesearch Tolerance', Stat )

      ForceDof = ListGetInteger( SolverParams, &
          'Nonlinear System Linesearch Force Index', Stat )
      IF(.NOT. Stat) ForceDof = Dofs

      FileName = ListGetString( SolverParams, &
          'Nonlinear System Linesearch Filename', SaveToFile )
      
      ! Computation of nonlinear change is now done with this routine
      ! so skip computing the change in the standard slot.
      !---------------------------------------------------------------
      CALL ListAddLogical(SolverParams,&
          'Skip Compute Nonlinear Change',.TRUE.)
    END IF

    !--------------------------------------------------------------------------
    ! This is the real residual: r=b-Ax
    ! We hope to roughly minimize L2 norm of r, or some related quantity
    !--------------------------------------------------------------------------
    A => Solver % Matrix
    b => Solver % Matrix % rhs

    ALLOCATE(r(n))
    IF (Parenv % Pes>1) THEN
      ALLOCATE(TempRHS(n))
      r = 0._dp
      TempRHS(1:n) = b(1:n)
      CALL ParallelInitSolve( A, x, TempRHS, r )

      MP => ParallelMatrix(A,mx,mb,mr)
      m = MP % NumberOfRows

      CALL ParallelMatrixVector( A, mx, r)
      r(1:m) = r(1:m) - TempRHS(1:m)
      Residual= ParallelNorm(n,r)
    ELSE
      CALL MatrixVectorMultiply( A, x, r)
      r(1:n) = r(1:n) - b(1:n)
      Residual = ComputeNorm(Solver, n, r)
    END IF

    ! Currently we compute all the costs to make it easier to study the 
    ! behavior of different measures when doing linesearch.
    IF( .TRUE. ) THEN
      Cost(1) = Residual
      Cost(2) = SUM( 0.5_dp * x(1:n) * ( r(1:n) - b(1:n) ) )
      Cost(3) = SUM( x(1:n) * r(1:n) )
      Cost(4) = SUM( r(ForceDof::Dofs) )
    ELSE      
      IF( CostMode == 1 ) THEN
        Cost(1) = Residual
      ELSE IF( CostMode == 2 ) THEN
        Cost(2) = SUM( 0.5_dp * x(1:n) * ( r(1:n) - b(1:n) ) )
      ELSE IF( CostMode == 3 ) THEN
        Cost(3) = SUM( x(1:n) * r(1:n) )
      ELSE IF( CostMode == 4 ) THEN
        Cost(4) = SUM( r(ForceDof::Dofs) )
      ELSE
        CALL Fatal('CheckStepSize','Unknown CostMode: '//TRIM(I2S(SearchMode)))
      END IF
      DEALLOCATE(r)
    END IF

    WRITE( Message,'(A,4ES15.7)') 'Cost: ',Cost
    CALL Info('CheckStepSize',Message,Level=8)

    ! At first iteration we cannot really do anything but go further 
    ! and save the reference residual for comparison.
    !-----------------------------------------------------------------------------
    IF( FirstIter ) THEN
      Tests = 0
      ReduceStep = .FALSE.
      x0(1:n) = x(1:n)
      Cost0 = Cost
      Residual0 = Residual

      IF( Debug ) THEN
        PRINT *,'x0 range: ',MINVAL(x0),MAXVAL(x0)
        PRINT *,'b0 range: ',MINVAL(b),MAXVAL(b)
        PRINT *,'Cost0: ',Cost0
      END IF

      IF( SaveToFile ) THEN
        CALL Info('CheckStepSize','Saving step information into file: '&
            //TRIM(FileName),Level=10)
        OPEN( 10, FILE = FileName, STATUS='UNKNOWN' )
        i = 0
        WRITE (10,'(2I6,5ES15.7)') Tests,i,Alpha,Cost
        CLOSE( 10 )
      END IF


      RETURN
    END IF

    Tests = Tests + 1

    IF( Tests == 1 ) THEN
      ! Started with no relaxation
      !---------------------------
      x1 = x
      Alpha = 1.0_dp
      Cost1 = Cost

      ! This is just debugging code waiting to be reused
      IF( .FALSE. ) THEN
        iter = iter + 1

        PRINT *,'Iter: ',iter
        NULLIFY( x2 ) 
        ALLOCATE( x2(n/2) ) 
        x2 = x(1::2)
        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
            'xiter '//TRIM(I2S(iter)),1,x2,Solver % Variable % Perm ) 
        PRINT *,'Xiter range:',MINVAL(x2),MAXVAL(x2)
        NULLIFY(x2)
        
!        NULLIFY( x2 ) 
!        ALLOCATE( x2(n/2) ) 
!        x2 = x(2::2)
!        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
!            'yiter '//TRIM(I2S(iter)),1,x2,Solver % Variable % Perm ) 
!        NULLIFY(x2)
      END IF
      
      IF( Debug ) THEN
        PRINT *,'b1 range: ',MINVAL(b),MAXVAL(b)
        PRINT *,'x1 range: ',MINVAL(x1),MAXVAL(x1)
        PRINT *,'Cost1: ',Cost1
      END IF

      ! Orthonormalization:
      ! The direction 'x0' has already been exhausted so remove that from 'x1'
      !-----------------------------------------------------------------------
      x0norm = ComputeNorm( Solver, n, x0 )
      IF( Ortho ) THEN
        IF( x0norm > EPSILON( x0norm ) ) THEN
          OrthoCoeff = SUM(x1*x0) / ( x0norm**2 )
          x1 = x1 - OrthoCoeff * x0
        END IF
      ELSE
        ! This basically checks whether the new and old solution is so 
        ! close that there is no point of finding better solution.
        x1 = x1 - x0 
        x1norm = ComputeNorm(Solver, n, x1)
        IF( x1norm < LinTol * x0norm ) THEN
          ReduceStep = .FALSE.
          GOTO 100
        END IF
      END IF

      IF( Debug ) THEN
        PRINT *,'x1 range orto: ',MINVAL(x1),MAXVAL(x1)
      END IF
    END IF

    ! Armijo-GoldStein Criterion for accepting stepsize
    !-----------------------------------------------------------------
    IF( SearchMode == 1 ) THEN
      ReduceStep = ArmijoGoldsteinSearch(Tests, Alpha )
    ELSE IF( SearchMode == 2 ) THEN
      ReduceStep = BisectMinimumSearch(Tests, Alpha) 
    ELSE IF( SearchMode == 3 ) THEN
      ReduceStep = BisectZeroSearch(Tests, Alpha)       
    ELSE
      CALL Fatal('CheckStepSize','Unknown SearchMode: '//TRIM(I2S(SearchMode)))
    END IF


    IF( SaveToFile ) THEN
      CALL Info('CheckStepSize','Saving step information into file: '&
          //TRIM(FileName),Level=10)
      OPEN( 10, FILE = FileName, POSITION='APPEND',STATUS='OLD' )
      IF( ReduceStep ) THEN
        i = 0
      ELSE
        i = 1
      END IF

      WRITE (10,'(2I6,5ES13.6)') Tests,i,Alpha,Cost
      CLOSE( 10 )
    END IF



100 IF( ReduceStep ) THEN
      IF( Tests >= MaxTests .AND. ReduceStep ) THEN
        CALL Fatal('CheckStepSize','Maximum number of linesearch steps taken without success!')
        ReduceStep = .FALSE.
      END IF
      
      ! New candidate 
      x(1:n) = x0(1:n) + Alpha * x1(1:n)

      WRITE(Message,'(A,I0,A,g15.6)') 'Step ',Tests,' rejected, trying new extent: ',Alpha
      CALL Info( 'CheckStepSize',Message,Level=6 )
    ELSE ! accept step
      WRITE(Message,'(A,I0,A,g15.6)') 'Step ',Tests,' accepted with extent: ',Alpha
      CALL Info( 'CheckStepSize',Message,Level=6 )
      
      ! Chosen candidate
      x(1:n) = x0(1:n) + Alpha * x1(1:n)

      PrevNorm = Norm
      Norm = ComputeNorm(Solver, n, x)

      IF( ConvergenceType == 'residual') THEN
        bNorm = ComputeNorm(Solver, n, b)
        IF( bNorm > 0.0_dp ) Change = Residual / bNorm
      ELSE
        Change = ABS( Norm-PrevNorm )
        IF( Norm + PrevNorm > 0.0) THEN
          Change = Change * 2.0 / ( Norm + PrevNorm )
        END IF
      END IF

      Solver % Variable % NonlinChange = Change
      Solver % Variable % Norm = Norm
     
      IF( Solver % Variable % NonlinChange <  NonlinTol ) THEN
        Solver % Variable % NonlinConverged = 1
      END IF
      
      SolverName = ListGetString( SolverParams, 'Equation',Stat)
      IF(.NOT. Stat) SolverName = Solver % Variable % Name
            
      IterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter')
      m = NINT(IterVar % Values(1))
      
      ! This replaces the standard error output usually written by the ComputeChange
      WRITE( Message, '(a,g15.8,g15.8,a)') &
          'NS (ITER='//TRIM(i2s(m))//') (NRM,RELC): (',Norm, Change, &
          ' ) :: '// TRIM(SolverName)
      CALL Info( 'CheckStepSize', Message, Level=3 )       
      
      WRITE(Message,'(A,I0,A,g15.6)') 'Step accepted after ',tests,' trials: ',Alpha
      CALL Info( 'CheckStepSize',Message,Level=5 )
      WRITE(Message,'(A,g15.6)') 'Previous cost:',Cost0(CostMode)
      CALL Info( 'CheckStepSize',Message,Level=6 )
      WRITE(Message,'(A,g15.6)') 'Initial cost: ',Cost1(CostMode)
      CALL Info( 'CheckStepSize',Message,Level=6 )
      WRITE(Message,'(A,g15.6)') 'Final cost:   ',Cost(CostMode)
      CALL Info( 'CheckStepSize',Message,Level=6 )
      
      Tests = 0
      x0 = x
      
      IF( Debug ) THEN
        PRINT *,'x0 range: ',MINVAL(x0),MAXVAL(x0)
        PRINT *,'Cost0: ',Cost0
        PRINT *,'Residual0: ',Residual0
      END IF

      IF( Newton ) FirstIter = .TRUE.

    END IF



  CONTAINS

!-----------------------------------------------------------------
!> Armijo-GoldStein Criterion for accepting stepsize
!-----------------------------------------------------------------
    FUNCTION ArmijoGoldsteinSearch(Tests,Alpha) RESULT ( ReduceStep )

      INTEGER :: Tests 
      REAL(KIND=dp) :: Alpha
      LOGICAL :: ReduceStep

      ReduceStep = ( Cost(CostMode) > ( 1.0_dp - Myy * Alpha ) * Cost0(CostMode) )
      IF( ReduceStep ) THEN
        Alpha = Alpha * Relaxation
      ELSE
        Cost0 = Cost
        Residual0 = Residual
      END IF

    END FUNCTION ArmijoGoldsteinSearch


!-------------------------------------------------------------------------------
!> Choose next parameter set from 1D bisection search
!-------------------------------------------------------------------------------

    FUNCTION BisectMinimumSearch(Tests, Alpha) RESULT ( ReduceStep ) 

      INTEGER :: Tests 
      REAL(KIND=dp) :: Alpha
      LOGICAL :: ReduceStep
      
      INTEGER :: i,j,k
      REAL(KIND=dp) :: step, p(3),c(3),r(3),raid,beta 
      
      SAVE step, p, c, r

      ReduceStep = .TRUE.
      
      IF(Tests == 1) THEN
        p(1) = 0.0_dp
        c(1) = Cost0(CostMode)
        r(1) = Residual0

        p(2) = 1.0_dp
        c(2) = Cost(CostMode)
        r(2) = Residual
        
        step = 0.25_dp
        Alpha = 0.5_dp
        RETURN
      ELSE 
        p(3) = Alpha
        c(3) = Cost(CostMode) 
        r(3) = Residual
      END IF

      
     ! Order the previous points so that p1 < p2 < p3
      DO k=1,2 
        DO i=k+1,3
          IF(p(i) < p(k)) THEN
            raid = p(k)
            p(k) = p(i)
            p(i) = raid

            raid = c(k)
            c(k) = c(i)
            c(i) = raid

            raid = r(k)
            r(k) = r(i)
            r(i) = raid
          END IF
        END DO
      END DO
      
      IF( Debug ) THEN
        PRINT *,'Bisect p:',p
        PRINT *,'Bisect c:',c
        PRINT *,'Bisect r:',r
      END IF

      ! The value of alpha already known accurately
      IF( MAXVAL(p)-MINVAL(p) < LineTol ) THEN
        ! PRINT *,'cond1'
        ReduceStep = .FALSE.
      END IF

      ! The value of cost function small compared to absolute value of it
      IF( MAXVAL(c)-MINVAL(c) < LineTol * MINVAL( ABS(c) ) ) THEN
        ! PRINT *,'cond2'
        ReduceStep = .FALSE.
      END IF

      ! We can also use the residual as criterion for stopping
      IF( Residual < LineTol * Residual0 ) THEN
        ! PRINT *,'cond3'
        ReduceStep = .FALSE.
      END IF

      ! Of these choose the one with smallest cost
      IF( .NOT. ReduceStep ) THEN
        i = 1
        DO k=2,3
          IF( c(k) < c(i) ) i = k
        END DO

        Alpha = p(i)
        Residual0 = r(i)
        Cost0(CostMode) = c(i)
        ! PRINT *,'Choosing i',i,Alpha,Residual0,Cost0

        RETURN
      END IF


      ! Monotonic line segment
      IF( (c(2)-c(1))*(c(3)-c(2)) > 0.0) THEN
        IF(c(3) < c(1)) THEN
          Alpha = p(3) + SIGN(step,p(3)-p(1))
          c(1) = c(3)
          p(1) = p(3)
          r(1) = r(3)
        ELSE
          Alpha = p(1) + SIGN(step,p(1)-p(3))
        END IF
      ELSE IF(c(2) < c(1) .OR. c(2) < c(3)) THEN 
        IF(c(3) < c(1)) THEN
          c(1) = c(3)
          p(1) = p(3)
          r(1) = r(3)
        END IF
        step = (p(2)-p(1))/2.0d0
        Alpha = p(1) + SIGN(step,p(2)-p(1))
      ELSE  
        IF( Debug ) THEN
          PRINT *,'p:',p
          PRINT *,'c:',c,Cost0
          PRINT *,'r:',r,Residual0
          PRINT *,'dc',c(2)-c(1),c(3)-c(2)
        END IF

        IF( MINVAL ( c ) < Cost0(CostMode) ) THEN
          i = 1
          DO k=2,3
            IF( c(k) < c(i) ) i = k
          END DO
          Alpha = p(i)
          Cost0(CostMode) = c(i)
          Residual0 = r(i)
         
          CALL Warn('BisectSearch','Bisection method improved but faced local maximium')
          ReduceStep = .FALSE.
        ELSE 
          IF( MINVAL ( r ) < Residual0 ) THEN
            CALL Warn('BisectSearch','Bisection method improved but faced local maximium')
          ELSE 
            CALL Warn('BisectSearch','Bisection method cannot handle local maxima')
          END IF

          i = 1
          DO k=2,3
            IF( r(k) < r(i) ) i = k
          END DO
          Alpha = p(i)
          Cost0(CostMode) = c(i)
          Residual0 = r(i)         
        END IF

        ReduceStep = .FALSE.
      END IF

      ! Because alpha should be in limit [0,1] make the corrections
      ! If the orthogonalization is used then we don't have the luxury
      ! of setting the extent as nicely.
      !------------------------------------------------------------
      IF( .NOT. Ortho ) THEN
        beta = alpha
        k = 0
        IF( Alpha < -EPSILON( Alpha ) ) THEN
          IF( p(1) < EPSILON( Alpha ) ) THEN
            step = (p(2)-p(1))/2.0_dp
            Alpha = p(1) + step
            k = 1
          ELSE
            Alpha = 0.0_dp
            k = 1
          END IF
        ELSE IF( Alpha > 1.0_dp + EPSILON( Alpha ) ) THEN
          IF( p(3) > 1.0_dp - EPSILON( Alpha ) ) THEN
            step = (p(3)-p(2))/2.0_dp
            Alpha = p(2) + step
            k = 2
          ELSE
            Alpha = 1.0_dp
            k = 3 
          END IF
        END IF
        
!        IF( ABS( beta-alpha) > TINY(alpha)) PRINT *,'Extent change',Beta,Alpha
      END IF

    END FUNCTION BisectMinimumSearch


!-------------------------------------------------------------------------------
!> Choose next parameter set from 1D bisection search
!-------------------------------------------------------------------------------
    FUNCTION BisectZeroSearch(Tests, Alpha) RESULT ( ReduceStep ) 

      INTEGER :: Tests 
      REAL(KIND=dp) :: Alpha
      LOGICAL :: ReduceStep
      
      INTEGER :: i,j,k
      REAL(KIND=dp) :: step, p(3),c(3),paid,caid,beta 
      
      SAVE step, p, c

      ReduceStep = .TRUE.
      
      IF(Tests == 1) THEN
        p(1) = 0.0_dp
        c(1) = Cost0(CostMode)
        
        p(2) = 1.0_dp
        c(2) = Cost1(CostMode)
        
        IF( Cost0(CostMode) * Cost1(CostMode) > 0.0_dp ) THEN
          CALL Warn('CostSearch','Lumped forces should have different sign!')
        END IF

        Alpha = 0.5_dp
        RETURN
      ELSE 
        p(3) = Alpha
        c(3) = Cost(CostMode) 
      END IF
      
     ! Order the previous points so that p1 < p2 < p3
      DO k=1,2 
        DO i=k+1,3
          IF(p(i) < p(k)) THEN
            paid = p(k)
            p(k) = p(i)
            p(i) = paid
            caid = c(k)
            c(k) = c(i)
            c(i) = caid
          END IF
        END DO
      END DO

      IF( Debug ) THEN
        PRINT *,'Cost p:',p
        PRINT *,'Cost c:',c
      END IF

      IF( p(3)-p(1) < LineTol ) THEN
        ReduceStep = .FALSE.
        RETURN
      END IF

      ! Zero value is between 1st interval
      IF( c(1)*c(2) < 0.0_dp ) THEN
        Alpha = (p(1)+p(2))/2.0_dp
      ELSE IF ( c(2)*c(3) < 0.0_dp ) THEN
        Alpha = (p(2)+p(3))/2.0_dp

        ! We don't need 1st values, but we do need 3rd 
        p(1) = p(3)
        c(1) = c(3)
      ELSE
        CALL Fatal('ForceSearch','Lumped forces should have different sign!')
      END IF
      
    END FUNCTION BisectZeroSearch

!------------------------------------------------------------------------------
  END FUNCTION CheckStepSize
!------------------------------------------------------------------------------

    

!------------------------------------------------------------------------------
!> Computing nodal weight may be good when one needs to transform nodal 
!> information back to continuous fields by dividing with the nodal weight. 
!> Active either for the permutation defined by the primary variable of the 
!> solver, or for a permutation vector defined by an optional flag that
!> is used as a mask to define the set of active nodes.
!------------------------------------------------------------------------------
  SUBROUTINE CalculateNodalWeights(Solver,WeightAtBoundary,&
      Perm,VarName,Var)
!------------------------------------------------------------------------------
    IMPLICIT NONE 
    TYPE(Solver_t) :: Solver
    LOGICAL :: WeightAtBoundary
    INTEGER, POINTER, OPTIONAL :: Perm(:)
    CHARACTER(*), OPTIONAL :: VarName
    TYPE(Variable_t), POINTER, OPTIONAL :: Var
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: IntVarName
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: WeightsVar
    TYPE(ValueList_t), POINTER :: ElemParams
    REAL(KIND=dp), POINTER :: Weights(:), Solution(:)    
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    INTEGER ::k, e, t, n, ElemStart, ElemFin, Coordsys
    INTEGER, POINTER :: IntPerm(:), Indexes(:),LocalIndexes(:)
    REAL(KIND=dp) :: u,v,w,s,detJ
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    LOGICAL :: GotIt, stat, VariableOutput, UseMask, RequireLogical, Hit
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)


    Mesh => Solver % Mesh
    CoordSys = CurrentCoordinateSystem()

    NULLIFY( WeightsVar ) 
    IF( PRESENT( VarName ) ) THEN
      IntVarName = VarName
    ELSE IF ( WeightAtBoundary ) THEN
      IntVarName = GetVarName(Solver % Variable) // ' Boundary Weights'
    ELSE
      IntVarName = GetVarName(Solver % Variable) // ' Weights'
    END IF
    WeightsVar => VariableGet( Mesh % Variables, IntVarName )

    IF( WeightAtBoundary ) THEN
      ElemStart = Mesh % NumberOfBulkElements + 1
      ElemFin = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      UseMask = ListCheckPresentAnyBC( CurrentModel, IntVarName )
    ELSE
      ElemStart = 1
      ElemFin = Mesh % NumberOfBulkElements 
      UseMask = ListCheckPresentAnyBodyForce( CurrentModel, IntVarName )
    END IF

    RequireLogical = .FALSE.
    NULLIFY( IntPerm ) 
    IF ( .NOT. ASSOCIATED(WeightsVar) ) THEN
      IF( PRESENT( Perm ) ) THEN
        IntPerm => Perm 
      ELSE
        IntPerm => Solver % Variable % Perm
      END IF
      IF( ASSOCIATED( IntPerm ) ) THEN
	NULLIFY( Solution )
	n = MAXVAL( IntPerm ) 
        ALLOCATE( Solution(n))
        Solution = 0.0d0        
        CALL VariableAdd( Mesh % Variables, Mesh, Solver,&
            IntVarName, 1, Solution, IntPerm )
        NULLIFY( Solution )
      ELSE
        CALL Warn('CalculateNodalWeights','Permutation vector not present?')
        RETURN
      END IF
      WeightsVar => VariableGet( Mesh % Variables, IntVarName )
    END IF

    IF( .NOT. ASSOCIATED( WeightsVar ) ) THEN
      CALL Fatal('CalculateNodalWeights','Solution variable not present?')
    END IF
    Weights => WeightsVar % Values
    IntPerm => WeightsVar % Perm
    IF ( .NOT. ASSOCIATED(Weights) ) THEN
      CALL Warn('CalculateNodalWeights','Solution vector not present?')
      RETURN
    ELSE
      IF( PRESENT( Var) ) Var => WeightsVar
    END IF

    CALL Info('ComputeNodalWeights',&
        'Computing weights for solver to variable: '//TRIM(IntVarName))
    n = Mesh % MaxElementNodes

    ALLOCATE(Basis(n), ElementNodes % x(n), ElementNodes % y(n), &
        ElementNodes % z(n), LocalIndexes(n) )
    Weights = 0.0_dp

    DO e=ElemStart,ElemFin

      Element => Mesh % Elements( e )
      Indexes => Element % NodeIndexes

      n = Element % TYPE % NumberOfNodes
      LocalIndexes(1:n) = IntPerm( Indexes ) 
      IF( ANY( LocalIndexes(1:n) == 0 ) ) CYCLE

      IF( UseMask ) THEN
        Hit = .FALSE.
        IF( WeightAtBoundary ) THEN
          DO k=1,CurrentModel % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(k) % Tag ) THEN
              Hit = .TRUE.
              EXIT
            END IF
          END DO
          IF( .NOT. Hit ) CYCLE
          ElemParams => CurrentModel % BCs(k) % Values
        ELSE
          ElemParams => CurrentModel % Bodies(Element % BodyId) % Values
        END IF
        IF( RequireLogical ) THEN
          IF( .NOT. ListGetLogical( ElemParams, IntVarName, Stat ) ) CYCLE
        ELSE
          IF( .NOT. ListCheckPresent( ElemParams, IntVarName ) ) CYCLE
        END IF
      END IF

      n = Element % TYPE % NumberOfNodes
      ElementNodes % x(1:n) = Mesh % Nodes % x(Indexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(Indexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(Indexes)

      IntegStuff = GaussPoints( Element )

      DO t=1,IntegStuff % n        
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)

        stat = ElementInfo( Element, ElementNodes, U, V, W, detJ, Basis )

        IF ( CoordSys /= Cartesian ) THEN
          X = SUM( ElementNodes % X(1:n) * Basis(1:n) )
          Y = SUM( ElementNodes % Y(1:n) * Basis(1:n) )
          Z = SUM( ElementNodes % Z(1:n) * Basis(1:n) )
          CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
          s = s * SqrtMetric
        END IF
        
        Weights( LocalIndexes(1:n) ) = &
            Weights( LocalIndexes(1:n) ) + s * detJ * Basis(1:n)
      END DO

    END DO

    DEALLOCATE(Basis, ElementNodes % x, ElementNodes % y, &
        ElementNodes % z, LocalIndexes )

    CALL Info('ComputeNodalWeights','All done')

  END SUBROUTINE CalculateNodalWeights
!------------------------------------------------------------------------------




  !> Calcualte the number of separature pieces in a serial mesh.
  !> This could be used to detect problems in mesh when suspecting 
  !> floating parts not fixed by any BC, for example.
  !---------------------------------------------------------------------------------
  SUBROUTINE CalculateMeshPieces( Mesh )

    TYPE(Mesh_t), POINTER :: Mesh

    LOGICAL :: Ready
    INTEGER :: i,j,k,n,t,MinIndex,MaxIndex,Loop,NoPieces
    INTEGER, ALLOCATABLE :: MeshPiece(:),PiecePerm(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: Indexes(:)
    TYPE(Variable_t), POINTER :: Var

    IF( ParEnv % PEs > 1 ) THEN
      CALL Warn('CalculateMeshPieces','Implemented only for serial meshes, doing nothing!')
      RETURN
    END IF

    n = Mesh % NumberOfNodes
    ALLOCATE( MeshPiece( n ) ) 
    MeshPiece = 0

    ! Only set the piece for the nodes that are used by some element
    ! For others the marker will remain zero. 
    DO t = 1, Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)        
      Indexes => Element % NodeIndexes
      MeshPiece( Indexes ) = 1
    END DO    
    j = 0
    DO i = 1, n
      IF( MeshPiece(i) > 0 ) THEN
        j = j + 1
        MeshPiece(i) = j
      END IF
    END DO

    CALL Info('CalculateMeshPieces','Number of non-body nodes in mesh is '//TRIM(I2S(n-j)))
    
    ! We go through the elements and set all the piece indexes to minimimum index
    ! until the mesh is unchanged. Thereafter the whole piece will have the minimum index
    ! of the piece.
    Ready = .FALSE.
    Loop = 0
    DO WHILE(.NOT. Ready) 
      Ready = .TRUE.
      DO t = 1, Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)        
        Indexes => Element % NodeIndexes

        MinIndex = MINVAL( MeshPiece( Indexes ) )
        MaxIndex = MAXVAL( MeshPiece( Indexes ) )
        IF( MaxIndex > MinIndex ) THEN
          MeshPiece( Indexes ) = MinIndex
          Ready = .FALSE.
        END IF
      END DO
      Loop = Loop + 1
    END DO
    CALL Info('CalculateMeshPieces','Mesh coloring loops: '//TRIM(I2S(Loop)))

    ! If the maximum index is one then for sure there is only one body
    IF( MaxIndex == 1 ) THEN
      CALL Info('CalculateMeshPieces','Mesh consists of single body!')
      RETURN
    END IF

    ! Compute the true number of different pieces
    ALLOCATE( PiecePerm( MaxIndex ) ) 
    PiecePerm = 0
    NoPieces = 0
    DO i = 1, n
      j = MeshPiece(i) 
      IF( j == 0 ) CYCLE
      IF( PiecePerm(j) == 0 ) THEN
        NoPieces = NoPieces + 1
        PiecePerm(j) = NoPieces 
      END IF
    END DO
    CALL Info('CalculateMeshPieces','Number of seperate pieces in mesh is '//TRIM(I2S(NoPieces)))


    ! Save the mesh piece field to > mesh piece < 
    Var => VariableGet( Mesh % Variables,'Mesh Piece' )
    IF(.NOT. ASSOCIATED( Var ) ) THEN      
      CALL VariableAddVector ( Mesh % Variables,Mesh, CurrentModel % Solver,'Mesh Piece' )
      Var => VariableGet( Mesh % Variables,'Mesh Piece' )
    END IF

    IF( .NOT. ASSOCIATED( Var ) ) THEN
      CALL Fatal('CalculateMeshPieces','Could not get handle to variable > Mesh Piece <')
    END IF

    DO i = 1, n
      j = Var % Perm( i ) 
      IF( j == 0 ) CYCLE
      IF( MeshPiece(i) > 0 ) THEN
        Var % Values( j ) = 1.0_dp * PiecePerm( MeshPiece( i ) )
      ELSE
        Var % Values( j ) = 0.0_dp
      END IF
    END DO
    CALL Info('CalculateMeshPieces','Saving mesh piece field to: mesh piece')
  
  END SUBROUTINE CalculateMeshPieces



!------------------------------------------------------------------------------
!> Compute weights of entities i.e. their area and volume in the mesh.
!------------------------------------------------------------------------------
  SUBROUTINE CalculateEntityWeights(Model, Mesh)
!------------------------------------------------------------------------------
    IMPLICIT NONE 
    TYPE(Model_t) :: Model 
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element, Left, Right
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    INTEGER ::i,j,k, e,t, n, Coordsys, TrueOwner, bc_id, bf_id, mat_id, body_id, &
        maxsize, ierr, PotOwner
    INTEGER :: NoBC, NoBodies, NoBF, NoMat
    INTEGER, POINTER :: Indexes(:)
    REAL(KIND=dp) :: u,v,w,s,detJ
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    REAL(KIND=dp), POINTER :: bc_weights(:),body_weights(:),&
        mat_weights(:),bf_weights(:),tmp_weights(:)
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3), &
        Coeff
    LOGICAL :: Found, Stat, BodyElem


    CoordSys = CurrentCoordinateSystem()

    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Warn('CalculateEntityWeights','Mesh not associated!')
      RETURN
    END IF

    CALL Info('ComputeNodalWeights','Computing weights for the mesh entities',Level=6)
    n = Mesh % MaxElementNodes

    NoBC = Model % NumberOfBCs
    NoBodies = Model % NumberOfBodies
    NoMat = Model % NumberOfMaterials
    NoBF = Model % NumberOfBodyForces
    

    ALLOCATE(Basis(n), &
        ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )

    IF( .NOT. Mesh % EntityWeightsComputed ) THEN
      IF( NoBC > 0 ) ALLOCATE( Mesh % BCWeight(NoBC) )
      IF( NoBodies > 0 ) ALLOCATE( Mesh % BodyWeight(NoBodies ) ) 
      IF( NoMat > 0 ) ALLOCATE( Mesh % MaterialWeight(NoMat) ) 
      IF( NoBF > 0 ) ALLOCATE( Mesh % BodyForceWeight(NoBF ) )       
    END IF

    IF( NoBC > 0 ) THEN
      bc_weights => Mesh % BCWeight
      bc_weights(1:NoBC ) = 0.0_dp
    END IF

    IF( NoBodies > 0 ) THEN
      body_weights => Mesh % BodyWeight
      body_weights(1:NoBodies ) = 0.0_dp
    END IF

    IF( NoMat > 0 ) THEN
      mat_weights => Mesh % MaterialWeight
      mat_weights(1:NoMat ) = 0.0_dp
    END IF

    IF( NoBF > 0 ) THEN
      bf_weights => Mesh % BodyForceWeight       
      bf_weights(1:NoBF ) = 0.0_dp
    END IF

    DO e=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      bf_id = 0
      mat_id = 0
      body_id = 0
      bc_id = 0
      Coeff = 1.0_dp

      BodyElem = ( e <= Mesh % NumberOfBulkElements ) 
      Element => Mesh % Elements( e )

      IF( BodyElem ) THEN
        body_id = Element % BodyId
        bf_id = ListGetInteger( Model % Bodies(body_id) % Values,&
            'Body Force',Found)
        mat_id = ListGetInteger( Model % Bodies(body_id) % Values,&
            'Material',Found)
      ELSE
        Found = .FALSE.
        DO bc_id = 1,Model % NumberOfBCs
          Found = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) 
          IF( Found ) EXIT
        END DO
        IF(.NOT. Found) CYCLE
      END IF

      Coeff = 1.0_dp
      
      ! In parallel compute the weight only at their true owners.
      ! Therefore cycle the halo elements. For shared BCs 
      ! take only half of the weight. 
      IF( ParEnv % PEs > 1 ) THEN
        IF( BodyElem ) THEN
          IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE
        ELSE
          TrueOwner = 0
          PotOwner = 0
          Left => Element % BoundaryInfo % Left
          IF( ASSOCIATED( Left ) ) THEN
            PotOwner = PotOwner + 1
            IF( Left % PartIndex == ParEnv % MyPe ) TrueOwner = TrueOwner + 1
          END IF
          Right => Element % BoundaryInfo % Right
          IF( ASSOCIATED( Right ) ) THEN
            PotOwner = PotOwner + 1
            IF( Right % PartIndex == ParEnv % MyPe ) TrueOwner = TrueOwner + 1
          END IF
          IF( PotOwner > 0 ) THEN
            IF( TrueOwner == 0 ) CYCLE
            Coeff = 1.0_dp * TrueOwner / PotOwner
          END IF
        END IF
      END IF

      Indexes => Element % NodeIndexes

      n = Element % TYPE % NumberOfNodes
      ElementNodes % x(1:n) = Mesh % Nodes % x(Indexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(Indexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(Indexes)

      IntegStuff = GaussPoints( Element )

      DO t=1,IntegStuff % n        
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)

        stat = ElementInfo( Element, ElementNodes, U, V, W, detJ, Basis )
        S = Coeff * DetJ * IntegStuff % s(t)

        IF ( CoordSys /= Cartesian ) THEN
          X = SUM( ElementNodes % X(1:n) * Basis(1:n) )
          Y = SUM( ElementNodes % Y(1:n) * Basis(1:n) )
          Z = SUM( ElementNodes % Z(1:n) * Basis(1:n) )
          CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
          s = s * 2 * PI * SqrtMetric
        END IF
                
        IF( bc_id > 0 .AND. bc_id <= NoBC) &
            bc_weights( bc_id ) = bc_weights( bc_id ) + s
        IF( body_id > 0 .AND. body_id <= NoBodies ) &
            body_weights( body_id ) = body_weights( body_id ) + s
        IF( mat_id > 0 .AND. mat_id <= NoMat ) &
            mat_weights( mat_id ) = mat_weights( mat_id ) + s
        IF( bf_id > 0 .AND. bf_id <= NoBF ) &
            bf_weights( bf_id ) = bf_weights( bf_id ) + s
      END DO

    END DO


    IF( ParEnv % PEs > 1 ) THEN
      maxsize = MAX( Model % NumberOfBCs, Model % NumberOfBodies ) 
      ALLOCATE( tmp_weights( maxsize ) ) 
      tmp_weights = 0.0_dp

      IF( NoBC > 0 ) THEN
        tmp_weights(1:NoBC ) = bc_weights
        CALL MPI_ALLREDUCE( tmp_weights, bc_weights, NoBC, &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
      END IF
      IF( NoBF > 0 ) THEN
        tmp_weights(1:NoBF ) = bf_weights
        CALL MPI_ALLREDUCE( tmp_weights, bf_weights, NoBF, &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
      END IF
      IF( NoBodies > 0 ) THEN
        tmp_weights(1:NoBodies ) = body_weights
        CALL MPI_ALLREDUCE( tmp_weights, body_weights, NoBodies, &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
      END IF
      IF( NoMat > 0 ) THEN
        tmp_weights(1:NoMat ) = mat_weights
        CALL MPI_ALLREDUCE( tmp_weights, mat_weights, NoMat, &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
      END IF
      DEALLOCATE( tmp_weights )
    END IF

    IF( ParEnv % MyPe == 0 ) THEN
      DO i = 1, NoBC
        PRINT *,'BC weight:',i,bc_weights(i)
      END DO
      DO i = 1, NoBF
        PRINT *,'BF weight:',i,bf_weights(i)
      END DO
      DO i = 1, NoBodies
        PRINT *,'Body weight:',i,body_weights(i)
      END DO
      DO i = 1, NoMat
        PRINT *,'Mat weight:',i,mat_weights(i)
      END DO
    END IF

    DEALLOCATE(Basis, &
        ElementNodes % x, ElementNodes % y, ElementNodes % z )

    Mesh % EntityWeightsComputed = .TRUE.

    CALL Info('CalculateEntityWeights','All done',Level=10)

  END SUBROUTINE CalculateEntityWeights
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>  Scale system Ax = b as:
!>  (DAD)y = Db, where D = 1/SQRT(Diag(A)), and y = D^-1 x
!------------------------------------------------------------------------------
  SUBROUTINE ScaleLinearSystem(Solver,A,b,x,DiagScaling, & 
          ApplyScaling,RhsScaling,ConstraintScaling)

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), OPTIONAL :: b(:),x(:)
    REAL(KIND=dp), OPTIONAL, TARGET :: DiagScaling(:)
    LOGICAL, OPTIONAL :: ApplyScaling, RhsScaling,ConstraintScaling
    INTEGER :: n,i,j
    REAL(KIND=dp) :: bnorm,s
    COMPLEX(KIND=dp) :: DiagC
    LOGICAL :: ComplexMatrix, DoRHS, DoCM
    REAL(KIND=dp), POINTER  :: Diag(:)

    TYPE(Matrix_t), POINTER :: CM

    n = A % NumberOfRows

    IF( PRESENT( DiagScaling ) ) THEN
      Diag => DiagScaling 
    ELSE
      IF(.NOT. ASSOCIATED(A % DiagScaling)) THEN
        ALLOCATE( A % DiagScaling(n) ) 
      END IF
      Diag => A % DiagScaling
      Diag = 0._dp
    
      ComplexMatrix = Solver % Matrix % COMPLEX
    
      IF ( ComplexMatrix ) THEN
        DO i=1,n,2
          j = A % Diag(i)
          IF(j>0) THEN
            Diag(i)   = A % Values(j)
            Diag(i+1) = A % Values(j+1)
          ELSE
            Diag(i)=0._dp;Diag(i+1)=0._dp
          END IF
        END DO
      ELSE
        DO i=1,n
          j=A % Diag(i)
          IF (j>0) Diag(i) = A % Values(j)
        END DO
      END IF
      
      IF ( ParEnv % PEs > 1 ) CALL ParallelSumVector(A, Diag)

      IF ( ComplexMatrix ) THEN
        DO i=1,n,2
          DiagC = CMPLX(Diag(i),-Diag(i+1),KIND=dp)
          IF (ABS(DiagC)/=0._dp) THEN
            Diag(i)   = 1.0_dp/SQRT(ABS(DiagC))
            Diag(i+1) = 1.0_dp/SQRT(ABS(DiagC))
          ELSE
            Diag(i)   = 1.0_dp; Diag(i+1) = 1.0_dp
          END IF
        END DO
      ELSE
        s = 0
        IF (ANY(ABS(Diag)<=TINY(bnorm))) s=1
        s = ParallelReduction(s,2) 

        IF(s/=0) THEN 
          DO i=1,n
            IF ( ABS(Diag(i)) <= TINY(bnorm) ) THEN
              Diag(i) = SUM( ABS(A % Values(A % Rows(i):A % Rows(i+1)-1)) )
            ELSE
              j = A % Diag(i)
              IF (j>0) Diag(i) = A % Values(j)
            END IF
          END DO
          IF ( ParEnv % PEs > 1 ) CALL ParallelSumVector(A, Diag)
        END IF

        DO i=1,n
          IF ( ABS(Diag(i)) > TINY(bnorm) ) THEN
            Diag(i) = 1.0_dp / SQRT(ABS(Diag(i)))
          ELSE
            Diag(i) = 1.0_dp
          END IF
        END DO
      END IF
    END IF    


    ! Optionally we may just create the diag and leave the scaling undone
    !--------------------------------------------------------------------
    IF( PRESENT( ApplyScaling ) ) THEN
      IF(.NOT. ApplyScaling ) RETURN
    END IF
    
    DO i=1,n
      DO j = A % Rows(i), A % Rows(i+1)-1
        A % Values(j) = A % Values(j) * &
            ( Diag(i) * Diag(A % Cols(j)) )
      END DO
    END DO
    
    IF ( ASSOCIATED( A % PrecValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % PrecValues(j) = A % PrecValues(j) * &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
      END IF
    END IF

    IF ( ASSOCIATED( A % MassValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % MassValues(j) = A % MassValues(j) * &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
      END IF
    END IF
    
    IF ( ASSOCIATED( A % DampValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % DampValues(j) = A % DampValues(j) * &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
      END IF
    END IF

    DoCM=.FALSE.
    IF(PRESENT(ConstraintScaling)) DoCm=ConstraintScaling

    IF(doCM) THEN
      CM => A % ConstraintMatrix
      IF (ASSOCIATED(CM)) THEN
        DO i=1,CM % NumberOFRows
          DO j=CM % Rows(i), CM % Rows(i+1)-1
            CM % Values(j) = CM % Values(j) * Diag(CM % Cols(j))
          END DO
        END DO
      END IF
    END IF

    ! Scale r.h.s. and initial guess
    !--------------------------------
    A % RhsScaling=1._dp
    IF( PRESENT( b ) ) THEN      
      b(1:n) = b(1:n) * Diag(1:n)
      DoRHS = .TRUE.
      IF (PRESENT(RhsScaling)) DoRHS = RhsScaling
      IF (DoRHS) THEN
        bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))
      ELSE
        bnorm = 1._dp
      END IF
      A % RhsScaling = bnorm

      Diag(1:n) = Diag(1:n) * bnorm
      b(1:n) = b(1:n) / bnorm
      IF( PRESENT( x) ) THEN
        x(1:n) = x(1:n) / Diag(1:n)
      END IF
    END IF
!-----------------------------------------------------------------------------
  END SUBROUTINE ScaleLinearSystem
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!>   Equilibrate the rows of the coefficient matrix A to
!>   minimize the condition number. The associated rhs vector f is also scaled.
!------------------------------------------------------------------------------
  SUBROUTINE RowEquilibration( A, f, Parallel )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: f(:)
    LOGICAL :: Parallel
!-----------------------------------------------------------------------------
    LOGICAL :: ComplexMatrix
    INTEGER :: i, j, n 
    REAL(kind=dp) :: norm, tmp
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:), Diag(:)
!-------------------------------------------------------------------------
    n = A % NumberOfRows
    ComplexMatrix = A % COMPLEX

    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    IF( .NOT. ASSOCIATED(A % DiagScaling) ) THEN
      ALLOCATE( A % DiagScaling(n) ) 
    END IF
    Diag => A % DiagScaling    
    
    Diag = 0.0d0
    norm = 0.0d0

    !---------------------------------------------
    ! Compute 1-norm of each row
    !---------------------------------------------
    IF (ComplexMatrix) THEN
      DO i=1,n,2
        tmp = 0.0d0
        DO j=Rows(i),Rows(i+1)-1,2
          tmp = tmp + ABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
        END DO

        IF ( .NOT. Parallel ) THEN
          IF (tmp > norm) norm = tmp        
        END IF

        IF (tmp > 0.0d0) THEN
          Diag(i) = tmp
          Diag(i+1) = tmp
        END IF
      END DO
    ELSE
      DO i=1,n
        tmp = 0.0d0
        DO j=Rows(i),Rows(i+1)-1        
          tmp = tmp + ABS(Values(j))          
        END DO

        IF ( .NOT. Parallel ) THEN
          IF (tmp > norm) norm = tmp        
        END IF

        IF (tmp > 0.0d0) Diag(i) = tmp       
      END DO
    END IF

    IF (Parallel) THEN
      CALL ParallelSumVector(A, Diag)
      norm = ParallelReduction(MAXVAL(Diag(1:n)),2)
    END IF

    !--------------------------------------------------
    ! Now, define the scaling matrix by inversion and 
    ! perform the actual scaling of the linear system
    !--------------------------------------------------
    IF (ComplexMatrix) THEN    
      DO i=1,n,2
        IF (Diag(i) > 0.0d0) THEN
          Diag(i) = 1.0d0/Diag(i)
          Diag(i+1) = 1.0d0/Diag(i+1)
        ELSE
          Diag(i) = 1.0d0
          Diag(i+1) = 1.0d0
        END IF
      END DO
    ELSE
      DO i=1,n      
        IF (Diag(i) > 0.0d0) THEN
          Diag(i) = 1.0d0/Diag(i)
        ELSE
          Diag(i) = 1.0d0
        END IF
      END DO
    END IF

    DO i=1,n    
      DO j=Rows(i),Rows(i+1)-1
        Values(j) = Values(j) * Diag(i)
      END DO
      f(i) = Diag(i) * f(i)
    END DO
    
    WRITE( Message, * ) 'Unscaled matrix norm: ', norm    
    CALL Info( 'OptimalMatrixScaling', Message, Level=5 )

!------------------------------------------------------------------------------
  END SUBROUTINE RowEquilibration
!------------------------------------------------------------------------------


  
!--------------------------------------------------------------
!>  Scale the system back to original.
!--------------------------------------------------------------
  SUBROUTINE BackScaleLinearSystem(Solver,A,b,x,DiagScaling,ConstraintScaling) 

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), OPTIONAL :: b(:),x(:)
    LOGICAL, OPTIONAL :: ConstraintScaling
    REAL(KIND=dp), OPTIONAL, TARGET :: DiagScaling(:)

    REAL(KIND=dp), POINTER :: Diag(:)
    REAL(KIND=dp) :: bnorm
    INTEGER :: n,i,j
    LOGICAL :: doCM

    TYPE(Matrix_t), POINTER :: CM
    
    n = A % NumberOfRows
    
    IF( PRESENT( DiagScaling ) ) THEN
      Diag => DiagScaling
    ELSE  
      Diag => A % DiagScaling
    END IF

    IF(.NOT. ASSOCIATED( Diag ) ) THEN
      CALL Warn('BackScaleLinearSystem','Diag not associated!')
      RETURN
    END IF
    IF( SIZE( Diag ) /= n ) THEN
      CALL Fatal('BackScaleLinearSystem','Diag of wrong size!')
    END IF 

    IF( PRESENT( b ) ) THEN
! 
!      Solve x:  INV(D)x = y, scale b back to orig
!      -------------------------------------------
      IF( PRESENT( x ) ) THEN
        x(1:n) = x(1:n) * Diag(1:n)
      END IF
      bnorm = A % RhsScaling
      Diag(1:n) = Diag(1:n) / bnorm
      b(1:n) = b(1:n) / Diag(1:n) * bnorm
    END IF
    
    DO i=1,Solver % NOFEigenValues
      !
      !           Solve x:  INV(D)x = y
      !           --------------------------
      IF ( Solver % Matrix % COMPLEX ) THEN
        Solver % Variable % EigenVectors(i,1:n/2) = &
            Solver % Variable % EigenVectors(i,1:n/2) * Diag(1:n:2)
      ELSE
        Solver % Variable % EigenVectors(i,1:n) = &
            Solver % Variable % EigenVectors(i,1:n) * Diag(1:n)
      END IF
    END DO
    
    DO i=1,n
      DO j=A % Rows(i), A % Rows(i+1)-1
        A % Values(j) = A % Values(j) / (Diag(i) * Diag(A % Cols(j)))
      END DO
    END DO
    
    IF ( ASSOCIATED( A % PrecValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % PrecValues(j) = A % PrecValues(j) / &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
      END IF
    END IF
    IF ( ASSOCIATED( A % MassValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % MassValues(j) = A % MassValues(j) / &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
      END IF
    END IF
    
    IF ( ASSOCIATED( A % DampValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % DampValues(j) = A % DampValues(j) / &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
      END IF
    END IF

    doCM=.FALSE.
    IF(PRESENT(ConstraintScaling)) doCM=ConstraintScaling
    IF(doCM) THEN
      CM => A % ConstraintMatrix
      IF (ASSOCIATED(CM)) THEN
        DO i=1,CM % NumberOFRows
          DO j=CM % Rows(i), CM % Rows(i+1)-1
            CM % Values(j) = CM % Values(j) / ( Diag(CM % Cols(j)) )
          END DO
        END DO
      END IF
    END IF

    A % RhsScaling=1._dp
    DEALLOCATE(A % DiagScaling); A % DiagScaling=>NULL()
    
  END SUBROUTINE BackScaleLinearSystem


!------------------------------------------------------------------------------
!> Scale the linear system back to original when the linear
!> system scaling has been done by row equilibration.
!------------------------------------------------------------------------------
  SUBROUTINE ReverseRowEquilibration( A, f )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: f(:)
!-----------------------------------------------------------------------------
    INTEGER :: i, j, n
    INTEGER, POINTER :: Rows(:)
    REAL(KIND=dp), POINTER :: Values(:), Diag(:)
!-----------------------------------------------------------------------------
    n = A % NumberOfRows
    Diag => A % DiagScaling   
    Values => A % Values
    Rows => A % Rows

    IF(.NOT. ASSOCIATED( Diag ) ) THEN
      CALL Fatal('ReverseRowEquilibration','Diag not associated!')
    END IF
    IF( SIZE( Diag ) /= n ) THEN
      CALL Fatal('ReverseRowEquilibration','Diag of wrong size!')
    END IF 

    f(1:n) = f(1:n) / Diag(1:n)
    DO i=1,n    
      DO j = Rows(i), Rows(i+1)-1
        Values(j) = Values(j) / Diag(i)
      END DO
    END DO

    DEALLOCATE(A % DiagScaling)
    A % DiagScaling => NULL()

!------------------------------------------------------------------------------
  END SUBROUTINE ReverseRowEquilibration
!------------------------------------------------------------------------------


  SUBROUTINE CalculateLoads( Solver, Aaid, x, DOFs, UseBulkValues, NodalLoads ) 

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), POINTER  :: Aaid
    REAL(KIND=dp) CONTIG :: x(:)
    INTEGER :: DOFs
    LOGICAL :: UseBulkValues
    TYPE(Variable_t), POINTER :: NodalLoads

    REAL(KIND=dp), POINTER :: LoadValues(:)
    INTEGER :: i,j,k,l,m,ii,This,DOF
    REAL(KIND=dp), POINTER :: TempRHS(:), TempVector(:), SaveValues(:), Rhs(:), TempX(:)
    REAL(KIND=dp) :: Energy, Energy_im
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL :: Found, Rotated


    REAL(KIND=dp), ALLOCATABLE :: BoundarySum(:), BufReal(:)
    INTEGER, ALLOCATABLE :: BoundaryShared(:),BoundaryActive(:),DofSummed(:),BufInteg(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: bc, ind, NoBoundaryActive, NoBCs, ierr
    LOGICAL :: OnlyGivenBCs


    IF( .NOT. ASSOCIATED(NodalLoads) ) RETURN
    ALLOCATE( TempVector(Aaid % NumberOfRows) )

    IF( UseBulkValues ) THEN
      SaveValues => Aaid % Values
      Aaid % Values => Aaid % BulkValues
      Rhs => Aaid % BulkRHS
    ELSE
      Rhs => Aaid % Rhs
    END IF


    IF ( ParEnv % PEs > 1 ) THEN
      ALLOCATE(TempRHS(SIZE(Rhs)))
      TempRHS = Rhs 
      CALL ParallelInitSolve( Aaid, x, TempRHS, Tempvector )
      CALL ParallelMatrixVector( Aaid, x, TempVector, .TRUE. )
    ELSE
      CALL MatrixVectorMultiply( Aaid, x, TempVector )
    END IF

    IF( ListGetLogical(Solver % Values, 'Calculate Energy Norm', Found) ) THEN
      Energy = 0._dp
      IF( ListGetLogical(Solver % Values, 'Linear System Complex', Found) ) THEN
        Energy_im = 0._dp
        DO i = 1, (Aaid % NumberOfRows / 2)
          IF ( ParEnv % Pes>1 ) THEN
            IF ( Aaid% ParMatrix % ParallelInfo % &
              NeighbourList(2*(i-1)+1) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
          END IF
          Energy    = Energy    + x(2*(i-1)+1) * TempVector(2*(i-1)+1) - x(2*(i-1)+2) * TempVector(2*(i-1)+2)
          Energy_im = Energy_im + x(2*(i-1)+1) * TempVector(2*(i-1)+2) + x(2*(i-1)+2) * TempVector(2*(i-1)+1) 
       END DO
       Energy    = ParallelReduction(Energy)
       Energy_im = ParallelReduction(Energy_im)

       CALL ListAddConstReal( Solver % Values, 'Energy norm', Energy)
       CALL ListAddConstReal( Solver % Values, 'Energy norm im', Energy_im)

       WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy Norm'
       CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy )

       WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy Norm im'
       CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy_im )

       WRITE( Message, * ) 'Energy Norm: ', Energy, Energy_im
       CALL Info( 'SolveLinearSystem', Message )
     ELSE 
       DO i=1,Aaid % NumberOfRows
         IF ( ParEnv % Pes>1 ) THEN
           IF ( Aaid % ParMatrix % ParallelInfo % &
                NeighbourList(i) % Neighbours(1) /= Parenv % MyPE ) CYCLE
         END IF
         Energy = Energy + x(i)*TempVector(i)
      END DO
      Energy = ParallelReduction(Energy)
      CALL ListAddConstReal( Solver % Values, 'Energy norm', Energy )

      WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy Norm'
      CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy )

      WRITE( Message, * ) 'Energy Norm: ', Energy
      CALL Info( 'SolveLinearSystem', Message )
    END IF
  END IF

    IF ( ParEnv % PEs>1 ) THEN
      DO i=1,Aaid % NumberOfRows
        IF ( AAid % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % Mype ) THEN
          TempVector(i) = TempVector(i) - TempRHS(i)
        ELSE
          TempVector(i) = 0
        END IF
      END DO
      CALL ParallelSumVector( AAid, Tempvector )
      DEALLOCATE( TempRhs ) 
    ELSE
      TempVector = TempVector - RHS
    END IF


    NoBCs = CurrentModel % NumberOfBCs
    DO This=1,NoBCs
      Projector => CurrentModel  % BCs(This) % PMatrix
      IF (ASSOCIATED(Projector))THEN
        DO DOF=1,DOFs
          DO i=1,Projector % NumberOfRows
            ii = Projector % InvPerm(i)
            IF( ii == 0 ) CYCLE
            k = Solver % Variable % Perm(ii)
            IF(k<=0) CYCLE
            k = DOFs * (k-1) + DOF
            TempVector(k)=0

            DO l = Projector % Rows(i), Projector % Rows(i+1)-1
              IF ( Projector % Cols(l) <= 0 ) CYCLE
              m = Solver % Variable % Perm( Projector % Cols(l) )
              IF ( m > 0 ) THEN
                m = DOFs * (m-1) + DOF
                TempVector(k) = TempVector(k) + Projector % Values(l)*TempVector(m)
              END IF
            END DO
          END DO
        END DO
      END IF
    END DO

    DO i=1,SIZE( NodalLoads % Perm )
      IF ( NodalLoads % Perm(i)>0 .AND. Solver % Variable % Perm(i)>0 ) THEN
        DO j=1,DOFs
          NodalLoads % Values(DOFs*(NodalLoads % Perm(i)-1)+j) =  &
              TempVector(DOFs*(Solver % Variable % Perm(i)-1)+j)
        END DO
      END IF
    END DO
    DEALLOCATE( TempVector )


    IF( ListGetLogical( Solver % Values,'Calculate Boundary Fluxes',Found ) ) THEN
      CALL Info('CalculateLoads','Computing boundary fluxes from nodal loads',Level=6)

      IF( Solver % Mesh % MaxEdgeDofs > 1 .OR. Solver % Mesh % MaxFaceDOFs > 1 ) THEN
        CALL Warn('CalculateLoads','Boundary flux computation implemented only for nodes for now!')
      END IF

      ALLOCATE( BoundarySum( NoBCs * DOFs ), &
          BoundaryActive( NoBCs ), &
          BoundaryShared( NoBCs ), &
          DofSummed( MAXVAL( NodalLoads % Perm ) ) )
      BoundarySum = 0.0_dp
      BoundaryActive = 0
      BoundaryShared = 0
      DofSummed = 0

      OnlyGivenBCs = ListCheckPresentAnyBC( CurrentModel,'Calculate Boundary Flux')
      
      k = Solver % Mesh % NumberOfBulkElements
      DO i = k+1,k + Solver % Mesh % NumberOfBoundaryElements
        Element => Solver % Mesh % Elements(i)
        bc = Element % BoundaryInfo % Constraint
           
        IF( bc == 0 ) CYCLE

        IF( OnlyGivenBCs ) THEN
          IF (.NOT. ListGetLogical( CurrentModel % BCs(bc) % Values,&
              'Calculate Boundary Flux',Found) ) CYCLE
        END IF

        DO j=1,Element % TYPE % NumberOfNodes
          ind = NodalLoads % Perm( Element % NodeIndexes(j) )
          IF( ind == 0 ) CYCLE

          ! In this partition sum up only the true owners
          IF ( ParEnv % PEs>1 ) THEN
            IF ( AAid % ParallelInfo % NeighbourList(ind) % Neighbours(1) &
                /= ParEnv % Mype ) CYCLE
          END IF

          ! Only sum each entry once. If there is a conflict we cannot 
          ! really resolve it with the chosen method so just warn. 
          IF( DofSummed(ind) == 0 ) THEN
            BoundarySum( DOFs*(bc-1)+1 :DOFs*bc ) = BoundarySum( DOFs*(bc-1)+ 1:DOFs*bc ) + &
                NodalLoads % Values( DOFs*(ind-1) + 1: DOFs * ind )
            DofSummed( ind ) = bc
            BoundaryActive( bc ) = 1
          ELSE IF( bc /= DofSummed(ind) ) THEN
            BoundaryShared(bc) = 1
            BoundaryShared(DofSummed(ind)) = 1
          END IF
        END DO
      END DO
      

      NoBoundaryActive = 0
      IF( ParEnv % PEs > 1 ) THEN
        ALLOCATE( BufInteg( NoBCs ), BufReal( NoBCs * DOFs ) )

        BufInteg = BoundaryActive
        CALL MPI_ALLREDUCE( BufInteg, BoundaryActive, NoBCs, MPI_INTEGER, &
            MPI_SUM, ParEnv % ActiveComm, ierr )
        
        BufInteg = BoundaryShared
        CALL MPI_ALLREDUCE( BufInteg, BoundaryShared, NoBCs, MPI_INTEGER, &
            MPI_SUM, ParEnv % ActiveComm, ierr )

        BufReal = BoundarySum 
        CALL MPI_ALLREDUCE( BufReal, BoundarySum, DOFs * NoBCs, MPI_DOUBLE_PRECISION, &
            MPI_SUM, ParEnv % ActiveComm, ierr )

        DEALLOCATE( BufInteg, BufReal ) 
      END IF


      DO i=1,CurrentModel % NumberOfBCs 
        IF( BoundaryActive(i) == 0 ) CYCLE
        IF( BoundaryShared(i) > 0) THEN
          CALL Warn('CalculateLoads','Boundary '//TRIM(I2S(i))//' includes inseparable dofs!')
        END IF
        NoBoundaryActive = NoBoundaryActive + 1

        DO j=1,DOFs
          IF( Dofs == 1 ) THEN
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' Flux over BC '//TRIM(I2S(i))
          ELSE
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' '//TRIM(I2S(j))//' Flux over BC '//TRIM(I2S(i))
          END IF
          CALL ListAddConstReal( CurrentModel % Simulation, 'res: '//TRIM(Message), &
              BoundarySum(DOFs*(i-1)+j) )
          WRITE( Message,'(A,ES12.5)') TRIM(Message)//': ',BoundarySum(DOFs*(i-1)+j)
          CALL Info('CalculateLoads',Message,Level=6)
        END DO
      END DO
      
      IF( NoBoundaryActive > 1 ) THEN
        DO j=1,DOFs
          IF( Dofs == 1 ) THEN
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' Flux over all BCs'
          ELSE
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' '//TRIM(I2S(j))//' Flux over all BCs'
          END IF
          WRITE( Message,'(A,ES12.5)') TRIM(Message)//': ',SUM(BoundarySum(j::DOFs))
          CALL Info('CalculateLoads',Message,Level=6)
        END DO
      END IF
      
      DEALLOCATE( DofSummed, BoundaryShared, BoundaryActive, BoundarySum )      
    END IF


    IF( UseBulkValues ) THEN
      Aaid % Values => SaveValues
    END IF

  END SUBROUTINE CalculateLoads


!------------------------------------------------------------------------------
!> Prints the values of the CRS matrix to standard output.
!------------------------------------------------------------------------------
  SUBROUTINE PrintMatrix( A, Parallel, CNumbering,SaveMass, SaveDamp, SaveStiff)
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A            !< Structure holding matrix
    LOGICAL :: Parallel    !< are we in parallel mode?
    LOGICAL :: CNumbering  !< Continuous numbering ?
    LOGICAL, OPTIONAL :: SaveMass  !< Should we save the mass matrix
    LOGICAL, OPTIONAL :: SaveDamp  !< Should we save the damping matrix
    LOGICAL, OPTIONAL :: SaveStiff !< Should we save the stiffness matrix
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,IndMass,IndDamp,IndStiff,IndMax,row,col
    LOGICAL :: DoMass, DoDamp, DoStiff, Found
    REAL(KIND=dp) :: Vals(3)
    INTEGER, ALLOCATABLE :: Owner(:)

    DoMass = .FALSE.
    IF( PRESENT( SaveMass ) ) DoMass = SaveMass
    IF( DoMass .AND. .NOT. ASSOCIATED( A % MassValues ) ) THEN
      CALL Warn('CRS_PrintMatrix','Cannot save nonexisting mass matrix')
      DoMass = .FALSE. 
    END IF

    DoDamp = .FALSE.
    IF( PRESENT( SaveDamp ) ) DoDamp = SaveDamp
    IF( DoDamp .AND. .NOT. ASSOCIATED( A % DampValues ) ) THEN
      CALL Warn('CRS_PrintMatrix','Cannot save nonexisting damp matrix')
      DoDamp = .FALSE. 
    END IF

    DoStiff = .TRUE.
    IF( PRESENT( SaveStiff ) ) DoStiff = SaveStiff
    IF( DoStiff .AND. .NOT. ASSOCIATED( A % Values ) ) THEN
      CALL Warn('CRS_PrintMatrix','Cannot save nonexisting stiff matrix')
      DoStiff = .FALSE. 
    END IF

    IF(.NOT. (DoStiff .OR. DoDamp .OR. DoMass ) ) THEN
      CALL Warn('CRS_PrintMatrix','Saving just the topology!')
    END IF
    
    IndStiff = 0
    IndDamp = 0
    IndMass = 0

    IF( DoStiff ) IndStiff = 1
    IF( DoDamp ) IndDamp = IndStiff + 1
    IF( DoMass ) IndMass = IndDamp + 1
    IndMax = MAX( IndStiff, IndDamp, IndMass )

    IF (Parallel.AND.Cnumbering) THEN
      n = SIZE(A % ParallelInfo % GlobalDOFs)
  
      ALLOCATE( A % Gorder(n), Owner(n) )
      CALL ContinuousNumbering( A % ParallelInfo, &
          A % Perm, A % Gorder, Owner )
    END IF

    DO i=1,A % NumberOfRows
      row = i
      IF(Parallel) THEN
        IF(Cnumbering) THEN
          row = A % Gorder(i)
        ELSE 
          row = A % ParallelInfo % GlobalDOFs(i)
        END IF
      END IF
      DO j = A % Rows(i),A % Rows(i+1)-1

        col = A % Cols(j)
        IF(Parallel) THEN
          IF(Cnumbering) THEN
            col = A % Gorder(col)
          ELSE 
            col = A % ParallelInfo % GlobalDOFs(col)
          END IF
        END IF

        WRITE(1,'(I0,A,I0,A)',ADVANCE='NO') row,' ',col,' '

        IF( DoStiff ) THEN
          Vals(IndStiff) = A % Values(j)
        END IF
        IF( DoDamp ) THEN
          Vals(IndDamp) = A % DampValues(j)
        END IF
        IF( DoMass ) THEN
          Vals(IndMass) = A % MassValues(j)
        END IF

        IF( IndMax > 0 ) THEN
          WRITE(1,*) Vals(1:IndMax)          
        ELSE
          WRITE(1,'(A)') ' '
        END IF
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE  PrintMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Prints the values of the right-hand-side vector to standard output.
!------------------------------------------------------------------------------
  SUBROUTINE PrintRHS( A, Parallel, CNumbering )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A  !< Structure holding matrix
    LOGICAL :: Parallel, CNumbering
!------------------------------------------------------------------------------
    INTEGER :: i, row
    REAL(KIND=dp) :: Val

    DO i=1,A % NumberOfRows
      row = i
      IF(Parallel) THEN
        IF(Cnumbering) THEN
          row = A % Gorder(i)
        ELSE 
          row = A % ParallelInfo % GlobalDOFs(i)
        END IF
      END IF

      Val = A % Rhs(i)
      WRITE(1,'(I0,A)',ADVANCE='NO') row,' '
      IF( ABS( Val ) <= TINY( Val ) ) THEN
        WRITE(1,'(A)') '0.0'
      ELSE
        WRITE(1,*) Val
      END IF
    END DO

  END SUBROUTINE PrintRHS
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Solves a linear system and also calls the necessary preconditioning routines.
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE SolveLinearSystem( A, b, &
       x, Norm, DOFs, Solver, BulkMatrix )
!------------------------------------------------------------------------------
    REAL(KIND=dp) CONTIG :: b(:), x(:)
    REAL(KIND=dp) :: Norm
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: DOFs
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Matrix_t), OPTIONAL, POINTER :: BulkMatrix
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var, NodalLoads
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Relax,GotIt,Stat,ScaleSystem, EigenAnalysis, HarmonicAnalysis,&
               BackRotation, ApplyRowEquilibration, ApplyLimiter, Parallel, &
               SkipZeroRhs, ComplexSystem, ComputeChangeScaled, ConstraintModesAnalysis, &
               RecursiveAnalysis
    INTEGER :: n,i,j,k,l,ii,m,DOF,istat,this,mn
    CHARACTER(LEN=MAX_NAME_LEN) :: Method, Prec, ProcName, SaveSlot
    INTEGER(KIND=AddrInt) :: Proc
    REAL(KIND=dp), ALLOCATABLE, TARGET :: Px(:), &
                TempVector(:), TempRHS(:), NonlinVals(:)
    REAL(KIND=dp), POINTER :: Diag(:)
    REAL(KIND=dp) :: s,Relaxation,Beta,Gamma,bnorm,Energy,xn,bn
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: Aaid, Projector, MP
    REAL(KIND=dp), POINTER :: mx(:), mb(:), mr(:)
    TYPE(Variable_t), POINTER :: IterV

    INTERFACE 
       SUBROUTINE VankaCreate(A,Solver)
          USE Types
          TYPE(Matrix_t) :: A
          TYPE(Solver_t) :: Solver
       END SUBROUTINE VankaCreate

       SUBROUTINE CircuitPrecCreate(A,Solver)
          USE Types
          TYPE(Matrix_t) :: A
          TYPE(Solver_t) :: Solver
       END SUBROUTINE CircuitPrecCreate

       SUBROUTINE FetiSolver(A,x,b,Solver)
          USE Types
          TYPE(Matrix_t), POINTER :: A
          TYPE(Solver_t) :: Solver
          REAL(KIND=dp) :: x(:), b(:)
       END SUBROUTINE FetiSolver

       SUBROUTINE BlockSolveExt(A,x,b,Solver)
          USE Types
          TYPE(Matrix_t), POINTER :: A
          TYPE(Solver_t) :: Solver
          REAL(KIND=dp) :: x(:), b(:)
       END SUBROUTINE BlockSolveExt
    END INTERFACE
!------------------------------------------------------------------------------

   ComplexSystem=ListGetLogical( Solver % Values, 'Linear System Complex', GotIt )
   IF ( GotIt ) A % COMPLEX = ComplexSystem

!------------------------------------------------------------------------------
!   If parallel execution, check for parallel matrix initializations
!------------------------------------------------------------------------------
    IF ( ParEnv % Pes>1.AND..NOT. ASSOCIATED(A % ParMatrix) ) THEN
      CALL ParallelInitMatrix( Solver, A )
    END IF

   IF ( ListGetLogical( Solver % Values, 'Linear System Save',GotIt )) THEN
      saveslot = ListGetString( Solver % Values,'Linear System Save Slot', GotIt )
      IF(SaveSlot == 'linear solve') CALL SaveLinearSystem( Solver, A )
    END IF

!------------------------------------------------------------------------------

    Params => Solver % Values
    n = A % NumberOfRows

    BackRotation = ListGetLogical(Params,'Back Rotate N-T Solution',GotIt)
    IF (.NOT.GotIt) BackRotation=.TRUE.

    IF ( Solver % Matrix % Lumped .AND. Solver % TimeOrder == 1 ) THEN
       Method = ListGetString( Params, 'Timestepping Method', GotIt)
       IF (  Method == 'runge-kutta' .OR. Method == 'explicit euler' ) THEN
         ALLOCATE(Diag(n), TempRHS(n))

         TempRHS= b(1:n)
         Diag = A % Values(A % Diag)

         IF(ParEnv % Pes>1) THEN
           CALL ParallelSumVector(A,Diag)
           CALL ParallelSumVector(A,TempRHS)
         END IF

         DO i=1,n
            IF ( ABS(Diag(i)) /= 0._dp ) x(i) = TempRHS(i) / Diag(i)
         END DO

         DEALLOCATE(Diag, TempRHS)

         IF (BackRotation) CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )
         Norm = ComputeNorm(Solver, n, x) 
         RETURN
       END IF
    END IF
    
!------------------------------------------------------------------------------
!  These definitions are needed if chanching the iterative solver on-the-fly

    Solver % MultiGridSolver = ( ListGetString( Params, &
        'Linear System Solver', GotIt ) == 'multigrid' )
    Solver % MultiGridTotal = MAX( Solver % MultiGridTotal, &
        ListGetInteger( Params,'MG Levels', GotIt, minv=1 ) )
    Solver % MultiGridTotal = MAX( Solver % MultiGridTotal, &
        ListGetInteger( Params,'Multigrid Levels', GotIt, minv=1 ) )
    Solver % MultiGridLevel = Solver % MultigridTotal
!------------------------------------------------------------------------------

    ScaleSystem = ListGetLogical( Params, 'Linear System Scaling', GotIt )
    IF ( .NOT. GotIt  ) ScaleSystem = .TRUE.

    EigenAnalysis = Solver % NOFEigenValues > 0 .AND. &
        ListGetLogical( Params, 'Eigen Analysis',GotIt )
    
    ConstraintModesAnalysis = ListGetLogical( Params, &
        'Constraint Modes Analysis',GotIt )

    HarmonicAnalysis = Solver % NOFEigenValues>0 .AND. &
        ListGetLogical( Params, 'Harmonic Analysis',GotIt )

    ! These analyses types may require recursive strategies and may also have zero rhs
    RecursiveAnalysis = HarmonicAnalysis .OR. EigenAnalysis .OR. ConstraintModesAnalysis


    ApplyLimiter = ListGetLogical( Params,'Apply Limiter',GotIt ) 
    SkipZeroRhs = ListGetLogical( Params,'Skip Zero Rhs Test',GotIt )

    IF ( .NOT. ( RecursiveAnalysis .OR. ApplyLimiter .OR. SkipZeroRhs ) ) THEN
      bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))      
      IF ( bnorm <= TINY( bnorm) ) THEN
        CALL Info('SolveSystem','Solution trivially zero!')
        x = 0.0d0

        ! Increase the nonlinear counter since otherwise some stuff may stagnate
        ! Normally this is done within ComputeChange
        iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
        Solver % Variable % NonlinIter = iterV % Values(1)
        iterV % Values(1) = iterV % Values(1) + 1 
        Solver % Variable % Norm = 0.0_dp
        Solver % Variable % NonlinConverged = 1
     
        RETURN
      END IF
    END IF

    IF ( Solver % MultiGridLevel == -1  ) RETURN

    ! Set the flags to false to allow recursive strategies for these analysis types, little dirty...
    IF( RecursiveAnalysis ) THEN
      IF( HarmonicAnalysis ) CALL ListAddLogical( Solver % Values,'Harmonic Analysis',.FALSE.)
      IF( EigenAnalysis ) CALL ListAddLogical( Solver % Values,'Eigen Analysis',.FALSE.)
      IF( ConstraintModesAnalysis ) CALL ListAddLogical( Solver % Values,'Constraint Modes Analysis',.FALSE.)
    END IF


!------------------------------------------------------------------------------
!   If solving harmonic analysis go there:
!   --------------------------------------
    IF ( HarmonicAnalysis ) THEN
      CALL SolveHarmonicSystem( A, Solver )
    END IF


!   If solving eigensystem go there:
!   --------------------------------
    IF ( EigenAnalysis ) THEN
      IF ( ScaleSystem ) CALL ScaleLinearSystem(Solver, A )

      CALL SolveEigenSystem( &
          A, Solver %  NOFEigenValues, &
          Solver % Variable % EigenValues,       &
          Solver % Variable % EigenVectors, Solver )
      
      IF ( ScaleSystem ) CALL BackScaleLinearSystem( Solver, A ) 
      IF ( BackRotation ) CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )

      Norm = ComputeNorm(Solver,n,x)
      Solver % Variable % Norm = Norm
      
      CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
          Solver % Variable % Name )
    END IF


!   If solving constraint modes analysis go there:
!   ----------------------------------------------
    IF ( ConstraintModesAnalysis ) THEN
      CALL SolveConstraintModesSystem( A, Solver )

      IF ( BackRotation ) CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )
      
      Norm = ComputeNorm(Solver,n,x)
      Solver % Variable % Norm = Norm
      
      CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
          Solver % Variable % Name )
    END IF

    ! We have solved {harmonic,eigen,constraint} system and no need to continue further
    IF( RecursiveAnalysis ) THEN
      IF( HarmonicAnalysis ) CALL ListAddLogical( Solver % Values,'Harmonic Analysis',.TRUE.)
      IF( EigenAnalysis ) CALL ListAddLogical( Solver % Values,'Eigen Analysis',.TRUE.)
      IF( ConstraintModesAnalysis ) CALL ListAddLogical( Solver % Values,'Constraint Modes Analysis',.TRUE.)
      RETURN
    END IF


! Check whether b=0 since then equation Ax=b has only the trivial solution, x=0. 
! In case of a limiter one still may need to check the limiter for contact.
!-----------------------------------------------------------------------------
    bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))
    IF ( bnorm <= TINY( bnorm) .AND..NOT.SkipZeroRhs) THEN
      CALL Info('SolveSystem','Solution trivially zero!')
      x = 0.0d0

      ! Increase the nonlinear counter since otherwise some stuff may stagnate
      ! Normally this is done within ComputeChange
      iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
      Solver % Variable % NonlinIter = iterV % Values(1)
      iterV % Values(1) = iterV % Values(1) + 1 
      Solver % Variable % Norm = 0.0_dp
      Solver % Variable % NonlinConverged = 1

      RETURN
    END IF

! 
!   Convert rhs & initial value to the scaled system:
!   -------------------------------------------------

    IF ( ScaleSystem ) THEN
      ApplyRowEquilibration = ListGetLogical(Params,'Linear System Row Equilibration',GotIt)
      IF ( ApplyRowEquilibration ) THEN
        Parallel = ParEnv % PEs > 1
        CALL RowEquilibration(A, b, Parallel)
      ELSE
        CALL ScaleLinearSystem(Solver, A, b, x, RhsScaling=bnorm/=0._dp,ConstraintScaling=.TRUE. )
      END IF
    END IF

    ComputeChangeScaled = ListGetLogical(Params,'Nonlinear System Compute Change in Scaled System',GotIt)
    IF(.NOT.GotIt) ComputeChangeScaled = .FALSE.

    IF(ComputeChangeScaled) THEN
       ALLOCATE(NonlinVals(SIZE(x)))
       NonlinVals = x
       CALL RotateNTSystemAll(NonlinVals, Solver % Variable % Perm, DOFs)
    END IF

    ! Sometimes the r.h.s. may abruptly diminish in value resulting to significant 
    ! convergence issues or it may be that the system scales linearly with the source. 
    ! This flag tries to improve on the initial guess of the linear solvers, and may 
    ! sometimes even result to the exact solution.
    IF( ListGetLogical( Params,'Linear System Normalize Guess',GotIt ) ) THEN
      ALLOCATE( TempVector(A % NumberOfRows) )

      IF ( ParEnv % PEs > 1 ) THEN
        IF( .NOT. ALLOCATED( TempRHS ) ) THEN
          ALLOCATE( TempRHS(A % NumberOfRows) ); TempRHS=0._dp
        END IF

        Tempvector = 0._dp
        TempRHS(1:n) = b(1:n)
        CALL ParallelInitSolve( A, x, TempRHS, Tempvector )

        MP => ParallelMatrix(A,mx,mb,mr)
        mn = MP % NumberOfRows

        TempVector = 0._dp
        CALL ParallelMatrixVector( A, mx, TempVector )

        bn = ParallelDot( mn, TempVector, mb )
        xn = ParallelDot( mn, TempVector, TempVector )
        DEALLOCATE( TempRHS )
      ELSE
        CALL MatrixVectorMultiply( A, x, TempVector )
        xn = SUM( TempVector(1:n)**2 )
        bn = SUM( TempVector(1:n) * b(1:n) )
      END IF

      IF( xn > TINY( xn ) ) THEN
        x(1:n) = x(1:n) * ( bn / xn )
        WRITE( Message,'(A,ES12.3)') 'Linear System Normalizing Factor: ',bn/xn
        CALL Info('SolveLinearSystem',Message,Level=8) 
      END IF
      DEALLOCATE( TempVector )
    END IF

    Method = ListGetString(Params,'Linear System Solver',GotIt)
    IF (Method=='multigrid' .OR. Method=='iterative' ) THEN
      Prec = ListGetString(Params,'Linear System Preconditioning',GotIt)
      IF ( Prec=='vanka' ) CALL VankaCreate(A,Solver)
      IF ( Prec=='circuit' ) CALL CircuitPrecCreate(A,Solver)
    END IF

    IF ( ParEnv % PEs <= 1 ) THEN

      SELECT CASE(Method)
      CASE('multigrid')
        CALL MultiGridSolve( A, x, b, &
            DOFs, Solver, Solver % MultiGridLevel )
      CASE('iterative')
        CALL IterSolver( A, x, b, Solver )
      CASE('feti')
        CALL Fatal('SolveLinearSystem', &
            'Feti solver available only in parallel.')
      CASE('block')
        CALL BlockSolveExt( A, x, b, Solver )
      CASE DEFAULT
        CALL DirectSolver( A, x, b, Solver )
      END SELECT
    ELSE
      SELECT CASE(Method)
      CASE('multigrid')
        CALL MultiGridSolve( A, x, b, &
            DOFs, Solver, Solver % MultiGridLevel )
      CASE('iterative')
        CALL ParallelIter( A, A % ParallelInfo, DOFs, &
            x, b, Solver, A % ParMatrix )
      CASE('feti')
        CALL FetiSolver( A, x, b, Solver )
      CASE('block')
        CALL BlockSolveExt( A, x, b, Solver )
     CASE DEFAULT
        CALL DirectSolver( A, x, b, Solver )
      END SELECT
    END IF

    IF(ComputeChangeScaled) THEN
      CALL ComputeChange(Solver,.FALSE.,n, x, NonlinVals, Matrix=A, RHS=b )
      DEALLOCATE(NonlinVals)
    END IF

    IF ( ScaleSystem ) THEN
      IF ( ApplyRowEquilibration ) THEN
        CALL ReverseRowEquilibration( A, b )
      ELSE
        CALL BackScaleLinearSystem( Solver, A, b, x, ConstraintScaling=.TRUE. )
      END IF
    END IF

    Aaid => A
    IF (PRESENT(BulkMatrix)) THEN
      IF (ASSOCIATED(BulkMatrix) ) Aaid=>BulkMatrix
    END IF
    
    NodalLoads => VariableGet( Solver % Mesh % Variables, &
        GetVarName(Solver % Variable) // ' Loads' )
    IF( ASSOCIATED( NodalLoads ) ) THEN
      CALL Info('SolveLinearSystem','Calculating nodal loads',Level=6)
      CALL CalculateLoads( Solver, Aaid, x, Dofs, .TRUE., NodalLoads ) 
    END IF

    IF (BackRotation) THEN
      CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )
      IF( ASSOCIATED( NodalLoads ) ) THEN
        CALL BackRotateNTSystem(NodalLoads % Values,NodalLoads % Perm,DOFs)
      END IF
    END IF

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Compute the change of the solution with different methods 
!------------------------------------------------------------------------------
    IF(.NOT.ComputeChangeScaled) THEN
      CALL ComputeChange(Solver,.FALSE.,n, x, Matrix=A, RHS=b )
    END IF
    Norm = Solver % Variable % Norm

!------------------------------------------------------------------------------
 
   Solver % Variable % PrimaryMesh => Solver % Mesh
   CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
         GetVarName(Solver % Variable) )
   
   IF ( ASSOCIATED( NodalLoads ) ) THEN
     NodalLoads % PrimaryMesh => Solver % Mesh
     CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
                  GetVarName(NodalLoads) )
   END IF

!------------------------------------------------------------------------------
! In order to be able to change the preconditoners or solvers the old matrix structures
! must be deallocated on request.

    IF( ListGetLogical( Params, 'Linear System Preconditioning Deallocate', GotIt) ) THEN
       ! ILU preconditioning
       IF( ASSOCIATED(A % ILUValues) ) THEN
          IF(  SIZE( A % ILUValues) /= SIZE(A % Values) ) &
             DEALLOCATE(A % ILUCols, A % ILURows, A % ILUDiag)
          DEALLOCATE(A % ILUValues)
       END IF
          
       ! Multigrid solver / preconditioner
       IF( Solver % MultigridLevel > 0 ) THEN
          Aaid => A 
          IF(ASSOCIATED( Aaid % Parent) ) THEN
             DO WHILE( ASSOCIATED( Aaid % Parent ) )
                Aaid => Aaid % Parent
             END DO
             DO WHILE( ASSOCIATED( Aaid % Child) )
                Aaid => Aaid % Child
                IF(ASSOCIATED(Aaid % Parent)) DEALLOCATE(Aaid % Parent )
                IF(ASSOCIATED(Aaid % Ematrix)) DEALLOCATE(Aaid % Ematrix )
             END DO
          END IF
       END IF
    END IF

  END SUBROUTINE SolveLinearSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Given a linear system Ax=b make a change of variables such that we will 
!> be solving for the residual Adx=b-Ax0 where dx=x-x0.
!------------------------------------------------------------------------------
  SUBROUTINE LinearSystemResidual( A, b, x, r )

    REAL(KIND=dp) CONTIG :: b(:)   
    REAL(KIND=dp) CONTIG :: x(:)   
    TYPE(Matrix_t), POINTER :: A   
    REAL(KIND=dp), POINTER :: r(:)
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: TmpXVec, TmpRVec, TmpRHSVec

    INTEGER :: i,n,nn 

    n = A % NumberOfRows 

    IF (Parenv % Pes>1) THEN
      CALL ParallelInitSolve(A,x,b,r)
      CALL ParallelMatrixVector(A,x,r,.TRUE.)
    ELSE
      CALL MatrixVectorMultiply( A, x, r)
    END IF

    DO i=1,n
      r(i) = b(i) - r(i)
    END DO

  END SUBROUTINE LinearSystemResidual


!------------------------------------------------------------------------------
!> Solve a system. Various additional utilities are included and 
!> naturally a call to the linear system solver.
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE SolveSystem( A,ParA,b,x,Norm,DOFs,Solver )
!------------------------------------------------------------------------------
    REAL(KIND=dp) CONTIG, TARGET :: b(:)   !< The RHS vector
    REAL(KIND=dp) CONTIG :: x(:)   !< Previous solution on entry, new solution on exit (hopefully)
    REAL(KIND=dp) :: Norm          !< L2 Norm of solution
    TYPE(Matrix_t), POINTER :: A   !< The coefficient matrix
    INTEGER :: DOFs                !< Number of degrees of freedom per node for this equation
    TYPE(Solver_t), TARGET :: Solver                 !< Holds various solver options.
    TYPE(SParIterSolverGlobalD_t), POINTER :: ParA   !< holds info for parallel solver, 
                                                     !< if not executing in parallel this is just a dummy.
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var, NodalLoads
    TYPE(Mesh_t), POINTER :: Mesh, SaveMEsh
    LOGICAL :: Relax, Found, NeedPrevSol, Timing, ResidualMode
    INTEGER :: n,i,j,k,l,m,istat,nrows,ncols,ConstrainedSolve,colsj,rowoffset
    CHARACTER(LEN=MAX_NAME_LEN) :: Method, ProcName, VariableName
    INTEGER(KIND=AddrInt) :: Proc
    REAL(KIND=dp) :: Relaxation,Beta,Gamma
    REAL(KIND=dp), ALLOCATABLE :: Diag(:), TempVector(:)
    REAL(KIND=dp), POINTER :: bb(:),Res(:)
    REAL(KIND=dp) :: t0,rt0,rst,st,ct
#ifndef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: CPUTime,RealTime
#endif
    TYPE(ValueList_t), POINTER :: Params

#ifndef USE_ISO_C_BINDINGS
    INTERFACE ExecLinSolveProcs
      INTEGER FUNCTION ExecLinSolveProcs( Proc,Model,Solver,A,b,x,n,DOFs,Norm )
        USE Types
        INTEGER(KIND=AddrInt) :: Proc
        TYPE(Model_t) :: Model
        TYPE(Solver_t) :: Solver
        TYPE(Matrix_t), POINTER :: A
        INTEGER :: n, DOFs
        REAL(KIND=dp) :: x(n),b(n), Norm
      END FUNCTION ExecLinSolveProcs
    END INTERFACE
#endif

!------------------------------------------------------------------------------
    Params => Solver % Values

    CALL Info('SolveSystem','Solving linear system',Level=10)

    Timing = ListCheckPrefix(Params,'Linear System Timing')
    IF( Timing ) THEN
      t0 = CPUTime(); rt0 = RealTime()
    END IF

    n = A % NumberOfRows

    ResidualMode = ListGetLogical( Params,'Linear System Residual Mode',Found )

!------------------------------------------------------------------------------
! The allocation of previous values has to be here in order to 
! work properly with the Dirichlet elimination.
!------------------------------------------------------------------------------
    NeedPrevSol = ResidualMode

    IF(.NOT. NeedPrevSol ) THEN
      Relaxation = ListGetConstReal( Params, &
          'Nonlinear System Relaxation Factor', Found )
      IF( Found ) NeedPrevSol = (Relaxation /= 1.0_dp)
    END IF

    IF(.NOT. NeedPrevSol ) THEN
      Method = ListGetString( Params, &
        'Nonlinear System Convergence Measure', Found ) 
      NeedPrevSol = ( Method == 'residual' .OR. Method == 'solution' )
    END IF

    IF( NeedPrevSol ) THEN
      CALL Info('SolveSystem','Previous solution must be stored before system is solved',Level=10)
      Found = ASSOCIATED(Solver % Variable % NonlinValues)
      IF( Found ) THEN
        IF ( SIZE(Solver % Variable % NonlinValues) /= n) THEN
          DEALLOCATE(Solver % Variable % NonlinValues)
          Found = .FALSE.
        END IF
      END IF
      IF(.NOT. Found) THEN
        ALLOCATE( Solver % Variable % NonlinValues(n), STAT=istat ) 
        IF ( istat /= 0 ) CALL Fatal( 'SolveSystem', 'Memory allocation error.' )
      END IF
      Solver % Variable % NonlinValues = x(1:n)
    END IF

    IF ( Solver % LinBeforeProc /= 0 ) THEN
       istat = ExecLinSolveProcs( Solver % LinBeforeProc,CurrentModel,Solver, &
                       A, b, x, n, DOFs, Norm )
       IF ( istat /= 0 ) GOTO 10
    END IF

    ! If residual mode is requested make change of variables:
    ! Ax=b -> Adx = b-Ax0 = r
    IF( ResidualMode ) THEN
      CALL Info('SolveSystem','Changing the equation to residual based mode',Level=10)
      ALLOCATE( Res(n) ) 

      ! If needed move the current solution to N-T coordinate system
      ! before computing the residual.     
      CALL RotateNTSystemAll(x, Solver % Variable % Perm, DOFs)

      CALL LinearSystemResidual( A, b, x, res )
      bb => res
      ! Set the initial guess for the redidual system to zero
      x = 0.0_dp
    ELSE
      bb => b
    END IF

    ConstrainedSolve = 0
    IF ( ASSOCIATED(A % ConstraintMatrix) )  THEN
      IF ( A % ConstraintMatrix % NumberOFRows >= 1 ) & 
        ConstrainedSolve = 1
    END IF

    IF ( ASSOCIATED(A % AddMatrix) )  THEN
      IF ( A % AddMatrix % NumberOFRows >= 1 ) ConstrainedSolve = 1
    END IF

    ConstrainedSolve = ParallelReduction(ConstrainedSolve*1._dp)
   
    IF ( ConstrainedSolve > 0 ) THEN
      CALL Info('SolveSystem','Solving linear system with constraint matrix',Level=10)
      IF( ListGetLogical( Params,'Save Constraint Matrix',Found ) ) THEN
        CALL SaveProjector(A % ConstraintMatrix,.TRUE.,'cm')
      END IF
      CALL SolveWithLinearRestriction( A,bb,x,Norm,DOFs,Solver )
    ELSE
      CALL Info('SolveSystem','Solving linear system without constraint matrix',Level=12)
      CALL SolveLinearSystem( A,bb,x,Norm,DOFs,Solver )
    END IF
    CALL Info('SolveSystem','Linear system solved',Level=12)

    ! Even in the residual mode the system is reverted back to complete vectors 
    ! and we may forget about the residual.
    IF( ResidualMode ) DEALLOCATE( Res ) 

!------------------------------------------------------------------------------

10  CONTINUE

    IF ( Solver % LinAfterProc /= 0 ) THEN
      istat = ExecLinSolveProcs( Solver % LinAfterProc, CurrentModel, Solver, &
              A, b, x, n, DOFs, Norm )
    END IF

    IF ( Solver % TimeOrder == 2 ) THEN
      CALL Info('SolveSystem','Setting up PrevValues for 2nd order transient equations',Level=12)

      IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        Gamma =  0.5d0 - Solver % Alpha
        Beta  = (1.0d0 - Solver % Alpha)**2 / 4.0d0
        DO i=1,n
          Solver % Variable % PrevValues(i,2) = &
             (1.0d0/(Beta*Solver % dt**2))* &
               (x(i)-Solver % Variable % PrevValues(i,3)) -  &
                  (1.0d0/(Beta*Solver % dt))*Solver % Variable % PrevValues(i,4)+ &
                        (1.0d0-1.0d0/(2*Beta))*Solver % Variable % PrevValues(i,5)

          Solver % Variable % PrevValues(i,1) = &
            Solver % Variable % PrevValues(i,4) + &
               Solver % dt*((1.0d0-Gamma)*Solver % Variable % PrevValues(i,5)+&
                  Gamma*Solver % Variable % PrevValues(i,2))
        END DO
      END IF
    END IF

    IF( Timing ) THEN
      st  = CPUTime() - t0;
      rst = RealTime() - rt0

      WRITE(Message,'(a,f8.2,f8.2,a)') 'Linear system time (CPU,REAL) for '&
          //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
      CALL Info('SolveSystem',Message)    
      
      IF( ListGetLogical(Params,'Linear System Timing',Found)) THEN
        CALL ListAddConstReal(CurrentModel % Simulation,'res: linsys cpu time '&
            //GetVarName(Solver % Variable),st)
        CALL ListAddConstReal(CurrentModel % Simulation,'res: linsys real time '&
            //GetVarName(Solver % Variable),rst)
      END IF
      
      IF( ListGetLogical(Params,'Linear System Timing Cumulative',Found)) THEN
        ct = ListGetConstReal(CurrentModel % Simulation,'res: cum linsys cpu time '&
            //GetVarName(Solver % Variable),Found)
        st = st + ct
        ct = ListGetConstReal(CurrentModel % Simulation,'res: cum linsys real time '&
            //GetVarName(Solver % Variable),Found)
        rst = rst + ct
        CALL ListAddConstReal(CurrentModel % Simulation,'res: cum linsys cpu time '&
            //GetVarName(Solver % Variable),st)
        CALL ListAddConstReal(CurrentModel % Simulation,'res: cum linsys real time '&
            //GetVarName(Solver % Variable),rst)
      END IF

    END IF

    CALL Info('SolveSystem','Finished solving the system',Level=12)

!------------------------------------------------------------------------------
END SUBROUTINE SolveSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solve a linear eigen system.
!------------------------------------------------------------------------------
SUBROUTINE SolveEigenSystem( StiffMatrix, NOFEigen, &
        EigenValues, EigenVectors,Solver )
!------------------------------------------------------------------------------
    USE EigenSolve
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: EigenValues(:),EigenVectors(:,:)
    REAL(KIND=dp) :: Norm
    TYPE(Matrix_t), POINTER :: StiffMatrix
    INTEGER :: NOFEigen
    TYPE(Solver_t) :: Solver
    !------------------------------------------------------------------------------
    INTEGER :: n
    !------------------------------------------------------------------------------
    n = StiffMatrix % NumberOfRows

    IF ( .NOT. Solver % Matrix % COMPLEX ) THEN
      IF ( ParEnv % PEs <= 1 ) THEN
        CALL ArpackEigenSolve( Solver, StiffMatrix, n, NOFEigen, &
                EigenValues, EigenVectors )
      ELSE
        CALL ParallelArpackEigenSolve( Solver, StiffMatrix, n, NOFEigen, &
                EigenValues, EigenVectors )
      END IF
    ELSE
      IF ( ParEnv % PEs <= 1 ) THEN
        CALL ArpackEigenSolveComplex( Solver, StiffMatrix, n/2, &
              NOFEigen, EigenValues, EigenVectors )
      ELSE
        CALL ParallelArpackEigenSolveComplex( Solver, StiffMatrix, n/2, NOFEigen, &
                EigenValues, EigenVectors )
      END IF
    END IF

!------------------------------------------------------------------------------
END SUBROUTINE SolveEigenSystem
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Solve a linear system with permutated constraints.
!------------------------------------------------------------------------------
SUBROUTINE SolveConstraintModesSystem( StiffMatrix, Solver )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: StiffMatrix
    TYPE(Solver_t) :: Solver
    TYPE(Variable_t), POINTER :: Var
    !------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,m
    LOGICAL :: PrecRecompute, Stat, Found, ComputeFluxes, Symmetric, ScaleSystem
    REAL(KIND=dp), POINTER :: PValues(:)
    REAL(KIND=dp), ALLOCATABLE :: Fluxes(:), FluxesMatrix(:,:)
    !------------------------------------------------------------------------------
    n = StiffMatrix % NumberOfRows

    Var => Solver % Variable
    IF( SIZE(Var % Values) /= n ) THEN
      CALL Fatal('SolveConstraintModesSystem','Conflicting sizes for matrix and variable!')
    END IF

    m = Var % NumberOfConstraintModes
    IF( m == 0 ) THEN
      CALL Fatal('SolveConstraintModesSystem','No constraint modes?!')
    END IF

    ComputeFluxes = ListGetLogical( Solver % Values,'Constraint Modes Fluxes',Found) 
    IF( ComputeFluxes ) THEN
      CALL Info('SolveConstraintModesSystem','Allocating for lumped fluxes',Level=10)
      ALLOCATE( Fluxes( n ) )
      ALLOCATE( FluxesMatrix( m, m ) )
      FluxesMatrix = 0.0_dp
    END IF

    DO i=1,m
      CALL Info('SolveConstraintModesSystem','Solving for mode: '//TRIM(I2S(i)))

      IF( i == 2 ) THEN
        CALL ListAddLogical( Solver % Values,'No Precondition Recompute',.TRUE.)
      END IF

      WHERE( Var % ConstraintModesIndeces == i ) StiffMatrix % Rhs = 1.0_dp      

      CALL SolveSystem( StiffMatrix,ParMatrix,StiffMatrix % rhs,&
          Var % Values,Var % Norm,Var % DOFs,Solver )

      WHERE( Var % ConstraintModesIndeces == i ) StiffMatrix % Rhs = 0.0_dp            

      Var % ConstraintModes(i,:) = Var % Values

      IF( ComputeFluxes ) THEN
        CALL Info('SolveConstraintModesSystem','Computing lumped fluxes',Level=8)
        PValues => StiffMatrix % Values
        StiffMatrix % Values => StiffMatrix % BulkValues
        Fluxes = 0.0_dp
        CALL MatrixVectorMultiply( StiffMatrix, Var % Values, Fluxes ) 
        StiffMatrix % Values => PValues

        DO j=1,n
          k = Var % ConstraintModesIndeces(j)
          IF( k > 0 ) THEN
            IF(.FALSE.) THEN
              FluxesMatrix(i,k) = FluxesMatrix(i,k) + Fluxes(j)
            ELSE
              IF( i /= k ) THEN
                FluxesMatrix(i,k) = FluxesMatrix(i,k) - Fluxes(j)
              END IF
              FluxesMatrix(i,i) = FluxesMatrix(i,i) + Fluxes(j)
            END IF
          END IF
        END DO
      END IF
    END DO

    IF( ComputeFluxes ) THEN
      Symmetric = ListGetLogical( Solver % Values,&
          'Constraint Modes Fluxes Symmetric', Found ) 
      IF( Symmetric ) THEN
        FluxesMatrix = 0.5_dp * ( FluxesMatrix + TRANSPOSE( FluxesMatrix ) )
      END IF
      DO i=1,m
        DO j=1,m
          IF( Symmetric .AND. j < i ) CYCLE
          WRITE( Message, '(I3,I3,ES15.5)' ) i,j,FluxesMatrix(i,j)
          CALL Info( 'SolveConstraintModesSystem', Message, Level=4 )
        END DO
      END DO
      DEALLOCATE( Fluxes )
    END IF


    CALL ListAddLogical( Solver % Values,'No Precondition Recompute',.FALSE.)
    
!------------------------------------------------------------------------------
  END SUBROUTINE SolveConstraintModesSystem
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> A parser of the variable name that returns the true variablename
!> where the inline options have been interpreted.
!------------------------------------------------------------------------------
SUBROUTINE VariableNameParser(var_name, NoOutput, Global, Dofs )

  CHARACTER(LEN=*)  :: var_name
  LOGICAL, OPTIONAL :: NoOutput, Global
  INTEGER, OPTIONAL :: Dofs

  INTEGER :: i,j,k

  IF(PRESENT(NoOutput)) NoOutput = .FALSE.
  IF(PRESENT(Global)) Global = .FALSE.
  IF(PRESENT(Dofs)) Dofs = 0

  DO WHILE( var_name(1:1) == '-' )
    IF ( SEQL(var_name, '-nooutput ') ) THEN
      IF(PRESENT(NoOutput)) NoOutput = .TRUE.
      var_name(1:LEN(var_name)-10) = var_name(11:)
    END IF
    
    IF ( SEQL(var_name, '-global ') ) THEN
      IF(PRESENT(Global)) Global = .TRUE.
      var_name(1:LEN(var_name)-8) = var_name(9:)
    END IF
    
    IF ( SEQL(var_name, '-dofs ') ) THEN
      IF(PRESENT(DOFs)) READ( var_name(7:), * ) DOFs     
      j = LEN_TRIM( var_name )
      k = 7
      DO WHILE( var_name(k:k) /= ' '  )
        k = k + 1
        IF ( k > j ) EXIT
      END DO
      var_name(1:LEN(var_name)-(k+2)) = var_name(k+1:)
    END IF
  END DO

END SUBROUTINE VariableNameParser


!------------------------------------------------------------------------------
!> Updates values for exported variables which are typically auxiliary variables derived
!> from the solution.
!------------------------------------------------------------------------------
  SUBROUTINE UpdateExportedVariables( Solver )  
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  
  INTEGER :: i,j,k,l,n,m,t,bf_id,dofs,nsize
  CHARACTER(LEN=MAX_NAME_LEN) :: str, var_name,tmpname,condname
  REAL(KIND=dp), POINTER :: Values(:), Solution(:), LocalSol(:), LocalCond(:)
  INTEGER, POINTER :: Indexes(:), Perm(:)
  LOGICAL :: Found, AllocationsDone, Conditional, GotIt, StateVariable
  LOGICAL, POINTER :: ActivePart(:),ActiveCond(:)
  TYPE(Variable_t), POINTER :: ExpVariable
  TYPE(ValueList_t), POINTER :: ValueList
  TYPE(Element_t),POINTER :: Element

  SAVE AllocationsDone

  AllocationsDone = .FALSE.

  l = 0
  DO WHILE( .TRUE. )
    l = l + 1
    str = ComponentName( 'exported variable', l )
    var_name = ListGetString( Solver % Values, str, GotIt )
    
    IF(.NOT. GotIt) EXIT
    
    CALL VariableNameParser( var_name ) 

    ExpVariable => VariableGet( Solver % Mesh % Variables, Var_name )
    IF( .NOT. ASSOCIATED(ExpVariable)) CYCLE
    
    WRITE(Message,*) 'Trying to set values for variable: '//TRIM(Var_name)
    CALL Info('UpdateExportedVariables',Message,Level=6)
  
    IF( .NOT. AllocationsDone) THEN      
      m = CurrentModel % NumberOFBodyForces
      ALLOCATE( ActivePart(m), ActiveCond(m) )

      m = Solver % Mesh % MaxElementDOFs
      ALLOCATE( LocalSol(m), LocalCond(m))

      AllocationsDone = .TRUE.
    END IF

    Dofs = ExpVariable % DOFs
    Values => ExpVariable % Values
    Perm => ExpVariable % Perm
    n = LEN_TRIM( var_name )

    StateVariable = ( SIZE( Values ) == DOFs )
    IF( StateVariable ) THEN
      IF( Dofs > 1 ) THEN
        tmpname = ComponentName( var_name(1:n), j )
        Solution => Values( j:j )
      ELSE
        tmpname = var_name(1:n)
        Solution => Values
      END IF
 
      DO bf_id=1,CurrentModel % NumberOFBodyForces
        IF( ListCheckPresent( &
            CurrentModel % BodyForces(bf_id) % Values,TmpName ) ) THEN
          CALL Info('UpdateExportedVariables',&
              'Found a proper definition for state variable',Level=6)
          Solution = ListGetCReal( CurrentModel % BodyForces(bf_id) % Values,TmpName)
          EXIT
        END IF
      END DO
      CYCLE
    END IF	

    
    DO j=1,DOFs
      
      IF( Dofs > 1 ) THEN
        tmpname = ComponentName( var_name(1:n), j )
        nSize = DOFs * SIZE(Solver % Variable % Values) / Solver % Variable % DOFs
        Perm => Solver % Variable % Perm
        Solution => Values( j:nSize-DOFs+j:DOFs )
      ELSE
        tmpname = var_name(1:n)
        Solution => Values
      END IF
      condname = TRIM(tmpname) //' Condition' 
        
      !------------------------------------------------------------------------------
      ! Go through the Dirichlet conditions in the body force lists
      !------------------------------------------------------------------------------
      
      ActivePart = .FALSE.
      ActiveCond = .FALSE.

      DO bf_id=1,CurrentModel % NumberOFBodyForces
        ActivePart(bf_id) = ListCheckPresent( &
            CurrentModel % BodyForces(bf_id) % Values,TmpName ) 
        ActiveCond(bf_id) = ListCheckPresent( &
            CurrentModel % BodyForces(bf_id) % Values,CondName )      
      END DO
      
      IF ( .NOT. ANY( ActivePart ) ) CYCLE

      CALL Info('UpdateExportedVariables','Found a proper definition in body forces',Level=6)

      DO t = 1, Solver % NumberOfActiveElements 
        Element => CurrentModel % Elements(Solver % ActiveElements(t) )
        bf_id = ListGetInteger( CurrentModel % Bodies(Element % BodyId) % Values,&
            'Body Force',GotIt)
        
        IF(.NOT. GotIt) CYCLE
        IF(.NOT. ActivePart(bf_id)) CYCLE
        Conditional = ActiveCond(bf_id)
        
        CurrentModel % CurrentElement => Element
        m = Element % TYPE % NumberOfNodes
        Indexes => Element % NodeIndexes
        ValueList => CurrentModel % BodyForces(bf_id) % Values
        
        LocalSol(1:m) = ListGetReal(ValueList, TmpName, m, Indexes(1:m) )
        IF( Conditional ) THEN
          LocalCond(1:m) = ListGetReal(ValueList, CondName, m, Indexes(1:m) )
          DO i=1,m
            IF( LocalCond(i) > 0.0_dp ) THEN
              Solution( Perm(Indexes(i)) ) = LocalSol(i)
            END IF
          END DO
        ELSE
          Solution( Perm(Indexes(1:m)) ) = LocalSol(1:m)
        END IF
      END DO
        
    END DO
  END DO

  IF( AllocationsDone ) THEN
    DEALLOCATE(ActivePart, ActiveCond, LocalSol, LocalCond)
    AllocationsDone = .FALSE.
  END IF

END SUBROUTINE UpdateExportedVariables


!------------------------------------------------------------------------------
!> Eliminates bubble degrees of freedom from a local linear system.
!> This version is suitable for flow models with velocity and pressure as 
!> unknowns.
!------------------------------------------------------------------------------
SUBROUTINE NSCondensate( N, Nb, dim, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb, dim
    REAL(KIND=dp) :: K(:,:),F(:),F1(:), Kbb(nb*dim,nb*dim), &
      Kbl(nb*dim,n*(dim+1)),Klb(n*(dim+1),nb*dim),Fb(nb*dim)

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

    Fb  = F1(Bdofs)
    F1(1:(dim+1)*n) = F1(1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
!------------------------------------------------------------------------------
END SUBROUTINE NSCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE Condensate( N, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N
    REAL(KIND=dp) :: K(:,:),F(:)
    REAL(KIND=dp), OPTIONAL :: F1(:)
!------------------------------------------------------------------------------    
    REAL(KIND=dp) :: Kbb(N,N), &
        Kbl(N,N),Klb(N,N),Fb(N)
    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(N)

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = Ldofs + n

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,n )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )

    IF( PRESENT( F1 ) ) THEN
      Fb  = F1(Bdofs)
      F1(1:n) = F1(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    END IF
!------------------------------------------------------------------------------
END SUBROUTINE Condensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>     Subroutine for condensation of p element bubbles from linear problem.
!>     Modifies given stiffness matrix and force vector(s) 
!------------------------------------------------------------------------------
SUBROUTINE CondensateP( N, Nb, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N               !< Sum of nodal, edge and face degrees of freedom.
    INTEGER :: Nb              !< Sum of internal (bubble) degrees of freedom.
    REAL(KIND=dp) :: K(:,:)    !< Local stiffness matrix.
    REAL(KIND=dp) :: F(:)      !< Local force vector.
    REAL(KIND=dp), OPTIONAL :: F1(:)  !< Local second force vector.
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Kbb(Nb,Nb), &
    Kbl(Nb,N), Klb(N,Nb), Fb(Nb)
    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    IF (PRESENT(F1)) THEN
      F1(1:n) = F1(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    END IF

    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
END SUBROUTINE CondensateP
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a harmonic system.
!------------------------------------------------------------------------------
SUBROUTINE SolveHarmonicSystem( G, Solver )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), TARGET :: G
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: BMatrix, A => NULL()
    INTEGER :: Rounds = 1000,i,j,k,n, ILUn, kr, ki, DOFs, ne, niter
    LOGICAL :: stat, Found, OptimizeBW, DirectLinearSolver,Real_given,Imag_given
    CHARACTER(LEN=MAX_NAME_LEN) :: Name
    REAL(KIND=dp) :: Omega = 10, norm, TOL=1.0d-6, s, ILUTol
    REAL(KIND=dp), POINTER :: Freqv(:,:)
    REAL(KIND=dp), ALLOCATABLE :: x(:), b(:)
    REAL(KIND=dp) :: frequency
    INTEGER :: Nfrequency
    TYPE(ValueList_t), POINTER :: BC

    DO j=1,Solver % Variable % DOFs
      Name = ComponentName( Solver % Variable % Name, j ) 
      DO i=1,CurrentModel % NumberOFBCs
        BC => CurrentModel % BCs(i) % Values
        real_given = ListCheckPresent( BC, Name )
        imag_given = ListCheckPresent( BC, TRIM(Name) // ' im' )

        IF ( real_given .AND. .NOT. imag_given ) THEN
            CALL ListAddConstReal( BC, TRIM(Name) // ' im', 0._dp)
        ELSE IF ( imag_given .AND. .NOT. real_given ) THEN
            CALL ListAddConstReal( BC, Name, 0._dp )
        END IF
      END DO
    END DO

    n = Solver % Matrix % NumberofRows
    DOFs = Solver % Variable % DOFs * 2

    OptimizeBW = ListGetLogical(Solver % Values, 'Optimize Bandwidth', Found)
    IF ( .NOT. Found ) OptimizeBW = .TRUE.

    A => G
    DO WHILE( ASSOCIATED(A) )
      BMatrix => A
      A => A % EMatrix
      IF ( ASSOCIATED(A) ) THEN
        IF ( A % COMPLEX ) EXIT
      END IF
    END DO

    IF ( .NOT. ASSOCIATED(A) ) THEN
      A => CreateMatrix( CurrentModel, Solver, Solver % Mesh,   &
              Solver % Variable % Perm, DOFs, MATRIX_CRS, OptimizeBW, &
              ListGetString( Solver % Values, 'Equation') )
      A % COMPLEX = .TRUE.
      BMatrix % EMatrix => A
    END IF

    ALLOCATE( x(2*n), b(2*n) )
    x = 0
    b(1:2*n:2) = G % RHS(1:n)
    b(2:2*n:2) = G % RHS_im(1:n)


    Nfrequency = ListGetInteger( Solver % Values,'Harmonic System Values',Found )
    IF( Nfrequency > 1 ) THEN
      freqv => ListGetConstRealArray( Solver % Values, 'Frequency' )
    ELSE
      Frequency = ListGetAngularFrequency( Solver % Values, Found ) / (2*PI)
      IF( .NOT. Found ) THEN
        CALL Fatal( 'AddEquation', '> Frequency < must be given for harmonic analysis.' )
      END IF
      Nfrequency = 1
    END IF

    niter = MIN(Nfrequency,Solver % NOFEigenValues)
    ne=Solver % NofEigenValues
    Solver % NofEigenValues=0

    DO i=1,niter
      IF( Nfrequency > 1 ) THEN
        Frequency = freqv(i,1)
        WRITE( Message, '(a,i5,e12.3)' ) 'Frequency sweep: ', i, frequency
      ELSE
        WRITE( Message, '(a,e12.3)' ) 'Frequency value: ', frequency
      END IF
      CALL Info( 'HarmonicSolve', ' ' )
      CALL Info( 'HarmonicSolve', Message )

      omega = 2 * PI * Frequency
      DO k=1,n
        kr = A % Rows(2*(k-1)+1)
        ki = A % Rows(2*(k-1)+2)
        DO j=G % Rows(k),G % Rows(k+1)-1
          A % Values(kr)   =  G % Values(j)
          IF (ASSOCIATED(G % MassValues)) A % Values(kr) = &
              A % Values(kr) - omega**2*G % MassValues(j)
          IF (ASSOCIATED(G % DampValues)) THEN
            A % Values(kr+1) = -G % Dampvalues(j) * omega
            A % Values(ki)   =  G % Dampvalues(j) * omega
          END IF
          A % Values(ki+1) =  G % Values(j)
          IF (ASSOCIATED(G % MassValues)) A % Values(ki+1) = &
            A % Values(ki+1) - omega**2*G % MassValues(j)
          kr = kr + 2
          ki = ki + 2
        END DO
      END DO

      DO j=1,Solver % Variable % DOFs
        Name = ComponentName( Solver % Variable % Name, j ) 

        CALL SetDirichletBoundaries( CurrentModel, A, b, Name, &
                2*j-1, DOFs, Solver % Variable % Perm )

        CALL SetDirichletBoundaries( CurrentModel, A, b, TRIM(Name) // ' im', &
                2*j, DOFs, Solver % Variable % Perm )
      END DO

      CALL SolveLinearSystem( A, b, x, Norm, DOFs, Solver )

      DO j=1,n
        Solver % Variable % EigenVectors(i,j) = &
                 CMPLX( x(2*(j-1)+1),x(2*(j-1)+2),KIND=dp )
      END DO
    END DO

    Solver % NOFEigenValues = ne

    DEALLOCATE( x, b )
!------------------------------------------------------------------------------
 END SUBROUTINE SolveHarmonicSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  This subroutine will solve the system with some linear restriction.
!>  The restriction matrix is assumed to be in the ConstraintMatrix-field of 
!>  the StiffMatrix. The restriction vector is the RHS-field of the
!>  ConstraintMatrix.
!>  NOTE: Only serial solver implemented so far ...
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE SolveWithLinearRestriction( StiffMatrix, ForceVector, Solution, &
        Norm, DOFs, Solver )
!------------------------------------------------------------------------------  
  IMPLICIT NONE
  TYPE(Matrix_t), POINTER :: StiffMatrix !< Linear equation matrix information. 
                                         !< The restriction matrix is assumed to be in the EMatrix-field
  REAL(KIND=dp),TARGET :: ForceVector(:)        !< The right hand side of the linear equation
  REAL(KIND=dp),TARGET :: Solution(:)           !< Previous solution as input, new solution as output.
  REAL(KIND=dp) :: Norm                  !< The L2 norm of the solution.
  INTEGER :: DOFs                        !< Number of degrees of freedon of the equation.
  TYPE(Solver_t), TARGET :: Solver       !< Linear equation solver options.
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: SolverPointer
  TYPE(Matrix_t), POINTER :: CollectionMatrix, RestMatrix, AddMatrix, &
       RestMatrixTranspose, TMat, XMat
  REAL(KIND=dp), POINTER CONTIG :: CollectionVector(:), RestVector(:),&
     MultiplierValues(:), AddVector(:), Tvals(:), Vals(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: CollectionSolution(:), TotValues(:)
  INTEGER :: NumberOfRows, NumberOfValues, MultiplierDOFs, istat, NoEmptyRows 
  INTEGER :: i, j, k, l, m, n, p,q, ix, Loop
  TYPE(Variable_t), POINTER :: MultVar
  REAL(KIND=dp) :: scl, rowsum
  LOGICAL :: Found, ExportMultiplier, NotExplicit, Refactorize, EnforceDirichlet, EliminateDiscont, &
              EmptyRow, ComplexSystem, ConstraintScaling, UseTranspose, EliminateConstraints, &
              SkipConstraints
  SAVE MultiplierValues, SolverPointer

  CHARACTER(LEN=MAX_NAME_LEN) :: MultiplierName
  TYPE(ListMatrix_t), POINTER :: cList
  TYPE(ListMatrixEntry_t), POINTER :: cPtr, cPrev, cTmp

  INTEGER, ALLOCATABLE, TARGET :: SlavePerm(:), SlaveIPerm(:), MasterPerm(:), MasterIPerm(:)
  INTEGER, POINTER :: UsePerm(:), UseIPerm(:)
  REAL(KIND=dp), POINTER :: UseDiag(:)
  TYPE(ListMatrix_t), POINTER :: Lmat(:)
  LOGICAL  :: EliminateFromMaster, EliminateSlave, Parallel
  REAL(KIND=dp), ALLOCATABLE, TARGET :: SlaveDiag(:), MasterDiag(:), DiagDiag(:)

!------------------------------------------------------------------------------
  CALL Info( 'SolveWithLinearRestriction ', ' ', Level=5 )
  SolverPointer => Solver

  Parallel = (ParEnv % PEs > 1 )

  NotExplicit = ListGetLogical(Solver % Values,'No Explicit Constrained Matrix',Found)
  IF(.NOT. Found) NotExplicit=.FALSE.

  RestMatrix => NULL()
  IF(.NOT.NotExplicit) &
        RestMatrix => StiffMatrix % ConstraintMatrix
  RestVector => Null()
  IF(ASSOCIATED(RestMatrix)) RestVector => RestMatrix % RHS

  AddMatrix => StiffMatrix % AddMatrix
  AddVector => NULL()
  IF(ASSOCIATED(AddMatrix)) &
    AddVector => AddMatrix % RHS

  NumberOfRows = StiffMatrix % NumberOfRows

  CollectionMatrix => StiffMatrix % CollectionMatrix
  Refactorize = ListGetLogical(Solver % Values,'Linear System Refactorize',Found)
  IF(.NOT.Found) Refactorize = .TRUE.

  IF(Refactorize.AND..NOT.NotExplicit) THEN
    IF(ASSOCIATED(CollectionMatrix)) THEN
      CALL FreeMatrix(CollectionMatrix)
      CollectionMatrix => NULL()
    END IF
  END IF

  IF(.NOT.ASSOCIATED(CollectionMatrix)) THEN
    CollectionMatrix => AllocateMatrix()
    CollectionMatrix % FORMAT = MATRIX_LIST
  ELSE
    DEALLOCATE(CollectionMatrix % RHS)
    CollectionMatrix % Values = 0.0_dp
  END IF
  IF(NotExplicit) CollectionMatrix % ConstraintMatrix => StiffMatrix % ConstraintMatrix

  NumberOfRows = StiffMatrix % NumberOfRows
  IF(ASSOCIATED(AddMatrix)) NumberOfRows = MAX(NumberOfRows,AddMatrix % NumberOfRows)
  EliminateConstraints = ListGetLogical( Solver % Values, 'Eliminate Linear Constraints', Found)
  IF(ASSOCIATED(RestMatrix)) THEN
    IF(.NOT.EliminateConstraints) &
      NumberOfRows = NumberOFRows + RestMatrix % NumberOfRows
  END IF

  ALLOCATE( CollectionMatrix % RHS( NumberOfRows ), &
       CollectionSolution( NumberOfRows ), STAT = istat )
  IF ( istat /= 0 ) CALL Fatal( 'SolveWithLinearRestriction', 'Memory allocation error.' )

  CollectionVector => CollectionMatrix % RHS
  CollectionVector = 0.0_dp
  CollectionSolution = 0.0_dp

!------------------------------------------------------------------------------
! If multiplier should be exported,  allocate memory and export the variable.
!------------------------------------------------------------------------------

  ExportMultiplier = ListGetLogical( Solver % Values, 'Export Lagrange Multiplier', Found )
  IF ( .NOT. Found ) ExportMultiplier = .FALSE.


  IF ( ExportMultiplier ) THEN
     MultiplierName = ListGetString( Solver % Values, 'Lagrange Multiplier Name', Found )
     IF ( .NOT. Found ) THEN
        CALL Info( 'SolveWithLinearRestriction', &
              'Lagrange Multiplier Name set to LagrangeMultiplier', Level=5 )
        MultiplierName = "LagrangeMultiplier"
     END IF

     MultVar => VariableGet(Solver % Mesh % Variables, MultiplierName)
     j = 0
     IF(ASSOCIATED(RestMatrix)) j = RestMatrix % NumberofRows
     IF(ASSOCIATED(AddMatrix))  j = j+MAX(0,AddMatrix % NumberofRows-StiffMatrix % NumberOfRows)

     IF ( .NOT. ASSOCIATED(MultVar) ) THEN
       ALLOCATE( MultiplierValues(j), STAT=istat )
       IF ( istat /= 0 ) CALL Fatal('SolveWithLinearRestriction','Memory allocation error.')

       MultiplierValues = 0.0_dp
       CALL VariableAdd(Solver % Mesh % Variables, Solver % Mesh, SolverPointer, &
                  MultiplierName, 1, MultiplierValues)
     END IF
     MultVar => VariableGet(Solver % Mesh % Variables, MultiplierName)

     MultiplierValues => MultVar % Values

     IF (j>SIZE(MultiplierValues)) THEN
       ALLOCATE(MultiplierValues(j)); MultiplierValues=0._dp
       MultiplierValues(1:SIZE(MultVar % Values)) = MultVar % Values
       DEALLOCATE(MultVar % Values)
       MultVar % Values => MultiplierValues
     END IF
  ELSE
     MultiplierValues => NULL()
  END IF


!------------------------------------------------------------------------------
! Put the RestMatrix to lower part of CollectionMatrix
!------------------------------------------------------------------------------

  EnforceDirichlet = ListGetLogical( Solver % Values, 'Enforce Exact Dirichlet BCs',Found)
  IF(.NOT.Found) EnforceDirichlet = .TRUE.
  EnforceDirichlet = EnforceDirichlet .AND. ALLOCATED(StiffMatrix % ConstrainedDOF)

  ComplexSystem = StiffMatrix % COMPLEX
  ComplexSystem = ComplexSystem .OR. ListGetLogical( Solver % Values, &
           'Linear System Complex', Found )

  UseTranspose = ListGetLogical( Solver % Values, 'Use Transpose values', Found)

  IF(ASSOCIATED(RestMatrix).AND..NOT.EliminateConstraints) THEN

    CALL Info('SolveWithLinearRestriction',&
        'Adding ConstraintMatrix into CollectionMatrix',Level=10)


    NoEmptyRows = 0
    ConstraintScaling = ListGetLogical(Solver % Values, 'Constraint Scaling',Found)
    IF(ConstraintScaling) THEN
      rowsum = ListGetConstReal( Solver % Values, 'Constraint Scale', Found)
      IF(Found) RestMatrix % Values = RestMatrix % Values * rowsum
    END IF

    DO i=RestMatrix % NumberOfRows,1,-1

      k=StiffMatrix % NumberOfRows
      IF(ASSOCIATED(AddMatrix)) k=MAX(k,AddMatrix % NumberOfRows)
      k=k+i

      CALL AddToMatrixElement( CollectionMatrix,k,k,0._dp )
      IF(ComplexSystem) THEN
        IF(MOD(k,2)==0) THEN
          CALL AddToMatrixElement( CollectionMatrix,k,k-1,0._dp )
        ELSE
          CALL AddToMatrixElement( CollectionMatrix,k,k+1,0._dp )
        END IF
      END IF
      EmptyRow = .TRUE.

      rowsum = 0._dp
      l = -1
      DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
        IF(RestMatrix % Cols(j)==k) l=j
        rowsum = rowsum + ABS(RestMatrix % Values(j))
      END DO

      IF(rowsum>EPSILON(1._dp)) THEN
        IF(ConstraintScaling) THEN
          IF(l<=0.OR.l>0.AND.RestMatrix % Values(l)==0) THEN
            DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
              RestMatrix % Values(j) = RestMatrix % values(j)/rowsum
            END DO
            RestMatrix % RHS(i) = RestMatrix % RHS(i) / rowsum
          END IF
        END IF

        DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
          Found = .TRUE.

          ! Skip non-positive column indexes
          IF( RestMatrix % Cols(j) <= 0 ) CYCLE
          IF ( .NOT. ComplexSystem ) THEN
            IF( ABS(RestMatrix % Values(j)) < EPSILON(1._dp)*rowsum ) CYCLE
          END IF

          EmptyRow = .FALSE.

          IF (EnforceDirichlet .AND. RestMatrix % Cols(j) <= StiffMatrix % NumberOfRows) &
                  Found = .NOT.StiffMatrix % ConstrainedDOF(RestMatrix % Cols(j))

          IF(Found) THEN
            IF (ASSOCIATED(RestMatrix % TValues)) THEN
              CALL AddToMatrixElement( CollectionMatrix, &
                 RestMatrix % Cols(j), k, RestMatrix % TValues(j))
            ELSE
              CALL AddToMatrixElement( CollectionMatrix, &
                 RestMatrix % Cols(j), k, RestMatrix % Values(j))
            END IF
          END IF

          IF (UseTranspose .AND. ASSOCIATED(RestMatrix % TValues)) THEN
            CALL AddToMatrixElement( CollectionMatrix, &
                     k, RestMatrix % Cols(j), RestMatrix % TValues(j))
          ELSE
            CALL AddToMatrixElement( CollectionMatrix, &
                    k, RestMatrix % Cols(j), RestMatrix % Values(j))
          END IF
        END DO
      END IF

      IF (EnforceDirichlet) THEN
        IF(ASSOCIATED(RestMatrix % InvPerm)) THEN
          l = RestMatrix % InvPerm(i)
          IF(l>0) THEN
            l = MOD(l-1,StiffMatrix % NumberOfRows)+1
            IF(StiffMatrix % ConstrainedDOF(l)) THEN
              CollectionVector(k) = 0
              CALL ZeroRow(CollectionMatrix,k)
              CALL SetMatrixElement(CollectionMatrix,k,k,1._dp)
            END IF
          END IF
        END IF
      END IF
      
      ! If there is no matrix entry, there can be no non-zero r.h.s.
      IF( EmptyRow ) THEN
        NoEmptyRows = NoEmptyRows + 1
        CollectionVector(k) = 0._dp
!        might not be the right thing to do in parallel!!
!       CALL SetMatrixElement( CollectionMatrix,k,k,1._dp )
      ELSE
        IF( ASSOCIATED( RestVector ) ) CollectionVector(k) = RestVector(i)
      END IF
    END DO

    IF( NoEmptyRows > 0 ) THEN
      CALL Info('SolveWithLinearRestriction',&
          'Constraint Matrix in partition '//TRIM(I2S(ParEnv % MyPe))// &
          ' has '//TRIM(I2S(NoEmptyRows))// &
          ' empty rows out of '//TRIM(I2S(RestMatrix % NumberOfRows)), &
	  Level=6 )
    END IF

    CALL Info('SolveWithLinearRestriction',&
        'Finished Adding ConstraintMatrix',Level=12)
  END IF

!------------------------------------------------------------------------------
! Put the AddMatrix to upper part of CollectionMatrix
!------------------------------------------------------------------------------
  IF(ASSOCIATED(AddMatrix)) THEN

    CALL Info('SolveWithLinearRestriction',&
        'Adding AddMatrix into CollectionMatrix',Level=10)

    DO i=AddMatrix % NumberOfRows,1,-1

      Found = .TRUE.
      IF (EnforceDirichlet .AND. i<=StiffMatrix % NumberOFRows) &
         Found = .NOT.StiffMatrix % ConstrainedDOF(i)

      IF(Found) THEN
        Found = .FALSE.
        DO j=AddMatrix % Rows(i+1)-1,AddMatrix % Rows(i),-1
            CALL AddToMatrixElement( CollectionMatrix, &
               i, AddMatrix % Cols(j), AddMatrix % Values(j))
            IF (i == AddMatrix % Cols(j)) Found = .TRUE.
        END DO

        CollectionVector(i) = CollectionVector(i) + AddVector(i)
        IF (.NOT.Found) THEN
          CALL AddToMatrixElement( CollectionMatrix, i, i, 0._dp )
          IF(ComplexSystem) THEN
            IF(MOD(i,2)==0) THEN
              CALL AddToMatrixElement( CollectionMatrix,i,i-1,0._dp )
            ELSE
              CALL AddToMatrixElement( CollectionMatrix,i,i+1,0._dp )
            END IF
          END IF
        END IF
      END IF
    END DO
    CALL Info('SolveWithLinearRestriction',&
        'Finished Adding AddMatrix',Level=12)
  END IF

!------------------------------------------------------------------------------
! Put the StiffMatrix to upper part of CollectionMatrix
!------------------------------------------------------------------------------
  CALL Info('SolveWithLinearRestriction',&
      'Adding Stiffness Matrix into CollectionMatrix',Level=10)

  DO i=StiffMatrix % NumberOfRows,1,-1
    DO j=StiffMatrix % Rows(i+1)-1,StiffMatrix % Rows(i),-1
      CALL AddToMatrixElement( CollectionMatrix, &
        i, StiffMatrix % Cols(j), StiffMatrix % Values(j) )
    END DO
    CollectionVector(i) = CollectionVector(i) + ForceVector(i)
  END DO

!------------------------------------------------------------------------------
! Eliminate constraints instead of adding the Lagrange coefficient equations.
! Assumes biorthogonal basis for Lagrange coefficient interpolation, but not
! necesserily biorthogonal constraint equation test functions.
!------------------------------------------------------------------------------
  IF (ASSOCIATED(RestMatrix).AND.EliminateConstraints) THEN
    CALL Info('SolveWithLinearRestriction',&
        'Eliminating Constraints from CollectionMatrix',Level=10)

    n = StiffMatrix % NumberOfRows
    m = RestMatrix % NumberOfRows

    ALLOCATE(SlaveDiag(m),MasterDiag(m),SlavePerm(n),MasterPerm(n),SlaveIPerm(m),MasterIPerm(m),DiagDiag(m))
    SlavePerm  = 0; SlaveIPerm  = 0; 
    MasterPerm = 0; MasterIPerm = 0

    Tvals => RestMatrix % TValues
    IF (.NOT.ASSOCIATED(Tvals)) Tvals => RestMatrix % Values 

    ! Extract diagonal entries for constraints:
    !------------------------------------------
    CALL Info('SolveWithLInearRestriction',&
        'Extracting diagonal entries for constraints',Level=15)
    DO i=1, RestMatrix % NumberOfRows
      m = RestMatrix % InvPerm(i)

      IF( m == 0 ) THEN
        PRINT *,'InvPerm is zero:',ParEnv % MyPe, i
        CYCLE
      END IF

      m = MOD(m-1,n) + 1
      SlavePerm(m)  = i
      SlaveIperm(i) = m

      DO j=RestMatrix % Rows(i), RestMatrix % Rows(i+1)-1
        k = RestMatrix % Cols(j)
        IF(k>n) THEN
           DiagDiag(i) = Tvals(j)
           CYCLE
        END IF

        IF( ABS( TVals(j) ) < TINY( 1.0_dp ) ) THEN
          PRINT *,'Tvals too small',ParEnv % MyPe,j,i,k,RestMatrix % InvPerm(i),Tvals(j)
        END IF

        IF(k == RestMatrix % InvPerm(i)) THEN
           SlaveDiag(i) = Tvals(j)
        ELSE
           MasterDiag(i) = Tvals(j)
           MasterPerm(k)  = i
           MasterIperm(i) = k
        END IF
      END DO
    END DO
  END IF

  IF (ASSOCIATED(RestMatrix).AND.EliminateConstraints) THEN
    EliminateSlave = ListGetLogical( Solver % values, 'Eliminate Slave',Found )
    EliminateFromMaster = ListGetLogical( Solver % values, 'Eliminate From Master',Found )

    IF(EliminateFromMaster) THEN
      CALL Info('SolveWithLInearRestriction',&
          'Eliminating from master',Level=15)      
      UsePerm  => MasterPerm 
      UseDiag  => MasterDiag
      UseIPerm => MasterIPerm 
    ELSE
      CALL Info('SolveWithLInearRestriction',&
          'Eliminating from slave',Level=15)            
      UsePerm  => SlavePerm
      UseDiag  => SlaveDiag
      UseIPerm => SlaveIPerm
    END IF

    IF(UseTranspose) THEN
      Vals => Tvals
    ELSE
      Vals => RestMatrix % Values
    END IF
  END IF

  IF ( ParEnv % Pes>1 ) THEN
    EliminateDiscont =  ListGetLogical( Solver % values, 'Eliminate Discont',Found )
    IF( EliminateDiscont ) THEN
      CALL totv( StiffMatrix, SlaveDiag, SlaveIPerm )
      CALL totv( StiffMatrix, DiagDiag, SlaveIPerm )
      CALL totv( StiffMatrix, MasterDiag, MasterIPerm )
      CALL tota( StiffMatrix, TotValues, SlavePerm )
    END IF
  ELSE
    EliminateDiscont = .FALSE.
  END IF

  IF (ASSOCIATED(RestMatrix).AND.EliminateConstraints) THEN
    ! Replace elimination equations by the constraints (could done be as a postprocessing
    ! step, if eq's totally eliminated from linsys.)
    ! ----------------------------------------------------------------------------------
    CALL Info('SolveWithLInearRestriction',&
        'Deleting rows from equation to be eliminated',Level=15)

    Lmat => CollectionMatrix % ListMatrix
    DO m=1,RestMatrix % NumberOfRows
      i = UseIPerm(m)
      CALL List_DeleteRow(Lmat, i, Keep=.TRUE.)
    END DO

    CALL Info('SolveWithLInearRestriction',&
        'Copying rows from constraint matrix to eliminate dofs',Level=15)
    DO m=1,RestMatrix % NumberOfRows
      i = UseIPerm(m)
      DO l=RestMatrix % Rows(m+1)-1, RestMatrix % Rows(m), -1
        j = RestMatrix % Cols(l)

        ! skip l-coeffient entries, handled separately afterwards:
        ! --------------------------------------------------------
        IF(j > n) CYCLE

        CALL List_AddToMatrixElement( Lmat, i, j, Vals(l) )
      END DO
      CollectionVector(i) = RestVector(m)
    END DO

    ! Eliminate slave dof cycles:
    ! ---------------------------
    Xmat => RestMatrix
    Found = .TRUE.
    Loop = 0
    DO WHILE(Found)
      DO i=Xmat % NumberofRows,1,-1
        q = 0
        DO j = Xmat % Rows(i+1)-1, Xmat % Rows(i),-1
          k = Xmat % Cols(j)
          IF(k>n) CYCLE
          IF(UsePerm(k)>0 .AND. ABS(TVals(j))>AEPS) q=q+1
        END DO
        IF(q>1) EXIT
      END DO
      Found = q>1

      Tmat => Xmat
      IF(Found) THEN
        Loop = Loop + 1
        CALL Info('SolveWithLInearRestriction',&
            'Recursive elimination round: '//TRIM(I2S(Loop)),Level=15)

        Tmat => AllocateMatrix()
        Tmat % Format = MATRIX_LIST

        DO i=Xmat % NumberofRows,1,-1
          DO j = Xmat % Rows(i+1)-1, Xmat % Rows(i),-1
            k = Xmat % Cols(j)
            IF ( ABS(Tvals(j))>AEPS ) &
              CALL List_AddToMatrixElement(Tmat % ListMatrix, i, k, TVals(j))
          END DO
        END DO


        DO m=1,Xmat % NumberOfRows
          i = UseIPerm(m)
          DO j=Xmat % Rows(m), Xmat % Rows(m+1)-1
            k = Xmat % Cols(j)

            ! The size of SlavePerm is often exceeded but I don't really undersrtand the operation...
            ! so this is just a dirty fix.
            IF( k > SIZE( SlavePerm ) ) CYCLE

            l = SlavePerm(k)

            IF(l>0 .AND. k/=i) THEN
              IF(ABS(Tvals(j))<AEPS) CYCLE
              scl = -TVals(j) / SlaveDiag(l)

              CALL List_DeleteMatrixElement( Tmat % ListMatrix, m, k )

              DO q=Xmat % Rows(l+1)-1, Xmat % Rows(l),-1
                IF(ABS(Tvals(q))<AEPS) CYCLE
                ix = Xmat % Cols(q)
                IF ( ix/=k ) &
                  CALL List_AddToMatrixElement( Tmat % ListMatrix, m, ix, scl * TVals(q) )
              END DO
            END IF
          END DO
        END DO

        CALL List_ToCRSMatrix(Tmat)
        Tvals => Tmat % Values
        IF(.NOT.ASSOCIATED(Xmat,RestMatrix)) CALL FreeMatrix(Xmat)
      END IF
      Xmat => TMat
    END DO

    ! Eliminate Lagrange Coefficients:
    ! --------------------------------

    CALL Info('SolveWithLInearRestriction',&
        'Eliminating Largrange Coefficients',Level=15)

    DO m=1,Tmat % NumberOfRows
      i = UseIPerm(m)
      IF( ABS( UseDiag(m) ) < TINY( 1.0_dp ) ) THEN
        PRINT *,'UseDiag too small:',m,ParEnv % MyPe,UseDiag(m)
        CYCLE
      END IF

      DO j=TMat % Rows(m), TMat % Rows(m+1)-1
        k = TMat % Cols(j)
        IF(k<=n) THEN
          IF(UsePerm(k)/=0) CYCLE

          IF ( EliminateDiscont ) THEN
            IF (EliminateFromMaster) THEN
              scl = -SlaveDiag(SlavePerm(k)) / UseDiag(m)
            ELSE
              scl = -MasterDiag(MasterPerm(k)) / UseDiag(m)
            END IF
          ELSE
            scl = -Tvals(j) / UseDiag(m)
          END IF
        ELSE
          k = UseIPerm(k-n)
          ! multiplied by 1/2 in GenerateConstraintMatrix()
          IF (EliminateDiscont) THEN
            scl = -2*DiagDiag(m) / UseDiag(m)
          ELSE
            scl = -2*Tvals(j) / UseDiag(m)
          END IF
        END IF

        DO l=StiffMatrix % Rows(i+1)-1, StiffMatrix % Rows(i),-1
          CALL List_AddToMatrixElement( Lmat, k, &
              StiffMatrix % Cols(l), scl * StiffMatrix % Values(l) )
        END DO
        CollectionVector(k) = CollectionVector(k) + scl * ForceVector(i)
      END DO
    END DO

    IF ( .NOT.ASSOCIATED(Tmat, RestMatrix ) ) CALL FreeMatrix(Tmat)

    ! Eliminate slave dofs, using the constraint equations:
    ! -----------------------------------------------------
    IF ( EliminateSlave ) THEN
      CALL Info('SolveWithLInearRestriction',&
          'Eliminate slave dofs using constraint equations',Level=15)

      IF(EliminateDiscont) THEN
        DO i=1,StiffMatrix % NumberOfRows
          IF ( UsePerm(i)/=0 ) CYCLE

          DO m=StiffMatrix % Rows(i), StiffMatrix % Rows(i+1)-1
             j = SlavePerm(StiffMatrix % Cols(m))
             IF ( j==0 ) CYCLE
             scl = -TotValues(m) / SlaveDiag(j)

             ! Delete elimination entry:
             ! -------------------------
             CALL List_DeleteMatrixElement(Lmat,i,StiffMatrix % Cols(m))

             k = UseIPerm(j)
             cTmp => Lmat(k) % Head
             DO WHILE(ASSOCIATED(cTmp))
                l = cTmp % Index
                IF ( l /= SlaveIPerm(j) ) &
                   CALL List_AddToMatrixElement( Lmat, i, l, scl*cTmp % Value )
              cTmp => cTmp % Next
            END DO
            CollectionVector(i) = CollectionVector(i) + scl * CollectionVector(k)
          END DO
        END DO
      ELSE

        CALL List_ToCRSMatrix(CollectionMatrix)
        Tmat => AllocateMatrix()
        Tmat % Format = MATRIX_LIST

        DO i=1,StiffMatrix % NumberOfRows
          IF(UsePerm(i)/=0) CYCLE

          DO m = CollectionMatrix % Rows(i), CollectionMatrix % Rows(i+1)-1
            j = SlavePerm(CollectionMatrix % Cols(m))

            IF(j==0) THEN
              CYCLE
            END IF
            IF( ABS( SlaveDiag(j) ) < TINY( 1.0_dp ) ) THEN
              PRINT *,'SlaveDiag too small:',j,ParEnv % MyPe,SlaveDiag(j)
              CYCLE
            END IF

            scl = -CollectionMatrix % Values(m) / SlaveDiag(j)
            CollectionMatrix % Values(m) = 0._dp

            ! ... and add replacement values:
            ! -------------------------------
            k = UseIPerm(j)
            DO p=CollectionMatrix % Rows(k+1)-1, CollectionMatrix % Rows(k), -1
               l = CollectionMatrix % Cols(p)
               IF ( l /= SlaveIPerm(j) ) &
                 CALL List_AddToMatrixElement( Tmat % listmatrix, i, l, scl*CollectionMatrix % Values(p) )
            END DO
            CollectionVector(i) = CollectionVector(i) + scl * CollectionVector(k)
          END DO
        END DO

        CALL List_ToListMatrix(CollectionMatrix)
        Lmat => CollectionMatrix % ListMatrix

        CALL List_ToCRSMatrix(Tmat)
        DO i=TMat % NumberOfRows,1,-1
          DO j=TMat % Rows(i+1)-1,TMat % Rows(i),-1
            CALL List_AddToMatrixElement( Lmat, i, TMat % cols(j), TMat % Values(j) )
          END DO
        END DO
        CALL FreeMatrix(Tmat)
      END IF
    END IF

    ! Optimize bandwidth, if needed:
    ! ------------------------------
    IF(EliminateFromMaster) THEN
      CALL Info('SolveWithLinearRestriction',&
          'Optimizing bandwidth after elimination',Level=15)
      DO i=1,RestMatrix % NumberOfRows
        j = SlaveIPerm(i)
        k = MasterIPerm(i)

        Ctmp => Lmat(j) % Head
        Lmat(j) % Head => Lmat(k) % Head
        Lmat(k) % Head => Ctmp

        l = Lmat(j) % Degree
        Lmat(j) % Degree = Lmat(k) % Degree
        Lmat(k) % Degree = l

        scl = CollectionVector(j)
        CollectionVector(j) = CollectionVector(k)
        CollectionVector(k) = scl
      END DO
    END IF

    CALL Info('SolveWithLinearRestriction',&
        'Finished Adding ConstraintMatrix',Level=12)
  END IF

  CALL Info('SolveWithLinearRestriction',&
      'Reverting CollectionMatrix back to CRS matrix',Level=10)
  IF(CollectionMatrix % FORMAT==MATRIX_LIST) &
    CALL List_toCRSMatrix(CollectionMatrix)

  CALL Info( 'SolveWithLinearRestriction', 'CollectionMatrix done', Level=5 )

!------------------------------------------------------------------------------
! Assign values to CollectionVector
!------------------------------------------------------------------------------

  j = StiffMatrix % NumberOfRows  
  CollectionSolution(1:j) = Solution(1:j)
  
  i = StiffMatrix % NumberOfRows+1
  j = SIZE(CollectionSolution)
  CollectionSolution(i:j) = 0._dp
  IF(ExportMultiplier) CollectionSolution(i:j) = MultiplierValues(1:j-i+1)

  CollectionMatrix % ExtraDOFs = CollectionMatrix % NumberOfRows - &
                  StiffMatrix % NumberOfRows

  CollectionMatrix % ParallelDOFs = 0
  IF(ASSOCIATED(AddMatrix)) &
    CollectionMatrix % ParallelDOFs = MAX(AddMatrix % NumberOfRows - &
                  StiffMatrix % NumberOfRows,0)

  CALL Info( 'SolveWithLinearRestriction', 'CollectionVector done', Level=5 )

!------------------------------------------------------------------------------
! Solve the Collection-system 
!------------------------------------------------------------------------------

! Collectionmatrix % Complex = StiffMatrix % Complex

  ! We may want to skip the constraints for norm if we use certain other options
  SkipConstraints = ListGetLogical( Solver % values, &
      'Nonlinear System Convergence Without Constraints',Found )
  IF(.NOT. Found ) THEN
    SkipConstraints = ListGetLogical( Solver % values, 'Linear System Residual Mode',Found ) 
    IF( SkipConstraints ) THEN
      CALL Info('SolveWithLinearRestriction','Linear system residual mode must skip constraints',Level=5)
    ELSE
      SkipConstraints = ListGetLogical( Solver % values, 'NonLinear System Consistent Norm',Found ) 
      IF( SkipConstraints ) THEN
        CALL Info('SolveWithLinearRestriction','Nonlinear system consistent norm must skip constraints',Level=5)
      END IF
    END IF
    IF( SkipConstraints ) THEN
      CALL Info('SolveWithLinearRestriction','Enforcing convergence without constraints to True',Level=5)
      CALL ListAddLogical( Solver % Values, &
           'Nonlinear System Convergence Without Constraints',.TRUE.)
    END IF
  END IF

  !------------------------------------------------------------------------------
  ! Look at the nonlinear system previous values again, not taking the constrained
  ! system into account...
  !------------------------------------------------------------------------------
  Found = ASSOCIATED(Solver % Variable % NonlinValues)
  IF( Found .AND. .NOT. SkipConstraints ) THEN
    k = CollectionMatrix % NumberOfRows
    IF ( SIZE(Solver % Variable % NonlinValues) /= k) THEN
      DEALLOCATE(Solver % Variable % NonlinValues)
      ALLOCATE(Solver % Variable % NonlinValues(k))
    END IF
    Solver % Variable % NonlinValues(1:k) = CollectionSolution(1:k)
  END IF

  CollectionMatrix % Comm = StiffMatrix % Comm

  CALL Info('SolveWithLinearRestriction',&
      'Now going for the coupled linear system',Level=10)

  CALL SolveLinearSystem( CollectionMatrix, CollectionVector, &
      CollectionSolution, Norm, DOFs, Solver, StiffMatrix )


!------------------------------------------------------------------------------
! Separate the solution from CollectionSolution
!------------------------------------------------------------------------------
    CALL Info('SolveWithLinearRestriction',&
      'Picking solution from collection solution',Level=10)

    Solution = 0.0_dp
    i = 1
    j = StiffMatrix % NumberOfRows
    Solution(i:j) = CollectionSolution(i:j)

    IF ( ExportMultiplier ) THEN
      i = StiffMatrix % NumberOfRows
      j=0
      IF(ASSOCIATED(RestMatrix)) j = RestMatrix % NumberOfRows
      IF(ASSOCIATED(AddMatrix)) &
        j=j+MAX(0,AddMatrix % NumberOfRows - StiffMatrix % NumberOFRows)

      MultiplierValues = 0.0_dp
      IF(ASSOCIATED(RestMatrix).AND.EliminateConstraints) THEN
        ! Compute eliminated l-coefficient values:
        ! ---------------------------------------
        DO i=1,RestMatrix % NumberOfRows
          scl = 1._dp / UseDiag(i)
          m = UseIPerm(i)
          MultiplierValues(i) = scl * ForceVector(m)
          DO j=StiffMatrix % Rows(m), StiffMatrix % Rows(m+1)-1
            MultiplierValues(i) = MultiplierValues(i) - &
              scl * StiffMatrix % Values(j) * Solution(StiffMatrix % Cols(j))
          END DO
        END DO
      ELSE
        MultiplierValues(1:j) = CollectionSolution(i+1:i+j)
      END IF

      IF(EliminateConstraints.AND.EliminateDiscont) THEN
        IF (EliminateFromMaster) THEN
          CALL totv(StiffMatrix,MultiplierValues,MasterIPerm)
        ELSE
          CALL totv(StiffMatrix,MultiplierValues,SlaveIPerm)
        END IF
      END IF
    END IF

!------------------------------------------------------------------------------

    StiffMatrix % CollectionMatrix => CollectionMatrix
    DEALLOCATE(CollectionSolution)
    CollectionMatrix % ConstraintMatrix => NULL()

    CALL Info( 'SolveWithLinearRestriction', 'All done', Level=5 )

CONTAINS

  SUBROUTINE totv( A, totvalues, perm )
    type(matrix_t), pointer :: A
    real(kind=dp) :: totvalues(:)
    integer, allocatable :: perm(:)

    real(kind=dp), ALLOCATABLE :: x(:),r(:)
    INTEGER :: i,j,ng

    ng = A % NumberOfRows
!   ng = ParallelReduction(1._dp*MAXVAL(A % ParallelInfo % GLobalDOfs))
    ALLOCATE(x(ng),r(ng))

    x = 0._dp
    IF(ALLOCATED(perm)) THEN
      DO i=1,SIZE(perm)
        j = Perm(i)
        !j = a % parallelinfo % globaldofs(j)
        x(j) = totvalues(i)
      END DO
    END IF

    CALL ParallelSumVector(A, x)
!   CALL MPI_ALLREDUCE( x,r, ng, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, i ); x=r

    IF(ALLOCATED(perm)) THEN
      DO i=1,SIZE(perm)
        j = Perm(i)
        !j = A % parallelinfo % globaldofs(j)
        totvalues(i) = x(j)
      END DO
    END IF
  END SUBROUTINE Totv
    

  SUBROUTINE Tota( A, TotValues, cperm )
     type(matrix_t), pointer :: A
     integer, allocatable :: cperm(:)
     real(kind=dp), ALLOCATABLE :: totvalues(:)

     INTEGER, POINTER :: Diag(:), Rows(:), Cols(:)
     LOGICAL ::  found
     INTEGER :: status(MPI_STATUS_SIZE)
     REAL(KIND=dp), ALLOCATABLE, TARGET :: rval(:)
     INTEGER, ALLOCATABLE :: cnt(:), rrow(:),rcol(:), perm(:)
     INTEGER :: i,j,k,l,m,ii,jj,proc,rcnt,nn, dof, dofs, Active, n, nm,ierr

     TYPE Buf_t
        REAL(KIND=dp), ALLOCATABLE :: gval(:)
        INTEGER, ALLOCATABLE :: grow(:),gcol(:)
     END TYPE Buf_t
     TYPE(Buf_t), POINTER :: buf(:)

     Diag => A % Diag
     Rows => A % Rows
     Cols => A % Cols

     n = A % NumberOfRows

     ALLOCATE(TotValues(SIZE(A % Values))); TotValues=A % Values

     IF (ParEnv  % PEs>1 ) THEN
       ALLOCATE(cnt(0:ParEnv % PEs-1))
       cnt = 0
       DO i=1,n
         DO j=Rows(i),Rows(i+1)-1
!          IF(Cols(j)<=nm .OR. Cols(j)>nm+n) CYCLE
           iF ( ALLOCATED(CPerm)) THEN
             IF(cperm(Cols(j))==0) CYCLE
           END IF
           IF(TotValues(j)==0) CYCLE

           IF ( A % ParallelInfo % Interface(Cols(j)) ) THEN
             DO k=1,SIZE(A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours)
               m = A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours(k)
               IF ( m==ParEnv % myPE ) CYCLE
               cnt(m) = cnt(m)+1
             END DO 
           END IF
         END DO
       END DO

       ALLOCATE( buf(0:ParEnv % PEs-1) )
       DO i=0,ParEnv % PEs-1
         IF ( cnt(i) > 0 ) &
           ALLOCATE( Buf(i) % gval(cnt(i)), Buf(i) % grow(cnt(i)), Buf(i) % gcol(cnt(i)) )
       END DO

       cnt = 0
       DO i=1,n
         DO j=Rows(i),Rows(i+1)-1
!          IF(Cols(j)<=nm .OR. Cols(j)>nm+n) CYCLE
           iF ( ALLOCATED(CPerm)) THEN
             IF(cperm(Cols(j))==0) CYCLE
           END IF
           IF(TotValues(j)==0) CYCLE

           IF ( A % ParallelInfo % Interface(Cols(j)) ) THEN
             DO k=1,SIZE(A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours)
               m = A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours(k)
               IF ( m==ParEnv % myPE ) CYCLE
               cnt(m) = cnt(m)+1
               Buf(m) % gcol(cnt(m)) = A % ParallelInfo % GlobalDOFs(Cols(j))
               Buf(m) % gval(cnt(m)) = TotValues(j)
               Buf(m) % grow(cnt(m)) = A % ParallelInfo % GlobalDOFs(i)
             END DO
           END IF
         END DO
       END DO

       DO i=0,ParEnv % PEs-1
         IF ( ParEnv % IsNeighbour(i+1) ) THEN
           CALL MPI_BSEND( cnt(i), 1, MPI_INTEGER, i, 7001, MPI_COMM_WORLD, status, ierr )
           IF ( cnt(i)>0 ) THEN
             CALL MPI_BSEND( Buf(i) % grow, cnt(i), MPI_INTEGER, &
                 i, 7002, MPI_COMM_WORLD, status, ierr )

             CALL MPI_BSEND( Buf(i) % gcol, cnt(i), MPI_INTEGER, &
                 i, 7003, MPI_COMM_WORLD, status, ierr )

             CALL MPI_BSEND( Buf(i) % gval, cnt(i), MPI_DOUBLE_PRECISION, &
                 i, 7004, MPI_COMM_WORLD, status, ierr )
           END IF
         END IF
       END DO

       DO i=0,ParEnv % PEs-1
         IF ( cnt(i)>0 ) &
           DEALLOCATE( Buf(i) % gval, Buf(i) % grow, Buf(i) % gcol )
       END DO
       DEALLOCATE( cnt,Buf )

       DO i=1,ParEnv % NumOfNeighbours
         CALL MPI_RECV( rcnt, 1, MPI_INTEGER, &
           MPI_ANY_SOURCE, 7001, MPI_COMM_WORLD, status, ierr )

         IF ( rcnt>0 ) THEN
           IF(.NOT.ALLOCATED(rrow)) THEN
             ALLOCATE( rrow(rcnt), rcol(rcnt), rval(rcnt) )
           ELSE IF(SIZE(rrow)<rcnt) THEN
             DEALLOCATE(rrow,rcol,rval)
             ALLOCATE( rrow(rcnt), rcol(rcnt), rval(rcnt) )
           ENDIF

           proc = status(MPI_SOURCE)
           CALL MPI_RECV( rrow, rcnt, MPI_INTEGER, &
              proc, 7002, MPI_COMM_WORLD, status, ierr )

           CALL MPI_RECV( rcol, rcnt, MPI_INTEGER, &
              proc, 7003, MPI_COMM_WORLD, status, ierr )

           CALL MPI_RECV( rval, rcnt, MPI_DOUBLE_PRECISION, &
              proc, 7004, MPI_COMM_WORLD, status, ierr )

           DO j=1,rcnt
             l = SearchNode(A % ParallelInfo,rcol(j),Order=A % Perm)
             IF ( l>0 ) THEN
               k = SearchNode(A % ParallelInfo,rrow(j),Order=A % Perm)
               IF ( k>0 ) THEN
                 IF ( l>=k ) THEN
                   DO m=Diag(k),Rows(k+1)-1
                     IF ( Cols(m) == l ) THEN
                       TotValues(m) = TotValues(m) + rval(j)
                       EXIT
                     ELSE IF( Cols(m)>l) THEN
                       EXIT
                     END IF
                   END DO
                 ELSE
                   DO m=Rows(k),Diag(k)-1
                     IF ( Cols(m)==l ) THEN
                       TotValues(m) = TotValues(m) + rval(j)
                       EXIT
                     ELSE IF( Cols(m)>l) THEN
                       EXIT
                     END IF
                   END DO
                 END IF
               END IF
             END IF
           END DO
         END IF
       END DO
     END IF
  End subroutine tota

!------------------------------------------------------------------------------
  END SUBROUTINE SolveWithLinearRestriction
!------------------------------------------------------------------------------
      

!------------------------------------------------------------------------------
  SUBROUTINE SaveLinearSystem( Solver, Ain )
!------------------------------------------------------------------------------
    TYPE( Solver_t ) :: Solver
    TYPE(Matrix_t), POINTER, OPTIONAL :: Ain
!------------------------------------------------------------------------------    
    TYPE(Matrix_t), POINTER :: A
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(LEN=MAX_NAME_LEN) :: dumpfile, dumpprefix
    INTEGER, POINTER :: Perm(:)
    INTEGER :: i
    LOGICAL :: SaveMass, SaveDamp, SavePerm, Found , Parallel, CNumbering
!------------------------------------------------------------------------------

    CALL Info('SaveLinearSystem','Saving linear system',Level=4)

    Parallel = ParEnv % PEs > 1

    Params => Solver % Values
    IF(.NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal('SaveLinearSystem','Parameter list not associated!')
    END IF

    CNumbering = ListGetLogical(Params, 'Linear System Save Continuous Numbering',Found)

    IF( PRESENT(Ain)) THEN
      A => Ain
    ELSE
      A => Solver % Matrix
    END IF

    IF(.NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal('SaveLinearSystem','Matrix not assciated!')
    END IF

    SaveMass = ListGetLogical( Params,'Linear System Save Mass',Found)
    SaveDamp = ListGetLogical( Params,'Linear System Save Damp',Found)
    
    dumpprefix = ListGetString( Params, 'Linear System Save Prefix', Found)
    IF(.NOT. Found ) dumpprefix = 'linsys'

    dumpfile = TRIM(dumpprefix)//'_a.dat'
    IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))
    CALL Info('SaveLinearSystem','Saving matrix to: '//TRIM(dumpfile))
    OPEN(1,FILE=dumpfile, STATUS='Unknown')
    CALL PrintMatrix(A,Parallel,Cnumbering,SaveMass=SaveMass,SaveDamp=SaveDamp)
    CLOSE(1)

    dumpfile = TRIM(dumpprefix)//'_b.dat'
    IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))
    CALL Info('SaveLinearSystem','Saving matrix rhs to: '//TRIM(dumpfile))
    OPEN(1,FILE=dumpfile, STATUS='Unknown')
    CALL PrintRHS(A, Parallel, CNumbering)
    CLOSE(1)
    
    SavePerm = ListGetLogical( Params,'Linear System Save Perm',Found)
    IF( SavePerm ) THEN
      Perm => Solver % Variable % Perm
      IF( .NOT. ASSOCIATED( Perm ) ) THEN
        CALL Warn('SaveLinearSystem','Permuation not associated!')
        SavePerm = .FALSE.
      ELSE
        dumpfile = TRIM(dumpprefix)//'_perm.dat'
        IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))
        CALL Info('SaveLinearSystem','Saving permutation to: '//TRIM(dumpfile))
        OPEN(1,FILE=dumpfile, STATUS='Unknown')
        DO i=1,SIZE(Perm)
          WRITE(1,'(I0,A,I0)') i,' ',Perm(i)
        END DO
        CLOSE( 1 ) 
      END IF
    END IF
    
    dumpfile = TRIM(dumpprefix)//'_sizes.dat'
    IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))
    CALL Info('SaveLinearSystem','Saving matrix sizes to: '//TRIM(dumpfile))
    OPEN(1,FILE=dumpfile, STATUS='Unknown')
    WRITE(1,*) A % NumberOfRows
    WRITE(1,*) SIZE(A % Values)
    IF( SavePerm ) WRITE(1,*) SIZE( Perm )
    CLOSE(1)

    IF(Parallel) THEN
      dumpfile = TRIM(dumpprefix)//'_sizes.dat'
      CALL Info('SaveLinearSystem','Saving matrix sizes to: '//TRIM(dumpfile))
      OPEN(1,FILE=dumpfile, STATUS='Unknown')
      WRITE(1,*) NINT(ParallelReduction(1._dP*A % ParMatrix % &
                           SplittedMatrix % InsideMatrix % NumberOfRows))
      WRITE(1,*) NINT(ParallelReduction(1._dp*SIZE(A % Values)))
      IF( SavePerm ) WRITE(1,*) NINT(ParallelReduction(1._dp*SIZE( Perm )))
      CLOSE(1)
    END IF
    
    IF( ListGetLogical( Params,'Linear System Save and Stop',Found ) ) THEN
      CALL Info('SaveLinearSystem','Just saved matrix and stopped!',Level=4)
      STOP
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SaveLinearSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Assemble Laplace matrix related to a solver and permutation vector. 
!------------------------------------------------------------------------------
  SUBROUTINE LaplaceMatrixAssembly( Solver, Perm, A )
    
    TYPE(Solver_t) :: Solver
    INTEGER, POINTER :: Perm(:)
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh
    !------------------------------------------------------------------------------

    INTEGER, POINTER :: BoundaryPerm(:), Indexes(:)
    INTEGER :: i,j,k,n,t,istat,BoundaryNodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    CHARACTER(LEN=MAX_NAME_LEN) :: BoundaryName
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp) :: detJ, val
    LOGICAL :: Stat
    
    
    Mesh => Solver % Mesh
        
    N = Mesh % MaxElementNodes 
    ALLOCATE( Basis(n), dBasisdx(n, 3), FORCE(N), STIFF(N,N), &
        Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        STAT=istat)
    
    IF(.FALSE.) THEN
      N = Mesh % NumberOfNodes
      ALLOCATE( BoundaryPerm(n) )
      BoundaryPerm = 0
      BoundaryNodes = 0
      BoundaryName = 'Laplace Boundary'
      CALL MakePermUsingMask( CurrentModel,Solver,Mesh,BoundaryName, &
          .FALSE., BoundaryPerm, BoundaryNodes )
    END IF


    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      n = Element % TYPE % NumberOfNodes
      Indexes => Element % NodeIndexes
      IF( ANY( Perm(Indexes) == 0 ) ) CYCLE

      Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
      Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
      Nodes % z(1:n) = Mesh % Nodes % z(Indexes)

      STIFF = 0.0d0
      FORCE = 0.0d0
      
      ! Numerical integration:
      !----------------------
      IP = GaussPoints( Element )
      DO k=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
            IP % W(k),  detJ, Basis, dBasisdx )
        
        ! Finally, the elemental matrix & vector:
        !----------------------------------------
        DO i=1,n
          val = IP % s(k) * DetJ 
          
          ! This condition removes the natural boundary condition that would 
          ! try to fix the normal gradient of the field to zero.
          !--------------------------------------------------------------------
          IF(.FALSE.) THEN
            IF( BoundaryPerm( Indexes(i) ) > 0 ) CYCLE
          END IF

          DO j=1,n
            STIFF(i,j) = STIFF(i,j) + val * &
                SUM( dBasisdx(i,:) * dBasisdx(j,:) ) 
          END DO
        END DO
      END DO
      
      CALL UpdateGlobalEquations( A,STIFF,A % rhs,FORCE,n,1,Perm(Indexes(1:n)) )
    END DO
    
    DEALLOCATE( Basis, dBasisdx, FORCE, STIFF, & 
        Nodes % x, Nodes % y, Nodes % z)

  END SUBROUTINE LaplaceMatrixAssembly

 
!------------------------------------------------------------------------------
!> Assemble mass matrix related to a solver and permutation vector. 
!------------------------------------------------------------------------------
  SUBROUTINE MassMatrixAssembly( Solver, Perm, A )
    
    TYPE(Solver_t) :: Solver
    INTEGER, POINTER :: Perm(:)
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh
    !------------------------------------------------------------------------------

    INTEGER, POINTER :: Indexes(:)
    INTEGER :: i,j,k,n,t,istat
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp), POINTER :: Basis(:),rhs(:)
    REAL(KIND=dp) :: detJ, val
    LOGICAL :: Stat
    
    
    Mesh => Solver % Mesh
        
    N = Mesh % MaxElementNodes 
    ALLOCATE( Basis(n), FORCE(N), STIFF(N,N), &
        Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        STAT=istat)

    ALLOCATE( rhs(A % NumberOfRows) )
    rhs = 0.0_dp

    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      n = Element % TYPE % NumberOfNodes
      Indexes => Element % NodeIndexes
      IF( ANY( Perm(Indexes) == 0 ) ) CYCLE

      Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
      Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
      Nodes % z(1:n) = Mesh % Nodes % z(Indexes)

      STIFF = 0.0d0
      FORCE = 0.0d0
      
      ! Numerical integration:
      !----------------------
      IP = GaussPoints( Element )

      DO k=1,IP % n

        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
            IP % W(k),  detJ, Basis )
        
        ! Finally, the elemental matrix & vector:
        !----------------------------------------
        DO i=1,n
          val = IP % s(k) * DetJ 
          DO j=1,n
            STIFF(i,j) = STIFF(i,j) + val * Basis(i) * Basis(j)
          END DO
        END DO
      END DO

      CALL UpdateGlobalEquations( A,STIFF,rhs,FORCE,n,1,Perm(Indexes(1:n)) )
    END DO

    DEALLOCATE( Basis, FORCE, STIFF, & 
        Nodes % x, Nodes % y, Nodes % z)
    DEALLOCATE( rhs )

  END SUBROUTINE MassMatrixAssembly


!---------------------------------------------------------------------------------
!> Multiply a linear system by a constant or a given scalar field.
!
!> There are three multiplication modes:
!> 1) Multiply matrix or rhs with a constant factor
!> 2) Multiply matrix or rhs with a constant factor but only blockwise
!> 3) Multiply matrix or rhs with a vector retrieved by a field variable
!
!> And also three things to multiply:
!> a) The right-hand-side of the linear system
!> b) The matrix part of the linear system
!> c) The diagonal entries of the matrix
!
!> Possible uses of the routine include cases where the user wants to introduce diagonal
!> implicit relaxation to the linear system, or to eliminate some coupling terms in 
!> monolithic systems that make the solution of the linear problems more difficult.
!----------------------------------------------------------------------------------
  SUBROUTINE LinearSystemMultiply( Solver )
!----------------------------------------------------------------------------------    
    TYPE(Solver_t) :: Solver
    !------------------------------------------------------------------------------
    INTEGER, POINTER :: Perm(:),Rows(:),Cols(:)
    REAL(KIND=dp), POINTER :: Values(:),Rhs(:)
    TYPE(Variable_t), POINTER :: ThisVar,CoeffVar
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Coeff,Coeff2
    INTEGER :: i,j,j2,k,l,jk,n,Mode,Dofs
    LOGICAL :: Found, UpdateRhs, Symmetric
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(LEN=MAX_NAME_LEN) :: str, VarName

    Params => Solver % Values
    Mesh => Solver % Mesh
    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('LinearSystemMultiply','Subroutine requires a Mesh!')
    END IF
    A => Solver % Matrix
    IF(.NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal('LinearSystemMultiply','Subroutine requires a matrix equation!')
    END IF
    ThisVar => Solver % Variable
    IF(.NOT. ASSOCIATED( ThisVar ) ) THEN
      CALL Fatal('LinearSystemMultiply','Subroutine requires a default variable to exist!')
    END IF

    Perm => ThisVar % Perm
    Dofs = ThisVar % Dofs
    n = A % NumberOfRows
    Cols => A % Cols
    Rows => A % Rows
    Rhs => A % Rhs
    Values => A % Values
        
    UpdateRhs = ListGetLogical( Params,'Linear System Multiply Consistent',Found)
    Symmetric = ListGetLogical( Params,'Linear System Multiply Symmetric',Found)

    ! First, Multiply the k:th piece of the r.h.s. vector if requested
    !-----------------------------------------------------------
    DO k=1,Dofs
      Mode = 0
      
      WRITE( str,'(A)') 'Linear System Rhs Factor'
      Coeff = ListGetCReal( Params, str, Found )
      IF( Found ) THEN
        Mode = 1
        WRITE( Message,'(A,ES12.3)') 'Multiplying the rhs with ',Coeff
        CALL Info('LinearSystemMultiply',Message, Level=6 )
      ELSE
        WRITE( str,'(A,I0)') TRIM(str)//' ',k
        Coeff = ListGetCReal( Params, str, Found )
        IF( Found ) THEN
          Mode = 2 
          WRITE( Message,'(A,I0,A,ES12.3)') 'Multiplying component ',k,' of the rhs with ',Coeff
          CALL Info('LinearSystemMultiply',Message, Level=6 )
        END IF
      END IF
      IF( Mode == 0 ) THEN
        str = 'Linear System Rhs Variable'
        VarName = ListGetString( Params,str,Found ) 
        NULLIFY( CoeffVar ) 
        IF( Found ) THEN
          CoeffVar => VariableGet( Mesh % Variables, VarName )
        ELSE
          WRITE( str,'(A,I0)') TRIM(str)//' ',k
          VarName = ListGetString( Params,str,Found )
          IF( Found ) CoeffVar => VariableGet( Mesh % Variables, VarName )
        END IF
        IF( ASSOCIATED( CoeffVar ) ) THEN
          IF( ANY( CoeffVar % Perm /= Perm ) ) THEN
            CALL Fatal('LinearSystemMultiply','Permutations should be the same')
          END IF
          Mode = 3
          WRITE( Message,'(A,I0,A)') 'Multiplying component ',k,' of the rhs with > '//TRIM(VarName)//' <'
          CALL Info('LinearSystemMultiply',Message, Level=6 )

          !PRINT *,'Range:',Mode,MINVAL(CoeffVar % Values),MAXVAL(CoeffVar % Values)
        END IF
      END IF
      IF( Mode == 0 ) CYCLE
 
      IF( Mode == 1 ) THEN
        IF( ABS( Coeff - 1.0_dp ) > EPSILON( Coeff ) ) THEN
          Rhs = Coeff * Rhs
        END IF
        EXIT
      ELSE 
        DO j=1,SIZE( Perm ) 
          jk = Dofs*(j-1)+k
          IF( Mode == 3 ) Coeff = CoeffVar % Values(j)
          Rhs( jk ) = Coeff * Rhs( jk )
        END DO
      END IF
    END DO
    ! End of r.h.s. multiplication

    ! Secondly, Multiply the kl block of the matrix
    !------------------------------------------------
    DO k=1,Dofs
      DO l=1,Dofs
        Mode = 0
        str = 'Linear System Matrix Factor'
        Coeff = ListGetCReal( Params, str, Found )
        IF( Found ) THEN
          Mode = 1
          WRITE( Message,'(A,ES12.3)') 'Multiplying the matrix with ',Coeff
          CALL Info('LinearSystemMultiply',Message, Level=6 )
        ELSE
          WRITE( str,'(A,I0,I0)') TRIM(str)//' ',k,l
          Coeff = ListGetCReal( Params, str, Found )
          IF( Found ) THEN
            Mode = 2
            WRITE( Message,'(A,I0,I0,A,ES12.3)') 'Multiplying block (',k,l,') of the matrix with ',Coeff
            CALL Info('LinearSystemMultiply',Message, Level=6 )
          END IF
        END IF
        IF( Mode == 0 ) THEN
          str = 'Linear System Matrix Variable'
          VarName = ListGetString( Params,str,Found )
          NULLIFY( CoeffVar ) 
          IF( Found ) THEN
            CoeffVar => VariableGet( Mesh % Variables, str )                                    
          ELSE
            WRITE( str,'(A,I0,I0)') TRIM(str)//' ',k,l
            VarName = ListGetString( Params,str,Found )
            IF( Found ) CoeffVar => VariableGet( Mesh % Variables, VarName )
          END IF
          IF( ASSOCIATED( CoeffVar ) ) THEN
            IF( ANY( CoeffVar % Perm /= Perm ) ) THEN
              CALL Fatal('LinearSystemMultiply','Permutations should be the same')
            END IF
            Mode = 3
            WRITE( Message,'(A,I0,I0,A)') 'Multiplying block (',k,l,') of the matrix with > '//TRIM(VarName)//' <'
            CALL Info('LinearSystemMultiply',Message, Level=6 )
          END IF
        END IF

        IF( Mode == 0 ) CYCLE

        IF( Mode == 1 ) THEN
          IF( ABS( Coeff - 1.0_dp ) > EPSILON( Coeff ) ) THEN
            Values = Coeff * Values
          END IF
        ELSE 
          DO j=1,SIZE( Perm ) 
            jk = Dofs*(j-1)+k
            IF( Mode == 3 ) Coeff = CoeffVar % Values(j)
            DO i=Rows(jk),Rows(jk+1)-1 
              IF( MODULO( Cols(i), Dofs ) == MODULO( l, Dofs ) ) THEN
                IF( Mode == 3 .AND. Symmetric ) THEN          
                  j2 = (Cols(i)-1)/Dofs + 1 
                  Coeff2 = CoeffVar % Values(j2)
                  Values( i ) = SQRT( Coeff * Coeff2 ) * Values( i )
                ELSE
                  Values( i ) = Coeff * Values( i )
                END IF
              END IF
            END DO
          END DO
        END IF
      END DO
      IF( Mode == 1 ) EXIT
    END DO
    ! end of matrix multiplication


    ! Finally, Multiply the diagonal of the matrix
    !------------------------------------------------
    DO k=1,Dofs
      Mode = 0

      str = 'Linear System Diagonal Factor'
      Coeff = ListGetCReal( Params, str, Found )
      IF( Found ) THEN
        Mode = 1
        WRITE( Message,'(A,ES12.3)') 'Multiplying the diagonal with ',Coeff
        CALL Info('LinearSystemMultiply',Message, Level=6 )
      ELSE
        WRITE( str,'(A,I0)') TRIM(str)//' ',k
        Coeff = ListGetCReal( Params, str, Found )
        IF( Found ) THEN
          Mode = 2 
          WRITE( Message,'(A,I0,A,ES12.3)') 'Multiplying component ',k,' of the matrix diagonal with ',Coeff
          CALL Info('LinearSystemMultiply',Message, Level=6 )          
        END IF
      END IF

      IF( Mode == 0 ) THEN
        str = 'Linear System Diagonal Variable'
        VarName = ListGetString( Params,str,Found )
        NULLIFY( CoeffVar )
        IF( Found ) THEN
          CoeffVar => VariableGet( Mesh % Variables, VarName )
        ELSE
          WRITE( str,'(A,I0)') TRIM(str)//' ',k
          VarName = ListGetString( Params,str,Found )
          IF( Found ) CoeffVar => VariableGet( Mesh % Variables, VarName )
        END IF
        IF( ASSOCIATED( CoeffVar ) ) THEN
          IF( ANY( CoeffVar % Perm /= Perm ) ) THEN
            CALL Fatal('LinearSystemMultiply','Permutations should be the same')
          END IF
          Mode = 3
          WRITE( Message,'(A,I0,A)') 'Multiplying component ',k,' of the diagonal with > '//TRIM(VarName)//' <'
          CALL Info('LinearSystemMultiply',Message, Level=6 )
        END IF
      END IF
      
      IF( Mode == 0 ) CYCLE

      IF( Mode == 1 ) THEN
        IF( ABS( Coeff - 1.0_dp ) > EPSILON( Coeff ) ) THEN
          IF( UpdateRhs ) Rhs = Rhs + ( Coeff - 1 ) * Values( A % Diag ) * ThisVar % Values
          Values( A % Diag ) = Coeff * Values( A % Diag )
        END IF
        EXIT
      ELSE 
        DO j=1,SIZE( Perm ) 
          jk = Dofs*(j-1)+k
          IF( Mode == 3 ) Coeff = CoeffVar % Values(j)

          IF( UpdateRhs ) Rhs( jk ) = Rhs( jk ) + (Coeff - 1) * Values(A % Diag(jk)) * ThisVar % Values(jk)
          Values( A % Diag(jk) ) = Coeff * Values( A % Diag(jk) )
        END DO
      END IF
    END DO
    ! end of diagonal multiplication


  END SUBROUTINE LinearSystemMultiply






!---------------------------------------------------------------------------------
!> Set the diagonal entry to given minumum.
!----------------------------------------------------------------------------------
  SUBROUTINE LinearSystemMinDiagonal( Solver )
!----------------------------------------------------------------------------------    
    TYPE(Solver_t) :: Solver
    !------------------------------------------------------------------------------
    INTEGER, POINTER :: Perm(:),Rows(:),Cols(:)
    REAL(KIND=dp), POINTER :: Values(:),Rhs(:)
    TYPE(Variable_t), POINTER :: ThisVar
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Coeff
    INTEGER :: i,j,j2,k,l,jk,n,Mode,Dofs
    LOGICAL :: Found, UpdateRhs, Symmetric
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    INTEGER :: NoSet
    REAL(KIND=dp) :: DiagSum, val, DiagMax

    Params => Solver % Values
    Mesh => Solver % Mesh
    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('LinearSystemMinDiagonal','Subroutine requires a Mesh!')
    END IF
    A => Solver % Matrix
    IF(.NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal('LinearSystemMinDiagonal','Subroutine requires a matrix equation!')
    END IF
    ThisVar => Solver % Variable
    IF(.NOT. ASSOCIATED( ThisVar ) ) THEN
      CALL Fatal('LinearSystemMinDiagonal','Subroutine requires a default variable to exist!')
    END IF

    Perm => ThisVar % Perm
    Dofs = ThisVar % Dofs
    n = A % NumberOfRows
    Cols => A % Cols
    Rows => A % Rows
    Rhs => A % Rhs
    Values => A % Values
        
    ! Set the mimimum value for each component, only nodel dofs consired
    !-------------------------------------------------------------------
    NoSet = 0
    DiagMax = 0.0_dp
    DiagSum = 0.0_dp
    n = MAXVAL( Perm ( 1:Mesh % NumberOfNodes ) )

    DO k=1,Dofs
      Mode = 0

      str = 'Linear System Diagonal Min'
      Coeff = ListGetCReal( Params, str, Found )
      IF( Found ) THEN
        Mode = 1
        WRITE( Message,'(A,ES12.3)') 'Setting minimum of the diagonal to ',Coeff
        CALL Info('LinearSystemMinDiagonal',Message, Level=6 )
      ELSE
        WRITE( str,'(A,I0)') TRIM(str)//' ',k
        Coeff = ListGetCReal( Params, str, Found )
        IF( Found ) THEN
          Mode = 2 
          WRITE( Message,'(A,I0,A,ES12.3)') 'Setting minimum of diagonal component ',k,' to ',Coeff
          CALL Info('LinearSystemMinDiagonal',Message, Level=6 )          
        END IF
      END IF
      
      IF( Mode == 0 ) CYCLE
      
      DO j=1,n
        jk = Dofs*(j-1)+k
        l = A % Diag(jk) 
        IF( l == 0 ) CYCLE

        ! Only add the diagonal to the owned dof
        IF( ParEnv % PEs > 1 ) THEN
          IF( A % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
        END IF

        val = ABS( Values( l ) )
        DiagSum = DiagSum + val
        DiagMax = MAX( DiagMax, val )
        IF( val < Coeff ) THEN
          Values( A % Diag(jk) ) = Coeff
          NoSet = NoSet + 1
        END IF
      END DO
    END DO

    CALL Info('LinearSystemMinDiagonal','Number of diagonal values set to mimimum: '//TRIM(I2S(NoSet)),Level=5)
    WRITE( Message,'(A,ES12.3)') 'Avarage abs(diagonal) entry: ',DiagSum / n    
    CALL Info('LinearSystemMinDiagonal',Message, Level=6 )
    
    WRITE( Message,'(A,ES12.3)') 'Maximum abs(diagonal) entry: ',DiagMax
    CALL Info('LinearSystemMinDiagonal',Message, Level=6 )


  END SUBROUTINE LinearSystemMinDiagonal





  !----------------------------------------------------------------------
  !> Make the high-order flux corrected transport (FCT) correction after 
  !> the low order approximation has been solved. 
  !
  !> For more information see, for example, 
  !> Dmitri Kuzmin (2008): "Explicit and implicit FEM-FCT algorithms with flux linearization"
  !----------------------------------------------------------------------
  SUBROUTINE FCT_Correction( Solver  )

    TYPE(Solver_t), POINTER :: Solver

    TYPE(Valuelist_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: n,i,j,k,k2
    INTEGER, POINTER :: Rows(:),Cols(:),Perm(:)
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp), POINTER :: BulkValues(:),u(:),M_L(:),M_C(:),udot(:),SaveValues(:),&
        pp(:),pm(:),qp(:),qm(:),corr(:),ku(:),ulow(:)
    REAL(KIND=dp) :: rsum, Norm,m_ij,k_ij,du,d_ij,f_ij,c_ij,Ceps,CorrCoeff,&
        rmi,rmj,rpi,rpj,dt
    TYPE(Variable_t), POINTER :: Var, Variables
    LOGICAL :: Found, Symmetric, SaveFields, SkipCorrection
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, TmpVarName

    Params => Solver % Values

    SkipCorrection = ListGetLogical( Params,'FCT Correction Skip',Found )
    Symmetric = ListGetLogical( Params,'FCT Correction Symmetric',Found )
    SaveFields = ListGetLogical( Params,'FCT Correction Save',Found )

    IF( SkipCorrection .AND. .NOT. SaveFields ) THEN
      CALL Info('FCT_Correction','Skipping the computation of FCT correction',Level=5)
    END IF

    CALL Info('FCT_Correction','Computing corrector for the low order solution',Level=5)

    ! PRINT *,'FCT Norm Before Correction:',SQRT( SUM( Solver % Variable % Values**2) )
 
    Mesh => Solver % Mesh
    Variables => Solver % Mesh % Variables
 
    ! Set pointers 
    A => Solver % Matrix
    n = A % NumberOfRows
    Rows => A % Rows
    Cols => A % Cols
    BulkValues => A % BulkValues
    M_C => A % MassValues
    M_L => A % MassValuesLumped 
    Perm => Solver % Variable % Perm

    Var => VariableGet( Variables,'timestep size')
    dt = Var % Values(1) 

    ! low order solution at the start, high order in the end
    u => Solver % Variable % Values
    VarName = GetVarName( Solver % Variable ) 

    ! Here a bunch of vectors are stored for visualization and debugging purposes
    !----------------------------------------------------------------------------

    ! low order solution is the solution without corrections
    ! This is created and saved only if requested
    !---------------------------------------------------------------------------
    IF( SaveFields ) THEN
      TmpVarName = TRIM( VarName )//' fctlow'    
      Var => VariableGet( Variables, TmpVarName )
      IF( .NOT. ASSOCIATED(Var) ) THEN
        CALL VariableAddVector( Variables, Mesh, Solver,&
            TmpVarName, Perm = Perm, Output = SaveFields )
        Var => VariableGet( Variables, TmpVarName )
      END IF
      ulow => Var % Values
      ulow = u
    END IF

    ! Create auxiliary vectors for the fct algorithm
    !---------------------------------------------------------------------------

    ! r.h.s. term
    TmpVarName = TRIM( VarName )//' fctku'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    ku => Var % Values

    ! time derivative from lower order analysis
    TmpVarName = TRIM( VarName )//' fctudot'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    udot => Var % Values

    ! Fields related to flux limiters
    TmpVarName = TRIM( VarName )//' fctpp'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    pp => Var % Values
    
    TmpVarName = TRIM( VarName )//' fctpm'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    pm => Var % Values
    
    TmpVarName = TRIM( VarName )//' fctqp'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    qp => Var % Values

    TmpVarName = TRIM( VarName )//' fctqm'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    qm => Var % Values

    ! higher order correction 
    TmpVarName = TRIM( VarName )//' fctcorr'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    corr => Var % Values


    ! 1) Compute the nodal time derivatives
    ! M_C*udot=K*ulow  (M_C is the consistent mass matrix)
    !----------------------------------------------------------------------
    CALL Info('FCT_Correction','Compute nodal time derivatives',Level=10)
    ! Compute: ku = K*ulow
    DO i=1,n
      rsum = 0.0_dp
      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)
        K_ij = BulkValues(k)
        rsum = rsum + K_ij * u(j) 
      END DO
      ku(i) = rsum
    END DO
    ! Solve the linear system for udot
    ! The stiffness matrix is momentarily replaced by the consistent mass matrix M_C
    ! Also the namespace is replaced to 'fct:' so that different strategies may 
    ! be applied to the mass matrix solution.
    CALL ListPushNameSpace('fct:')
    SaveValues => A % Values
    A % Values => M_C
    CALL SolveLinearSystem( A, ku, udot, Norm, 1, Solver )
    A % Values => SaveValues
    CALL ListPopNamespace()

    ! Computation of correction factors (Zalesak's limiter)
    ! Code derived initially from Kuzmin's subroutine   
    !---------------------------------------------------------
    CALL Info('FCT_Correction','Compute correction factors',Level=10)
    pp = 0 
    pm = 0
    qp = 0 
    qm = 0

    DO i=1,n
      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)

        ! This is symmetric so lower section will do
        IF( i >= j .AND. Symmetric ) CYCLE
         
        ! Compute the raw antidiffusive fluxes
        ! f_ij=m_ij*[udot(i)-udot(j)]+d_ij*[ulow(i)-ulow(j)]
        !-----------------------------------------------------
        ! d_ij and m_ij are both symmetric
        ! Hence F_ji = -F_ij
        f_ij = M_C(k)*(udot(i)-udot(j)) + A % FCT_D(k) *(u(i)-u(j))
  
        ! Compared to Kuzmin's paper F_ij=-F_ij since d_ij and 
        ! udot have different signs. 
        f_ij = -f_ij

        ! Antidiffusive fluxes to be limited
        du = u(j)-u(i)

        ! Prelimiting of antidiffusive fluxes i.e. du and the flux have different signs
        IF (f_ij*du >= TINY(du)) THEN
          f_ij = 0
        ELSE        
          ! Positive/negative edge contributions
          pp(i) = pp(i)+MAX(0d0,f_ij)
          pm(i) = pm(i)+MIN(0d0,f_ij)
          ! symmetric part
          IF( Symmetric ) THEN
            pp(j) = pp(j)+MAX(0d0,-f_ij)
            pm(j) = pm(j)+MIN(0d0,-f_ij)
          END IF
        END IF

        ! Maximum/minimum solution increments
        qp(i) = MAX(qp(i),du)
        qm(i) = MIN(qm(i),du)
        ! symmetric part
        IF( Symmetric ) THEN
          qp(j) = MAX(qp(j),-du)
          qm(j) = MIN(qm(j),-du)
        END IF
      END DO
    END DO

    ! Computation of nodal correction factors
    ! These are eliminated as vectors to save some space
    ! and replaced by rpi, rpj, rmi, rmj
    ! DO i=1,n
    !  IF( pp(i) > Ceps ) THEN
    !    rp(i) = MIN( 1.0_dp, M_L(i)*qp(i)/pp(i) )
    !  END IF
    !  IF( pm(i) < -Ceps ) THEN
    !    rm(i) = MIN( 1.0_dp, M_L(i)*qm(i)/pm(i) )
    !  END IF
    ! END DO

    ! Correct the low-order solution
    ! (M_L*ufct)_i=(M_L*ulow)_i+dt*sum(alpha_ij*f_ij)
    !-------------------------------------------------
    ! Correction of the right-hand side
    ! Symmetric flux limiting


    Ceps = TINY( Ceps )
    corr = 0.0_dp
    DO i=1,n

      IF( pp(i) > Ceps ) THEN
        rpi = MIN( 1.0_dp, M_L(i)*qp(i)/pp(i) )
      ELSE
        rpi = 0.0_dp
      END IF

      IF( pm(i) < -Ceps ) THEN
        rmi = MIN( 1.0_dp, M_L(i)*qm(i)/pm(i) )
      ELSE
        rmi = 0.0_dp
      END IF

      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)
        IF( i >= j .AND. Symmetric ) CYCLE
               
        f_ij = M_C(k)*(udot(i)-udot(j)) + A % FCT_D(k) *(u(i)-u(j))
        f_ij = -f_ij

        IF (f_ij > 0) THEN 
          IF( pm(j) < -Ceps ) THEN
            rmj = MIN( 1.0_dp, M_L(j)*qm(j)/pm(j) )
          ELSE
            rmj = 0.0_dp
          END IF
          c_ij = MIN(rpi,rmj)
        ELSE 
          IF( pp(j) > Ceps ) THEN
            rpj = MIN( 1.0_dp, M_L(j)*qp(j)/pp(j) )
          ELSE
            rpj = 0.0_dp
          END IF
          c_ij = MIN(rmi,rpj)
        END IF

        corr(i) = corr(i) + c_ij * f_ij
        IF( Symmetric ) THEN
          corr(j) = corr(j) - c_ij * f_ij
        END IF
      END DO
    END DO

    CorrCoeff = ListGetCReal( Params,'FCT Correction Coefficient',Found )
    IF( .NOT. Found ) CorrCoeff = 1.0_dp

    ! This was suggestd in some code but results to invalida units, and poor results
    ! CorrCoeff = CorrCoeff * dt

    corr = CorrCoeff * corr / M_L

    ! Optinally skip applying the correction, just for debugging purposes
    IF( SkipCorrection ) THEN
      CALL Info('FCT_Correction','Skipping Applying corrector',Level=4)
    ELSE     
      CALL Info('FCT_Correction','Applying corrector for the low order solution',Level=10)
      u = u + corr
      ! PRINT *,'FCT Norm After Correction:',SQRT( SUM( Solver % Variable % Values**2) )
    END IF

  END SUBROUTINE FCT_Correction



  ! Create Linear constraints from mortar BCs:
  ! -------------------------------------------   
  SUBROUTINE GenerateProjectors(Model,Solver,Nonlinear,SteadyState) 
    
     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver
     LOGICAL, OPTIONAL :: Nonlinear, SteadyState

     LOGICAL :: IsNonlinear,IsSteadyState,Timing, RequireNonlinear, ContactBC
     LOGICAL :: ApplyMortar, ApplyContact, Found
     INTEGER :: i,j,k,l,n,dsize,size0,col,row,dim
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Matrix_t), POINTER :: CM, CMP, CM0, CM1
     TYPE(Variable_t), POINTER :: DispVar
     REAL(KIND=dp) :: t0,rt0,rst,st,ct
#ifndef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: CPUTime,RealTime
#endif

     ApplyMortar = ListGetLogical(Solver % Values,'Apply Mortar BCs',Found) 
     ApplyContact = ListGetLogical(Solver % Values,'Apply Contact BCs',Found) 

     IF( .NOT. ( ApplyMortar .OR. ApplyContact) ) RETURN
     
     i = ListGetInteger( Solver % Values,'Mortar BC Master Solver',Found ) 
     IF( Found ) THEN
       Solver % MortarBCs => CurrentModel % Solvers(i) % MortarBCs
       IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) THEN
         CALL Fatal('GenerateProjectors','Could not reuse projectors from solver: '//TRIM(I2S(i)))
       END IF
       CALL Info('GenerateProjectors','Reusing projectors from solver: '//TRIM(I2S(i)),Level=8)
       RETURN
     END IF

     CALL Info('GenerateProjectors','Generating mortar projectors',Level=8)

     Timing = ListCheckPrefix(Solver % Values,'Projector Timing')
     IF( Timing ) THEN
       t0 = CPUTime(); rt0 = RealTime()      
     END IF

     IsNonlinear = .FALSE.
     IF( PRESENT( Nonlinear ) ) IsNonlinear = Nonlinear
     IsSteadyState = .NOT. IsNonlinear

     IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) THEN
       ALLOCATE( Solver % MortarBCs( Model % NumberOfBCs ) )
       DO i=1, Model % NumberOfBCs
         Solver % MortarBCs(i) % Projector => NULL()
       END DO
     END IF
     
     dim = CoordinateSystemDimension()

     DO i=1,Model % NumberOFBCs
       BC => Model % BCs(i) % Values
       
       ContactBC = .FALSE.
       j = ListGetInteger( BC,'Mortar BC',Found)       
       IF( .NOT. Found ) THEN
         j = ListGetInteger( BC,'Contact BC',Found)       
         ContactBC = Found
       END IF
       IF( .NOT. Found ) CYCLE

       RequireNonlinear = ListGetLogical( BC,'Mortar BC Nonlinear',Found)
       IF( .NOT. Found ) THEN
         RequireNonlinear = ContactBC .AND. .NOT. ListGetLogical( BC,'Tie Contact',Found )
       END IF

       IF( IsNonlinear ) THEN
         IF( .NOT. RequireNonlinear ) CYCLE
       ELSE
         IF( RequireNonlinear ) CYCLE
       END IF             

       IF( ASSOCIATED( Solver % MortarBCs(i) % Projector ) ) THEN
         IF( ListGetLogical( BC,'Mortar BC Static',Found) ) CYCLE         
         
         IF( ASSOCIATED( Solver % MortarBCs(i) % Projector % Ematrix ) ) THEN
           CALL FreeMatrix( Solver % MortarBCs(i) % Projector % Ematrix )
         END IF
         CALL FreeMatrix( Solver % MortarBCs(i) % Projector )
       END IF
       
       Solver % MortarBCs(i) % Projector => &
           PeriodicProjector(Model,Solver % Mesh,i,j,dim,.TRUE.)
       
       IF( ASSOCIATED( Solver % MortarBCs(i) % Projector ) ) THEN
         Solver % MortarBCsChanged = .TRUE.
       END IF

     END DO


     IF( Timing ) THEN
       st  = CPUTime() - t0;
       rst = RealTime() - rt0
       
       WRITE(Message,'(a,f8.2,f8.2,a)') 'Projector creation time (CPU,REAL) for '&
           //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
       CALL Info('GenerateProjectors',Message)    
       
       IF( ListGetLogical(Solver % Values,'Projector Timing',Found)) THEN
         CALL ListAddConstReal(CurrentModel % Simulation,'res: projector cpu time '&
             //GetVarName(Solver % Variable),st)
         CALL ListAddConstReal(CurrentModel % Simulation,'res: projector real time '&
             //GetVarName(Solver % Variable),rst)
       END IF

       IF( ListGetLogical(Solver % Values,'Projector Timing Cumulative',Found)) THEN
         ct = ListGetConstReal(CurrentModel % Simulation,'res: cum projector cpu time '&
             //GetVarName(Solver % Variable),Found)
         st = st + ct
         ct = ListGetConstReal(CurrentModel % Simulation,'res: cum projector real time '&
             //GetVarName(Solver % Variable),Found)
         rst = rst + ct
         CALL ListAddConstReal(CurrentModel % Simulation,'res: cum projector cpu time '&
             //GetVarName(Solver % Variable),st)
         CALL ListAddConstReal(CurrentModel % Simulation,'res: cum projector real time '&
             //GetVarName(Solver % Variable),rst)
       END IF
     END IF
     
   END SUBROUTINE GenerateProjectors



   ! Generate constraint matrix from mortar projectors. 
   ! This routine takes each boundary projector and applies it 
   ! to the current field variable (scalar or vector) merging 
   ! all into one single projector. 
   !---------------------------------------------------------
   SUBROUTINE GenerateConstraintMatrix( Model, Solver )

     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver

     INTEGER, POINTER :: Perm(:)
     INTEGER :: i,j,j2,k,k2,dofs,maxperm,permsize,bc_ind,constraint_ind,row,col,col2,mcount,bcount
     TYPE(Matrix_t), POINTER :: Atmp,Btmp, Ctmp
     LOGICAL :: AllocationsDone, CreateSelf, ComplexMatrix, TransposePresent, Found, &
         SetDof, SomeSet, SomeSkip, SumProjectors, NewRow
     INTEGER, ALLOCATABLE :: SumPerm(:),SumCount(:)
     LOGICAL, ALLOCATABLE :: ActiveComponents(:), SetDefined(:)
     TYPE(ValueList_t), POINTER :: BC
     TYPE(MortarBC_t), POINTER :: MortarBC
     REAL(KIND=dp) :: wsum, Scale
     INTEGER :: rowoffset, arows, sumrow, EliminatedRows, NeglectedRows
     CHARACTER(LEN=MAX_NAME_LEN) :: Str
     LOGICAL :: AnyPriority
     INTEGER :: Priority, PrevPriority
     INTEGER, ALLOCATABLE :: BCOrdering(:), BCPriority(:)


     CALL Info('GenerateConstraintMatrix','Building constraint matrix',Level=12)

     IF( Solver % MortarBCsOnly .AND. .NOT. Solver % MortarBCsChanged ) THEN
       CALL Info('GenerateConstraintMatrix','Nothing to do!',Level=12)
       RETURN
     ELSE IF ( Solver % MortarBCsOnly ) THEN
      CALL Info('GenerateConstraintMatrix','Releasing constraint matrix',Level=15)
      CALL ReleaseConstraintMatrix(Solver)
     END IF
     
     ! Compute the size of the initial boundary matrices.
     !------------------------------------------------------
     row    = 0
     mcount = 0
     bcount = 0
     IF(.NOT. Solver % MortarBcsOnly) THEN
       Ctmp => Solver % Matrix % ConstraintMatrix
       DO WHILE(ASSOCIATED(Ctmp))
         mcount = mcount + 1
         row = row + Ctmp % NumberOfRows
       END DO
     END IF

     IF(Solver % MortarBCsChanged) THEN
       DO bc_ind=1,Model % NumberOFBCs
         Atmp => Solver % MortarBCs(bc_ind) % Projector
         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE
         bcount = bcount + 1
         row = row + Atmp % NumberOfRows
       END DO
     END IF

     IF( row==0 .OR. bcount==0 .AND. mcount<=1 )  RETURN
     
     SumProjectors = ListGetLogical( Solver % Values,&
         'Mortar BCs Additive', Found )
     IF( .NOT. Found .AND. bcount > 1 ) THEN
       IF( ListGetLogical( Solver % Values, &
           'Eliminate Linear Constraints',Found ) ) THEN
         CALL Info('GenerateConstraintMatrix',&
             'Enforcing > Mortar BCs Additive < to True to enable elimination',Level=6)
         SumProjectors = .TRUE.
       END IF
     END IF
     EliminatedRows = 0

     CALL Info('GenerateConstraintMatrix','There are '&
         //TRIM(I2S(row))//' initial rows in constraint matrices',Level=10)
     
     dofs = Solver % Variable % DOFs
     Perm => Solver % Variable % Perm
     permsize = SIZE( Perm )
     maxperm  = MAXVAL( Perm )
     AllocationsDone = .FALSE.
     arows = Solver % Matrix % NumberOfRows
     
     ALLOCATE( ActiveComponents(dofs), SetDefined(dofs) ) 
     
     IF( SumProjectors ) THEN
       ALLOCATE( SumPerm( dofs * permsize ) )
       SumPerm = 0
       ALLOCATE( SumCount( arows ) )
       SumCount = 0
     END IF

     ComplexMatrix = Solver % Matrix % Complex
     IF( ComplexMatrix ) THEN
       IF( MODULO( Dofs,2 ) /= 0 ) CALL Fatal('GenerateConstraintMatrix',&
           'Complex matrix should have even number of components!')
     ELSE
       ! Currently complex matrix is enforced if there is an even number of 
       ! entries since it seems that we cannot rely on the flag to be set.
       ComplexMatrix = ListGetLogical( Solver % Values,'Linear System Complex',Found )
       IF( .NOT. Found ) ComplexMatrix = ( MODULO( Dofs,2 ) == 0 )
     END IF

     
     AnyPriority = ListCheckPresentAnyBC( Model,'Projector Priority') 
     IF( AnyPriority ) THEN
       IF(.NOT. SumProjectors ) THEN
         CALL Warn('GenerateConstraintMatrix','Priority has effect only in additive mode!')
         AnyPriority = .FALSE.
       ELSE
         CALL Info('GenerateConstraintMatrix','Using priority for projector entries',Level=7)
         ALLOCATE( BCPriority(Model % NumberOfBCs), BCOrdering( Model % NumberOfBCs) )
         BCPriority = 0; BCOrdering = 0
         DO bc_ind=1, Model % NumberOFBCs
           Priority = ListGetInteger( Model % BCs(bc_ind) % Values,'Projector Priority',Found)
           BCPriority(bc_ind) = -bc_ind + Priority * Model % NumberOfBCs 
           BCOrdering(bc_ind) = bc_ind
         END DO
         CALL SortI( Model % NumberOfBCs, BCPriority, BCOrdering )
       END IF
     END IF
     NeglectedRows = 0


100  sumrow = 0
     k2 = 0
     rowoffset = 0
     Priority = -1
     PrevPriority = -1
    
     TransposePresent = .FALSE.
     Ctmp => Solver % Matrix % ConstraintMatrix
     DO constraint_ind=Model % NumberOFBCs+mcount,1,-1
       
       ! This is the default i.e. all components are applied mortar BCs
       ActiveComponents = .TRUE.
       
       IF(constraint_ind>Model % NumberOfBCs) THEN
         Atmp => Ctmp
         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE
         Ctmp => Ctmp % ConstraintMatrix
       ELSE
         IF(.NOT. Solver % MortarBCsChanged) EXIT
         
         IF( AnyPriority ) THEN
           bc_ind = BCOrdering(constraint_ind)
         ELSE
           bc_ind = constraint_ind 
         END IF

         MortarBC => Solver % MortarBCs(bc_ind) 
         Atmp => MortarBC % Projector
         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE

         CALL Info('GenerateConstraintMatrix','Adding projector for BC: '//TRIM(I2S(bc_ind)),Level=8)

         IF( .NOT. ASSOCIATED( Atmp % InvPerm ) ) THEN
           CALL Fatal('GenerateConstraintMatrix','InvPerm is required!')
         END IF

         Priority = ListGetInteger( Model % BCs(bc_ind) % Values,'Projector Priority',Found)

         ! Enable that the user can for vector valued cases either set some 
         ! or skip some field components. 
         SomeSet = .FALSE.
         SomeSkip = .FALSE.
         DO i=1,Dofs
           str = ComponentNameVar( Solver % Variable, i )
           SetDof = ListGetLogical( Model % BCs(bc_ind) % Values,'Mortar BC '//TRIM(str),Found )
           SetDefined(i) = Found
           IF(Found) THEN
             ActiveComponents(i) = SetDof
             IF( SetDof ) THEN
               SomeSet = .TRUE.
             ELSE
               SomeSkip = .TRUE.
             END IF
           END IF
         END DO
         
         ! By default all components are applied mortar BC and some are turned off.
         ! If the user does the opposite then the default for other components is True. 
         IF( SomeSet .AND. .NOT. ALL(SetDefined) ) THEN
           IF( SomeSkip ) THEN
             CALL Fatal('GenerateConstraintMatrix','Dont know what to do with all components')
           ELSE 
             CALL Info('GenerateConstraintMatrix',&
                 'Unspecified components will not be set for BC '//TRIM(I2S(bc_ind)),Level=10)
             DO i=1,Dofs
               IF( .NOT. SetDefined(i) ) ActiveComponents(i) = .FALSE.
             END DO
           END IF
         END IF

       END IF
       TransposePresent = TransposePresent .OR. ASSOCIATED(Atmp % Child)
       IF( TransposePresent ) THEN
         CALL Info('GenerateConstraintMatrix','Transpose matrix is present')
       END IF

       ! If the projector is of type x_s=P*x_m then generate a constraint matrix
       ! of type [D-P]x=0.
       CreateSelf = ( Atmp % ProjectorType == PROJECTOR_TYPE_NODAL ) 
       
       IF( SumProjectors .AND. CreateSelf ) THEN
         CALL Fatal('GenerateConstraintMatrix','It is impossible to sum up nodal projectors!')
       END IF
       
       IF( Dofs == 1 ) THEN         

         IF( .NOT. ActiveComponents(1) ) CYCLE

         DO i=1,Atmp % NumberOfRows           

           IF( Atmp % Rows(i) >= Atmp % Rows(i+1) ) CYCLE ! skip empty rows

           ! If the mortar boundary is not active at this round don't apply it
           IF( ASSOCIATED( MortarBC % Active ) ) THEN
             IF( .NOT. MortarBC % Active(i) ) CYCLE
           END IF
           
           ! Node does not have an active dof to be constrained
           k = Atmp % InvPerm(i)
           IF( k == 0 ) CYCLE
           IF( Perm(k) == 0 ) CYCLE

           IF( SumProjectors ) THEN
             NewRow = ( SumPerm(k) == 0 )
             IF( NewRow ) THEN
               sumrow = sumrow + 1                
               IF( Priority /= 0 ) THEN
                 SumPerm(k) = -sumrow
               ELSE
                 SumPerm(k) = sumrow 
               END IF
             ELSE IF( Priority /= PrevPriority .AND. SumPerm(k) < 0 ) THEN
               IF(.NOT. AllocationsDone ) THEN
                 NeglectedRows = NeglectedRows + 1
                 !PRINT *,'Neglecting:',bc_ind,k
               END IF
               CYCLE
             ELSE
               IF(.NOT. AllocationsDone ) THEN
                 EliminatedRows = EliminatedRows + 1
                 !PRINT *,'Eliminating:',bc_ind,k
               END IF
             END IF
             row = ABS( SumPerm(k) )
           ELSE
             sumrow = sumrow + 1
             row = sumrow
           END IF
           
           IF( AllocationsDone ) THEN
             Btmp % InvPerm(row) = rowoffset + Perm( k ) 
           END IF

           wsum = 0.0_dp

           DO k=Atmp % Rows(i),Atmp % Rows(i+1)-1
             
             col = Atmp % Cols(k) 
             
             IF( col <= permsize ) THEN
               col2 = Perm(col)
               IF( col2 == 0 ) CYCLE
             ELSE 
               PRINT *,'col too large',col,permsize
               CYCLE
             END IF
             
             IF( AllocationsDone ) THEN
               IF( CreateSelf ) THEN
                 ! We want to create [D-P] hence the negative sign
                 Scale = MortarBC % MasterScale
                 wsum = wsum + Atmp % Values(k)
               ELSE IF( ASSOCIATED( MortarBC % Perm ) ) THEN
                 ! Look if the component refers to the slave
                 IF( MortarBC % Perm( col ) > 0 ) THEN
                   Scale = MortarBC % SlaveScale 
                   wsum = wsum + Atmp % Values(k) 
                 ELSE
                   Scale = MortarBC % MasterScale
                 END IF
               ELSE
                 ! Ok, we don't have any complex physics, just use the scaling as it is
                 Scale = 1.0_dp
               END IF

               ! Add a new column index to the summed up row
               IF( SumProjectors ) THEN
                 k2 = Btmp % Rows(row)
                 DO WHILE( Btmp % Cols(k2) > 0 )
                   k2 = k2 + 1
                 END DO
               ELSE
                 k2 = k2 + 1
               END IF

               Btmp % Cols(k2) = col2
               Btmp % Values(k2) = Scale * Atmp % Values(k)
               IF(ASSOCIATED(Btmp % TValues)) THEN
                 IF(ASSOCIATED(Atmp % Child)) THEN
                   Btmp % TValues(k2) = Scale * Atmp % Child % Values(k)
                 ELSE
                   Btmp % TValues(k2) = Scale * Atmp % Values(k)
                 END IF
               END IF
             ELSE
               k2 = k2 + 1
               IF( SumProjectors ) THEN
                 SumCount(row) = SumCount(row) + 1
               END IF
             END IF
           END DO
           
           ! Add the self entry as in 'D'
           IF( CreateSelf ) THEN
             k2 = k2 + 1
             IF( AllocationsDone ) THEN
               Btmp % Cols(k2) = Perm( Atmp % InvPerm(i) )
               Btmp % Values(k2) = MortarBC % SlaveScale * wsum
             END IF
           END IF


           ! Add a diagonal entry if requested. When this is done at the final stage
           ! all the hazzle with the right column index is easier. 
           IF( ASSOCIATED( MortarBC % Diag ) ) THEN
             IF( MortarBC % LumpedDiag ) THEN
               k2 = k2 + 1             
               IF( AllocationsDone ) THEN
                 Btmp % Cols(k2) = row + arows 
                 ! The factor 0.5 comes from the fact that the 
                 ! contribution is summed twice, 2nd time as transpose
                 ! For Nodal projector the entry is 1/(weight*coeff)
                 ! For Galerkin projector the is weight/coeff 
                 Btmp % Values(k2) = -0.5_dp * MortarBC % Diag(i) * wsum
               END IF
             ELSE
               DO k=Atmp % Rows(i),Atmp % Rows(i+1)-1                 
                 col = Atmp % Cols(k) 
                 
                 IF( col > permsize ) CYCLE
                 col2 = Perm(col)

                 IF( CreateSelf ) THEN
                   Scale = -MortarBC % MasterScale
                 ELSE 
                   IF( MortarBC % Perm( col ) > 0 ) THEN
                     Scale = MortarBC % SlaveScale 
                   ELSE
                     CYCLE                     
                   END IF
                 END IF

                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN                   
                   Btmp % Cols(k2) = MortarBC % Perm( col ) + arows + rowoffset
                   Btmp % Values(k2) = -0.5_dp * Atmp % Values(k) * &
                       MortarBC % Diag(i) 
                 END IF
               END DO
              
             END IF
           END IF

           IF( AllocationsDone ) THEN
             IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
               Btmp % Rhs(row) = Btmp % Rhs(row) + wsum * MortarBC % rhs(i)
             END IF
             IF( .NOT. SumProjectors ) THEN
               Btmp % Rows(row+1) = k2 + 1
             END IF
           END IF
         END DO
       ELSE ! dofs > 1
         ! In case of a vector valued problem create a projector that acts on all 
         ! components of the vector. Otherwise follow the same logic.
         DO i=1,Atmp % NumberOfRows           
           DO j=1,Dofs
             
             IF( .NOT. ActiveComponents(j) ) CYCLE

             ! For complex matrices both entries mist be created
             ! since preconditioning benefits from 
             IF( ComplexMatrix ) THEN
               IF( MODULO( j, 2 ) == 0 ) THEN
                 j2 = j-1
               ELSE 
                 j2 = j+1
               END IF
             ELSE
               j2 = 0
             END IF

             IF( ASSOCIATED( MortarBC % Active ) ) THEN
               IF( .NOT. MortarBC % Active(Dofs*(i-1)+j) ) CYCLE
             END IF
                          
             k = Atmp % InvPerm(i)
             IF( k == 0 ) CYCLE
             IF( Perm(k) == 0 ) CYCLE

             IF( SumProjectors ) THEN
               NewRow = ( SumPerm(Dofs*(k-1)+j) == 0 )
               IF( NewRow ) THEN
                 sumrow = sumrow + 1                
                 IF( Priority /= 0 ) THEN
                   ! Use negative sign to show that this has already been set by priority
                   SumPerm(Dofs*(k-1)+j) = -sumrow 
                 ELSE
                   SumPerm(Dofs*(k-1)+j) = sumrow 
                 END IF
               ELSE IF( Priority /= PrevPriority .AND. SumPerm(Dofs*(k-1)+j) < 0 ) THEN
                 IF(.NOT. AllocationsDone ) THEN
                   NeglectedRows = NeglectedRows + 1
                 END IF                 
                 CYCLE
               ELSE
                 IF(.NOT. AllocationsDone ) THEN
                   EliminatedRows = EliminatedRows + 1
                 END IF
               END IF
               row = ABS( SumPerm(Dofs*(k-1)+j) )
             ELSE
               sumrow = sumrow + 1
               row = sumrow
             END IF

             IF( AllocationsDone ) THEN
               Btmp % InvPerm(row) = rowoffset + Dofs * ( Perm( k ) - 1 ) + j
             END IF

             wsum = 0.0_dp

             DO k=Atmp % Rows(i),Atmp % Rows(i+1)-1             

               col = Atmp % Cols(k)                
               
               IF( col <= permsize ) THEN
                 col2 = Perm(col)
                 IF( col2 == 0 ) CYCLE
               ELSE 
                 PRINT *,'col too large',col,permsize
                 CYCLE
               END IF              
               
               IF( AllocationsDone ) THEN
                 
                 IF( CreateSelf ) THEN
                   Scale = MortarBC % MasterScale
                   wsum = wsum + Atmp % Values(k)
                 ELSE IF( ASSOCIATED( MortarBC % Perm ) ) THEN
                   IF( MortarBC % Perm(col) > 0 ) THEN
                     Scale = MortarBC % SlaveScale 
                     wsum = wsum + Atmp % Values(k) 
                   ELSE
                     Scale = MortarBC % MasterScale
                   END IF
                 ELSE
                   Scale = 1.0_dp
                 END IF

                 IF( SumProjectors ) THEN
                   k2 = Btmp % Rows(row)
                   DO WHILE( Btmp % Cols(k2) > 0 )
                     k2 = k2 + 1
                   END DO
                 ELSE
                   k2 = k2 + 1
                 END IF
                 
                 Btmp % Cols(k2) = Dofs * ( col2 - 1) + j
                 Btmp % Values(k2) = Scale * Atmp % Values(k)
                 IF(ASSOCIATED(Btmp % Tvalues)) THEN
                   IF(ASSOCIATED(Atmp % Child)) THEN
                     Btmp % TValues(k2) = Scale * Atmp % Child % Values(k)
                   ELSE
                     Btmp % TValues(k2) = Scale * Atmp % Values(k)
                   END IF
                 END IF
               ELSE
                 k2 = k2 + 1 
                 IF( SumProjectors ) THEN
                   SumCount(row) = SumCount(row) + 1
                 END IF                 
               END IF
             END DO
              
             ! Add the self entry as in 'D'
             IF( CreateSelf ) THEN
               k2 = k2 + 1
               IF( AllocationsDone ) THEN
                 Btmp % Cols(k2) = Dofs * ( Perm( Atmp % InvPerm(i) ) -1 ) + j
                 Btmp % Values(k2) = MortarBC % SlaveScale * wsum
               END IF
             END IF
             
             ! Create the imaginary part (real part) corresponding to the 
             ! real part (imaginary part) of the projector. 
             IF( j2 /= 0 ) THEN
               DO k=Atmp % Rows(i),Atmp % Rows(i+1)-1             

                 col = Atmp % Cols(k)                
                 
                 IF( col <= permsize ) THEN
                   col2 = Perm(col)
                   IF( col2 == 0 ) CYCLE
                 END IF
                 
                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   Btmp % Cols(k2) = Dofs * ( col2 - 1) + j2
                 END IF
               END DO

               IF( CreateSelf ) THEN
                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   Btmp % Cols(k2) = Dofs * ( Perm( Atmp % InvPerm(i) ) -1 ) + j2
                 END IF
               END IF
             END IF

            
             IF( ASSOCIATED( MortarBC % Diag ) ) THEN
               IF( MortarBC % LumpedDiag ) THEN
                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   Btmp % Cols(k2) = row + arows
                   Btmp % Values(k2) = -0.5_dp * wsum * MortarBC % Diag(Dofs*(i-1)+j)
                 END IF
               ELSE
                 DO k=Atmp % Rows(i),Atmp % Rows(i+1)-1                 
                   col = Atmp % Cols(k) 
                   
                   IF( col > permsize ) CYCLE
                   col2 = Perm(col)
                   
                   IF( CreateSelf ) THEN
                     Scale = -MortarBC % MasterScale
                   ELSE 
                     IF( MortarBC % Perm( col ) > 0 ) THEN
                       Scale = MortarBC % SlaveScale 
                     ELSE
                       CYCLE                     
                     END IF
                   END IF
                   
                   k2 = k2 + 1
                   IF( AllocationsDone ) THEN                   
                     Btmp % Cols(k2) = Dofs*(MortarBC % Perm( col )-1)+j + arows + rowoffset
                     Btmp % Values(k2) = -0.5_dp * Atmp % Values(k) * &
                         MortarBC % Diag(i) 
                   END IF
                 END DO

               END IF
             END IF

               
             IF( AllocationsDone ) THEN
               IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
                 Btmp % Rhs(row) = wsum * MortarBC % rhs(Dofs*(i-1)+j)
               END IF
               IF(.NOT. SumProjectors ) THEN
                 Btmp % Rows(row+1) = k2 + 1
               END IF
             END IF

           END DO
         END DO
       END IF ! dofs > 1
       
       IF( .NOT. SumProjectors ) THEN
         rowoffset = rowoffset + Arows
       END IF

       PrevPriority = Priority 
     END DO

     IF( k2 == 0 ) THEN
       CALL Info('GenerateConstraintMatrix','No entries in constraint matrix!',Level=6)
!      Solver % Matrix % ConstraintMatrix => NULL()
       RETURN
     END IF

     ! Allocate the united matrix of all the boundary matrices
     !-------------------------------------------------------
     IF( .NOT. AllocationsDone ) THEN
       CALL Info('GenerateConstraintMatrix','Allocating '//&
           TRIM(I2S(sumrow))//' rows and '//TRIM(I2S(k2))//' nonzeros',&
           Level=6)

       Btmp => AllocateMatrix()
       ALLOCATE( Btmp % RHS(sumrow), Btmp % Rows(sumrow+1), &
           Btmp % Cols(k2), Btmp % Values(k2), &
           Btmp % InvPerm(sumrow) )

       Btmp % Rhs = 0.0_dp
       Btmp % Rows = 0
       Btmp % Cols = 0
       Btmp % Values = 0.0_dp
       Btmp % NumberOFRows = sumrow 
       Btmp % InvPerm = 0
       Btmp % Rows(1) = 1

       IF(TransposePresent) THEN
         ALLOCATE(Btmp % TValues(k2))
         Btmp % Tvalues = 0._dp
       END IF

       IF( SumProjectors ) THEN
         DO i=2,sumrow+1
           Btmp % Rows(i) = Btmp % Rows(i-1) + SumCount(i-1)
         END DO
         SumPerm = 0
         DEALLOCATE( SumCount ) 
       END IF

       AllocationsDone = .TRUE.

       GOTO 100
     END IF

     ! Eliminate entries
     IF( EliminatedRows > 0 ) THEN
       CALL Info('GenerateConstraintMatrix','Number of eliminated rows: '&
           //TRIM(I2S(EliminatedRows)))
       CALL CRS_PackMatrix( Btmp ) 
     END IF

     IF( NeglectedRows > 0 ) THEN
       CALL Info('GenerateConstraintMatrix','Number of neglected rows: '&
           //TRIM(I2S(NeglectedRows)))
     END IF

     Solver % Matrix % ConstraintMatrix => Btmp
     
     Solver % MortarBCsChanged = .FALSE.

     CALL Info('GenerateConstraintMatrix','Finished creating constraint matrix',Level=12)

   END SUBROUTINE GenerateConstraintMatrix
     

   SUBROUTINE ReleaseConstraintMatrix(Solver) 
     TYPE(Solver_t) :: Solver

     CALL FreeMatrix(Solver % Matrix % ConstraintMatrix)
     Solver % Matrix % ConstraintMatrix => Null()

   END SUBROUTINE ReleaseConstraintMatrix


   SUBROUTINE ReleaseProjectors(Model, Solver) 

     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver

     TYPE(ValueList_t), POINTER :: BC
     TYPE(Matrix_t), POINTER :: Projector
     INTEGER :: i
     

     IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) RETURN

     DO i=1,Model % NumberOFBCs
       BC => Model % BCs(i) % Values
       Projector => Solver % MortarBCs(i) % Projector 
       IF( ASSOCIATED( Projector ) ) THEN
         IF( ASSOCIATED( Projector % EMatrix ) ) THEN
           CALL FreeMatrix( Projector % Ematrix ) 
         END IF
         CALL FreeMatrix( Projector )
         Solver % MortarBCs(i) % Projector => NULL()
       END IF
     END DO

   END SUBROUTINE ReleaseProjectors


END MODULE SolverUtils

!> \}
