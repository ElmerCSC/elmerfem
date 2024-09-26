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
   USE PElementMaps, ONLY : isActivePElement, getEdgeDOFs, getFaceDOFs, getBubbleDOFs

   
   IMPLICIT NONE

   INTERFACE CondensateP
     MODULE PROCEDURE CondensatePR, CondensatePC
   END INTERFACE CondensateP

CONTAINS
   
!> Sets the matrix element to a desired value. 
!------------------------------------------------------------------------------
   SUBROUTINE SetMatrixElement( A, i, j, val )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A  !< Structure holding the matrix
     INTEGER :: i                            !< Row index
     INTEGER :: j                            !< Column index
     REAL(KIND=dp) :: val                  !< Value to be obtained
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_SetMatrixElement( A, i, j, val )
         IF(A % FORMAT == MATRIX_LIST) THEN
           CALL List_toListMatrix(A)
           CALL List_SetMatrixElement( A % ListMatrix, i, j, val )
         END IF

       CASE( MATRIX_LIST )
         CALL List_SetMatrixElement( A % ListMatrix, i, j, val )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_SetMatrixElement( A, i, j, val )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE SetMatrixElement
!------------------------------------------------------------------------------

!> Gets a matrix element. 
!------------------------------------------------------------------------------
   FUNCTION GetMatrixElement( A, i, j ) RESULT ( val )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A  !< Structure holding the matrix
     INTEGER :: i                            !< Row index
     INTEGER :: j                            !< Column index
     REAL(KIND=dp) :: val                  !< Value to be obtained
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         val = CRS_GetMatrixElement( A, i, j )

      CASE( MATRIX_LIST )
         val = List_GetMatrixElement( A % ListMatrix, i, j )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         val = Band_GetMatrixElement( A, i, j )
     END SELECT
!------------------------------------------------------------------------------
   END FUNCTION GetMatrixElement
!------------------------------------------------------------------------------

!> Changes the value of a given matrix element.
!------------------------------------------------------------------------------
   FUNCTION ChangeMatrixElement( A, i, j, NewVal ) RESULT ( OldVal )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i,j
     REAL(KIND=dp) :: NewVal, OldVal
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         OldVal = CRS_ChangeMatrixElement( A, i, j, NewVal )

       CASE DEFAULT
         CALL Warn('ChangeMatrixElement','Not implemented for this type')

     END SELECT
!------------------------------------------------------------------------------
   END FUNCTION ChangeMatrixElement
!------------------------------------------------------------------------------


!> Adds to the value of a given matrix element.
!------------------------------------------------------------------------------
   SUBROUTINE AddToMatrixElement( A, i, j,val )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i,j
     REAL(KIND=dp) :: val
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_AddToMatrixElement( A, i, j, val )
         IF(A % FORMAT == MATRIX_LIST) THEN
           CALL List_toListMatrix(A)
           CALL List_AddToMatrixElement( A % ListMatrix, i, j, val )
         END IF

      CASE( MATRIX_LIST )
         CALL List_AddToMatrixElement( A % ListMatrix, i, j, val )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_AddToMatrixElement( A, i, j, val )
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
     REAL(KIND=dp) :: val

     val = ChangeMatrixElement(A, i1, j1, 0.0_dp)
     CALL AddToMatrixElement(A, i2, j2, val )
     
!------------------------------------------------------------------------------
   END SUBROUTINE MoveMatrixElement
!------------------------------------------------------------------------------


!> Zeros a row in matrix.
!------------------------------------------------------------------------------
   SUBROUTINE ZeroRow( A, n )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A  !< Structure holding the matrix 
     INTEGER :: n         !< Row to be zeroed.
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
           IF( PerFlip(Indexes(i)) .NEQV. PerFlip(Indexes(j)) ) THEN
             A(i,j) = -A(i,j)
           END IF
         END DO
       END DO
     ELSE
       DO i=1,n
         DO j=1,n
           IF( PerFlip(Indexes(i)) .NEQV. PerFlip(Indexes(j)) ) THEN
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


!---------------------------------------------------------------------------
!> Skip matrix assembly because the elements are similar.
!---------------------------------------------------------------------------
   FUNCTION UseLocalMatrixCopy( Solver, elemind, activeind) RESULT ( Skip ) 
     TYPE(Solver_t) :: Solver
     INTEGER, OPTIONAL :: elemind, activeind
     LOGICAL :: Skip 
     
     INTEGER :: eind, vind
     
     Skip = .FALSE.
     IF( Solver % LocalSystemMode <= 0 ) RETURN
     
     IF( PRESENT(activeind) ) THEN
       eind = activeind
     ELSE
       IF(eind > Solver % NumberOfActiveElements ) RETURN
       eind = Solver % InvActiveElements(elemind)
       IF(eind==0) RETURN
     END IF

     vind = Solver % LocalSystem(eind) % eind
     IF( vind > 0 .AND. vind /= eind) Skip = .TRUE.
     
   END FUNCTION UseLocalMatrixCopy
     

   
!---------------------------------------------------------------------------
!> Store local matrix, e.g. for topology optimization.
!---------------------------------------------------------------------------
   SUBROUTINE UseLocalMatrixStorage( Solver, n, K, F, elemind, activeind) 
     TYPE(Solver_t) :: Solver
     INTEGER :: n
     REAL(KIND=dp) :: K(:,:), F(:)
     INTEGER, OPTIONAL :: elemind, activeind
     
     TYPE(LocalSystemStorage_t), POINTER :: pLocal
     TYPE(Variable_t), POINTER :: cVar => NULL(), rVar => NULL()
     CHARACTER(:), ALLOCATABLE :: multname
     REAL(KIND=dp) :: cmult, rmult
     INTEGER :: prevSolverId = -1
     INTEGER :: eind, vind
     LOGICAL :: Found
     LOGICAL :: DoMultiply = .FALSE., DoMultiplyRhs = .FALSE.
     
     SAVE cVar, DoMultiply, DoMultiplyRhs, prevSolverId
     
     IF( PRESENT(activeind) ) THEN
       eind = activeind
     ELSE IF( PRESENT( elemind ) ) THEN
       ! This is probably related to boundary elements for which we haven't saved stuff. 
       IF( elemind > SIZE(Solver % InvActiveElements ) ) RETURN
       !IF( elemind < 1 .OR. elemind > SIZE(Solver % InvActiveElements ) ) THEN
       !  CALL Fatal('GetLocalMatrixStorage','Parameter "elemind" is out of bounds!')
       !END IF
       ! Given element index, return its position in the active element list.
       eind = Solver % InvActiveElements(elemind)
       IF(eind==0) RETURN
     ELSE
       CALL Fatal('UseLocalMatrixStorage','Give element index as a parameter in some way!')
     END IF

     ! Size of local system is same as size of ActiveElements table. 
     pLocal => Solver % LocalSystem(eind) 
     IF( pLocal % eind == eind .OR. pLocal % eind < 1 ) THEN
       ! Save local system for this element.
       IF(pLocal % n < n ) THEN
         ! The size of local system allocated. 
         IF( plocal % n > 0 ) THEN
           DEALLOCATE(pLocal % K, pLocal % F) 
         END IF
         pLocal % n = n
         ALLOCATE(pLocal % K(n,n), pLocal % F(n)) 
       END IF
       ! This elements is saved in its location. 
       pLocal % eind = eind
       pLocal % K(1:n,1:n) = K(1:n,1:n)
       pLocal % F(1:n) = F(1:n)
     ELSE
       ! Obtain local system for this element which is copy of some other element
       pLocal => Solver % LocalSystem(pLocal % eind)
       K(1:n,1:n) = pLocal % K(1:n,1:n) 
       F(1:n) = pLocal % F(1:n)        
     END IF

     IF(Solver % SolverId /= PrevSolverId) THEN
       ! For the 1st element obtain the multiplier vector, for other elements it will be the same. 
       cvar => NULL()
       DoMultiply = .FALSE.
       DoMultiplyRhs = .FALSE.
       multname = ListGetString( Solver % Values,'Matrix Multiplier Name',Found )
       IF(Found ) THEN
         cvar => VariableGet( Solver % Mesh % Variables, multname, UnfoundFatal = .TRUE.)
         IF(.NOT. ASSOCIATED(cVar) ) THEN
           CALL Fatal('UseLocalMatrixStorage','Field not found: '//TRIM(multname))
         END IF
         IF( cvar % TYPE /= Variable_on_elements ) THEN
           CALL Fatal('UseLocalMatrixStorage','"Field should be elemental: '//TRIM(multname))
         END IF
         DoMultiply = .TRUE.

         ! We may multiply the r.h.s. with a different multiplier. 
         multname = ListGetString( Solver % Values,'Rhs Multiplier Name',Found )
         IF(Found ) THEN
           rvar => VariableGet( Solver % Mesh % Variables, multname, UnfoundFatal = .TRUE.)
           IF(.NOT. ASSOCIATED(rVar) ) THEN
             CALL Fatal('UseLocalMatrixStorage','Field not found: '//TRIM(multname))
           END IF
           IF( rvar % TYPE /= Variable_on_elements ) THEN
             CALL Fatal('UseLocalMatrixStorage','"Field should be elemental: '//TRIM(multname))
           END IF
           DoMultiplyRhs = .TRUE.
         END IF
       END IF
       prevSolverId = Solver % SolverId       
     END IF

     ! Multiply locally stored matrix. Possible use is, for example, density in topology optimization. 
     IF( DoMultiply ) THEN
       vind = cvar % Perm(eind)
       IF(vind > 0 ) THEN
         cmult = cvar % Values(vind)       
         K(1:n,1:n) = cmult * K(1:n,1:n)
       END IF
       IF( DoMultiplyRhs ) THEN
         vind = rvar % Perm(eind)
         IF( vind > 0 ) THEN
           rmult = rvar % Values(vind)
           F(1:n) = rmult * F(1:n)
         END IF
       END IF         
     END IF
           
   END SUBROUTINE UseLocalMatrixStorage

   
!---------------------------------------------------------------------------
!> Obtain local matrix, e.g. for topology optimization.
!> If the elements are alike, the element index may point to a different
!> element than itself. 
!---------------------------------------------------------------------------
   SUBROUTINE GetLocalMatrixStorage( Solver, n, K, F, Found, elemind, activeind ) 
     TYPE(Solver_t) :: Solver
     INTEGER :: n
     REAL(KIND=dp) :: K(:,:), F(:)
     LOGICAL :: Found
     INTEGER, OPTIONAL :: elemind, activeind
  
     TYPE(LocalSystemStorage_t), POINTER :: pLocal
     INTEGER :: eind

     Found = .FALSE.
     IF( PRESENT(activeind) ) THEN
       eind = activeind
     ELSE IF( PRESENT( elemind ) ) THEN
       IF( elemind > SIZE(Solver % InvActiveElements ) ) RETURN
       eind = Solver % InvActiveElements(elemind)
       IF(eind==0) RETURN
     ELSE
       CALL Fatal('GetLocalMatrixStorage','Give element index as a parameter in some way!')
     END IF
     
     pLocal => Solver % LocalSystem(eind) 
     IF(eind /= pLocal % eind ) THEN
       pLocal => Solver % LocalSystem(pLocal % eind)
     END IF
     IF(pLocal % eind < 1 ) RETURN
     
     IF(pLocal % n == n ) THEN
       K(1:n,1:n) = pLocal % K(1:n,1:n)
       F(1:n) = pLocal % F(1:n)
       Found = .TRUE.
     END IF
       
   END SUBROUTINE GetLocalMatrixStorage

   
   
   
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
    INTEGER :: N               !< The count of nodal, edge and face degrees of freedom.
    INTEGER :: Nb              !< The count of internal (bubble) degrees of freedom.
    REAL(KIND=dp) :: K(:,:)    !< Local stiffness matrix.
    REAL(KIND=dp), OPTIONAL :: F(:)   !< Local force vector.
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

    IF (PRESENT(F)) F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
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
    INTEGER :: N               !< The count of nodal, edge and face degrees of freedom.
    INTEGER :: Nb              !< The count of internal (bubble) degrees of freedom.
    COMPLEX(KIND=dp) :: K(:,:)    !< Local stiffness matrix.
    COMPLEX(KIND=dp), OPTIONAL :: F(:)   !< Local force vector.
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

    IF (PRESENT(F)) F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    IF (PRESENT(F1)) THEN
      Fb  = F1(Bdofs)
      F1(1:n) = F1(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    END IF

    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE CondensatePC
!------------------------------------------------------------------------------


END MODULE MatrixAssembly
