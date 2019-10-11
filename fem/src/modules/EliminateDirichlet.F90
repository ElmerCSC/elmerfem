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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 13 Sep 2002
! *
! ****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!------------------------------------------------------------------------------
!>  Subroutine for reducing a linsystem DOFs by eliminating DOFs corresponding
!>  to known DOF values.
!------------------------------------------------------------------------------
 RECURSIVE INTEGER FUNCTION EliminateDirichlet( Model, Solver, A, b, x, n, DOFs, Norm )
!------------------------------------------------------------------------------
  USE Types
  USE Lists
  USE SolverUtils
  USE CRSmatrix
  USE GeneralUtils

  IMPLICIT NONE
  
  TYPE(model_t)  :: Model        !> All model information (mesh,materials,BCs,etc...)
  TYPE(solver_t) :: Solver       !> Linear equation solver options
  TYPE(matrix_t), POINTER :: A   !> Linear equation matrix information
  INTEGER :: DOFs                !> Number of degrees of freedom of the equation
  INTEGER :: n                   !> Length of unknown vector
  REAL(KIND=dp) :: b(n)          !> The unknown in the linear equation
  REAL(KIND=dp) :: x(n)          !> The right hand side of the linear equation
  REAL(KIND=dp) :: Norm          !> The norm to determine the convergence of the nonlinear system
!------------------------------------------------------------------------------

  TYPE(Matrix_t), POINTER :: Bmatrix, Cmatrix, SaveMatrix
  REAL(KIND=dp), POINTER :: F(:), U(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) ::  TotTime, LinTime, at,s,r, br, bi
#else
  REAL(KIND=dp) ::  CPUtime, TotTime, LinTime, at,s,r, br, bi
#endif

  LOGICAL :: stat, GotIt, EigAnal, DoChange

  INTEGER :: istat, ii, i, j, k, l, m, t, Dirichlet, Level = 0, DOF, zz
  INTEGER, POINTER :: Permutation(:), DCount(:), PermSave(:)
  LOGICAL :: AllocationsDone = .FALSE., Eigen = .FALSE., &
             Mass = .FALSE., Damp=.FALSE., Add

  COMPLEX(KIND=dp), ALLOCATABLE :: TempEigenVectors(:,:)

  SAVE Bmatrix, AllocationsDone

!------------------------------------------------------------------------------

   TotTime = CPUtime()

!------------------------------------------------------------------------------
   BMatrix => A % EMatrix

   Mass = .FALSE.
   IF ( ASSOCIATED( A % MassValues ) ) THEN
      IF ( SIZE( A % MassValues ) == SIZE( A % Values ) ) Mass = .TRUE.
   END IF

   Damp = .FALSE.
   IF ( ASSOCIATED( A % DampValues ) ) THEN
      IF ( SIZE( A % DampValues ) == SIZE( A % Values ) ) Damp = .TRUE.
   END IF

   IF ( ASSOCIATED( BMatrix ) ) THEN
      IF ( Mass .AND. .NOT. ASSOCIATED(BMatrix % MassValues) ) THEN
         ALLOCATE( Bmatrix % MassValues( SIZE(BMatrix % Values) ) )
      END IF
   END IF

   IF ( ASSOCIATED( BMatrix ) ) THEN
      IF ( Damp .AND. .NOT. ASSOCIATED(BMatrix % DampValues) ) THEN
         ALLOCATE( Bmatrix % DampValues( SIZE(BMatrix % Values) ) )
      END IF
   END IF

   IF ( .NOT. ASSOCIATED( BMatrix ) ) THEN
      IF ( ASSOCIATED( A % Perm ) ) THEN
         Permutation => A % Perm
      ELSE
        ALLOCATE( Permutation( A % NumberOfRows ) )
      END IF

      Permutation = 0
      Dirichlet = 0

      l = 0
      DO i=1,A % NumberOFRows
         s = A % Values( A % Diag(i) )
         A % Values( A % Diag(i) ) = 0.0d0 
         j = A % Rows(i)
         k = A % Rows(i+1)-1
         IF ( ALL( ABS(A % Values(j:k))<1.d-12 ) ) THEN
            Dirichlet = Dirichlet + 1
         ELSE
            l = l + 1
            Permutation(i) = l
         END IF
         A % Values( A % Diag(i) ) = s
      END DO

!     IF ( Dirichlet <= 0 ) THEN
!         DEALLOCATE( Permutation )
!         EliminateDirichlet = 0
!         RETURN
!     END IF

      Bmatrix => AllocateMatrix()

      Bmatrix % NumberOFRows = A % NumberOFRows - Dirichlet
      ALLOCATE( Bmatrix % Rows( Bmatrix % NumberOfRows + 1 ), &
                Bmatrix % Diag( Bmatrix % NumberOfRows ) )

      j = 0
      k = 1
zz = 0
      DO i=1, A % NumberOFRows 
         IF ( Permutation(i) /= 0 ) THEN
            j = j + 1
            Bmatrix % Rows(j) = k
            DO l = A % Rows(i),A % Rows(i+1)-1
               t = A % Cols(l)
               Add = .TRUE.
               IF ( .NOT. A % Complex ) THEN
                 Add = A % Values(l) /= 0
                 IF ( Mass ) Add = Add .OR. A % MassValues(l) /= 0
                 IF ( Damp ) Add = Add .OR. A % DampValues(l) /= 0
if ( .not. add ) zz=zz+1
               END IF
               Add = Add .OR. l==A % Diag(i)
               Add = Add .AND. Permutation(t) /= 0
               IF ( Add ) k = k + 1
            END DO
         END IF
      END DO
!print*,'removing zeros: ', zz
      Bmatrix % Rows(Bmatrix % NumberOFRows+1) = k

      IF ( Mass ) THEN
         ALLOCATE( Bmatrix % MassValues(k-1) )
      END IF
      IF ( Damp ) THEN
         ALLOCATE( Bmatrix % DampValues(k-1) )
      END IF
      ALLOCATE( Bmatrix % Values(k-1), Bmatrix % Cols(k-1), &
          BMatrix % RHS( Bmatrix % NumberOfRows ) )

      ALLOCATE( Dcount( A % NumberOfRows ) )
      IF ( Permutation(1) == 0 ) THEN
         DCount(1) = 1
      ELSE
         DCount(1) = 0
      END IF

      DO i=2, A % NumberOFRows 
         IF ( Permutation(i) /= 0 ) THEN
            DCount(i) = DCount(i-1)
         ELSE
            DCount(i) = DCount(i-1) + 1
         END IF
      END DO

      j = 0
      k = 1
      DO i=1, A % NumberOFRows 
        IF ( Permutation(i) /= 0 ) THEN
           j = j + 1
           DO l = A % Rows(i),A % Rows(i+1)-1
              t = A % Cols(l)
              Add = .TRUE.
              IF ( .NOT. A % Complex ) THEN
                Add = A % Values(l) /= 0
                IF ( Mass ) Add = Add .OR. A % MassValues(l) /= 0
                IF ( Damp ) Add = Add .OR. A % DampValues(l) /= 0
              END IF
              Add = Add .OR. l==A % Diag(i)
              Add = Add .AND. Permutation(t) /= 0
              IF ( Add ) THEN
                 Bmatrix % Cols(k) = t - DCount(t)
                 k = k + 1
              END IF
           END DO
        END IF
      END DO

     DEALLOCATE( DCount )

      DO i=1,Bmatrix % NumberOfRows
         DO j=Bmatrix % Rows(i), Bmatrix % Rows(i+1)-1
            IF ( Bmatrix % Cols(j) == i ) THEN
               Bmatrix % Diag(i) = j
            END IF
         END DO
      END DO
      CALL CRS_SortMatrix( BMatrix )

      A % EMatrix => BMatrix
      A % Perm => Permutation
   END IF

   Permutation => A % Perm
   F => BMatrix % RHS
   ALLOCATE( U( Bmatrix % NumberOfRows ) )

   j  = 1
   k  = 1
   IF ( Solver % Matrix % Complex ) THEN
      DO i=1, A % NumberOFRows,2
        IF ( Permutation(i) /= 0 ) THEN
           F(i)   = B(j)
           F(i+1) = B(j+1)
           U(i)   = X(j)
           U(i+1) = X(j+1)

           k  = BMatrix % Rows(j)
           DO l = A % Rows(i), A % Rows(i+1)-1,2
              m = A % Cols(l)

              IF ( Permutation(m) == 0 ) THEN
                br = b(m)   / A  % Values(A % Diag(m))
                bi = b(m+1) / A  % Values(A % Diag(m+1))
                F(j)=F(j)-(A % Values(l)*br+A % Values(l+1)*bi)
                F(j+1)=F(j+1)-(A % Values(l)*bi-A % Values(l+1)*br)
              ELSE
                 Bmatrix % Values(k)   = A % Values(l)
                 Bmatrix % Values(k+1) = A % Values(l+1)
                 IF ( Mass ) THEN
                    Bmatrix % MassValues(k)   = A % MassValues(l)
                    Bmatrix % MassValues(k+1) = A % MassValues(l+1)
                 END IF
                 IF ( Damp ) THEN
                    Bmatrix % DampValues(k)   = A % DampValues(l)
                    Bmatrix % DampValues(k+1) = A % DampValues(l+1)
                 END IF
                 k = k + 2
              END IF
           END DO
           j = j + 2
        END IF
      END DO

      DO i=1,BMatrix % NumberOfRows,2
         k = BMatrix % Rows(i+1)
         DO j=BMatrix % Rows(i),BMatrix % Rows(i+1)-1,2
            BMatrix % Values(k)   = -BMatrix % Values(i+1)
            BMatrix % Values(k+1) =  BMatrix % Values(i)
            IF ( Mass ) THEN
               Bmatrix % MassValues(k)   = -BMatrix % MassValues(i+1)
               Bmatrix % MassValues(k+1) =  BMatrix % MassValues(i)
            END IF
            IF ( Damp ) THEN
               Bmatrix % DampValues(k)   = -BMatrix % DampValues(i+1)
               Bmatrix % DampValues(k+1) =  BMatrix % DampValues(i)
            END IF
            k = k + 2
         END DO
      END DO
   ELSE
      j = 0
      DO i=1, A % NumberOFRows 
        IF ( Permutation(i) /= 0 ) THEN
           j = j + 1
           F(j) = B(i) 
           U(j) = X(i)
           DO l = A % Rows(i), A % Rows(i+1)-1
              t = A % Cols(l)
              IF ( Permutation(t) == 0 ) THEN
                 F(j) = F(j) - A % Values(l) * b(t) / A % Values( A % Diag(t) )
              ELSE
                 Add = A % Values(l) /= 0
                 IF ( Mass ) Add = Add .OR. A % MassValues(l) /= 0
                 IF ( Damp ) Add = Add .OR. A % DampValues(l) /= 0
                 Add = Add .OR. l==A % Diag(i)
                 IF ( Add ) THEN
                   Bmatrix % Values(k) = A % Values(l)
                   IF ( Mass ) THEN
                      Bmatrix % MassValues(k) = A % MassValues(l)
                   END IF
                   IF ( Damp ) THEN
                      Bmatrix % DampValues(k) = A % DampValues(l)
                   END IF
                   k = k + 1
                 END IF
              END IF
           END DO
        END IF
      END DO
   END IF

  WRITE( Message, * ) 'Eliminated ', TRIM(I2S(A % NumberOfRows - &
        Bmatrix % NumberOFrows)), ' unknowns out of ', TRIM(i2S(A % NumberOfRows))
  CALL Info('Dirichlet',Message,Level=5)
!------------------------------------------------------------------------------
!   Solve the system
!------------------------------------------------------------------------------
  Bmatrix % Lumped    = A % Lumped
  Bmatrix % Complex   = A % Complex
  Bmatrix % Symmetric = A % Symmetric

!print*,k, size(a % values), a % numberofrows, bmatrix % numberofrows

  IF ( ParEnv % PEs > 1 ) THEN
     PermSave => Solver % Variable % Perm
     CMatrix => Solver % Matrix
     Solver % Matrix => BMatrix

     ALLOCATE( Solver % Variable % Perm( SIZE(PermSave) ) )
     Solver % Variable % Perm = 0
     DO i=1,SIZE( PermSave )
       j = PermSave(i)
       IF ( j > 0 ) THEN
          DOF = DOFs * (j-1) + 1
          IF ( Permutation(DOF) > 0 ) THEN
              Solver % Variable % Perm(i) = ( Permutation(DOF)-1 ) / DOFs + 1
          END IF
       END IF
     END DO

     IF ( .NOT. ASSOCIATED( BMatrix % ParMatrix ) ) THEN
        CALL ParallelInitMatrix( Solver, BMatrix )
      END IF
  END IF


  at = CPUtime()

! Eliminate calling of ComputeChange since the it is better to make it for the original system matrix
!----------------------------------------------------------------------------------------------
  DoChange = ListGetLogical( Solver % Values,'Skip Compute Nonlinear Change',GotIt)
  CALL ListAddLogical( Solver % Values,'Skip Compute Nonlinear Change',.TRUE.)

  CALL SolveLinearSystem( Bmatrix, F, U, Norm, DOFs, Solver )

  CALL ListAddLogical( Solver % Values,'Skip Compute Nonlinear Change',DoChange)

  LinTime = CPUtime() - at
  WRITE( Message, * ) 'Time used in SolveLinearSystem (CPU): ', LinTime
  CALL Info( 'Dirichlet', Message, Level=5 )

  IF ( ParEnv % PEs > 1 ) THEN
     Solver % Matrix => CMatrix
     DEALLOCATE( Solver % Variable % Perm )
     Solver % Variable % Perm => PermSave
  END IF

!------------------------------------------------------------------------------
!   Map the result back to original nodes
!------------------------------------------------------------------------------

  IF ( Solver % NOFEigenValues > 0 ) THEN
     ALLOCATE( TempEigenVectors( Solver % NOFEigenValues, A % NumberOfRows ) )
     TempEigenVectors = Solver % Variable % EigenVectors

     j = 0
     IF ( Solver % Matrix % Complex ) THEN
        DO i=1,A % NumberOFRows / 2
           IF ( Permutation(i) == 0 ) THEN
              DO k=1,Solver % NOFEigenValues
                 Solver % Variable % EigenVectors(k,i) = CMPLX( &
                     B(2*i) / A % Values( A % Diag(2*i) ), B(2*i+1) / A % Values( A % Diag(2*i+1) ),KIND=dp )
              END DO
           ELSE
              j = j + 1
              DO k=1,Solver % NOFEigenValues
                 Solver % Variable % EigenVectors(k,i) = TempEigenVectors(k,j)
              END DO
           END IF
        END DO
     ELSE
        DO i=1,A % NumberOFRows
           IF ( Permutation(i) == 0 ) THEN
              DO k=1,Solver % NOFEigenValues
                 Solver % Variable % EigenVectors(k,i) = B(i) / A % Values( A % Diag(i) )
              END DO
           ELSE
              j = j + 1
              DO k=1,Solver % NOFEigenValues
                 Solver % Variable % EigenVectors(k,i) = TempEigenVectors(k,j)
              END DO
           END IF
        END DO
     END IF
  ELSE
     j = 0
     DO i=1,A % NumberOFRows
        IF ( Permutation(i) == 0 ) THEN
           X(i) = B(i) / A % Values( A % Diag(i) )
        ELSE
           j = j + 1
           X(i) = U(j)
        END IF
     END DO
  END IF

  IF ( Solver % NOFEigenValues > 0 ) DEALLOCATE( TempEigenVectors )
!------------------------------------------------------------------------------

! CALL FreeMatrix( Bmatrix )
  DEALLOCATE( U )

  at = CPUTime()

  CALL ComputeChange(Solver, .FALSE.)
  LinTime = LinTime + (CPUTime() - at)
  ! TotTime = CPUTime(Solver,.FALSE.) - TotTime
  TotTime = CPUTime() - TotTime

  WRITE( Message, * ) 'Total time spent in DirichletReduction (CPU): ', TotTime 
  CALL Info( 'Dirichlet', Message, Level=5 )
  WRITE( Message, * ) 'Additional time spent in DirichletReduction (CPU): ', TotTime - LinTime
  CALL Info( 'Dirichlet', Message, Level=5 )
  CALL Info( 'Dirichlet', ' ', Level=5 )

  CALL FreeMatrix( A % EMatrix )
  NULLIFY( A % EMatrix )
  EliminateDirichlet = 1

END FUNCTION EliminateDirichlet

!> \} ElmerLib
