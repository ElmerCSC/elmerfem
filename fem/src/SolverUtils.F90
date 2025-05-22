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

!> Basic utilities used by individual solvers. 
!------------------------------------------------------------------------------

!> \ingroup ElmerLib
!> \{


MODULE SolverUtils

#include "../config.h"

   USE LoadMod
   USE Multigrid
   USE ElementUtils
   USE LumpingUtils
   USE TimeIntegrate
   USE ModelDescription
   USE MeshUtils, ONLY : CommunicateParallelSystemTag, CylinderFit, &
       DisplaceMesh, FindExtremumNodes, FindMeshEdges, GetLagrangeIndexes, &
       IntegralProjector, PeriodicProjector, SaveProjector, MakePermUsingMask
   USE ParallelUtils
   USE ParallelEigenSolve
   USE MatrixAssembly
   
   IMPLICIT NONE

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
     CHARACTER(:), ALLOCATABLE :: str
     TYPE(NormalTangential_t), POINTER :: NT
!------------------------------------------------------------------------------
     
     CALL Info('InitializeToZero','Initializing the linear system to zero',Level=12)
     
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

     IF ( Solver % Variable % DOFs <= 1 ) RETURN

     str = 'Normal-Tangential'
     IF ( SEQL(Solver % Variable % Name, 'flow solution') ) THEN
       str = TRIM(str) // ' Velocity'
     ELSE
       str = TRIM(str) // ' ' // GetVarName(Solver % Variable)
     END IF
     AnyNT = ListGetLogicalAnyBC( CurrentModel, str )

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
       NT => CurrentModel % Solver % NormalTangential
       NT % NormalTangentialNOFNodes = 0
       NT % NormalTangentialName = TRIM(str)

       CALL CheckNormalTangentialBoundary( CurrentModel, NT % NormalTangentialName, &
           NT % NormalTangentialNOFNodes, NT % BoundaryReorder, &
           NT % BoundaryNormals, NT % BoundaryTangent1, NT % BoundaryTangent2, dim )
       
       CALL AverageBoundaryNormals( CurrentModel, NT % NormalTangentialName, &
           NT % NormalTangentialNOFNodes, NT % BoundaryReorder, &
           NT % BoundaryNormals, NT % BoundaryTangent1, NT % BoundaryTangent2, &
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
   SUBROUTINE MaskedMatrixVectorMultiply( A,u,v,ActiveRow,ActiveCol )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n
     REAL(KIND=dp), DIMENSION(:) CONTIG :: u,v
     LOGICAL, DIMENSION(:) :: ActiveRow
     LOGICAL, DIMENSION(:) :: ActiveCol
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_MaskedMatrixVectorMultiply( A,u,v,ActiveRow, ActiveCol )

     CASE DEFAULT
       CALL Fatal('MaskedMatrixVectorMultiply','Not implemented for List matrix type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE MaskedMatrixVectorMultiply
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
     REAL(KIND=dp) :: s, t, zeta
     REAL(KIND=dp) :: PrevSol(DOFs*n,Solver % Order), CurSol(DOFs*n), LForce(n*DOFs)
     TYPE(Variable_t), POINTER :: DtVar
     REAL(KIND=dp) :: Dts(Solver % Order)
     LOGICAL :: ConstantDt
     TYPE(Element_t), POINTER :: Element
     CHARACTER(:), ALLOCATABLE :: Method
!------------------------------------------------------------------------------
     INTEGER :: PredCorrOrder       !< Order of predictor-corrector scheme

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
       
     CASE('adams-bashforth')
       zeta = ListGetConstReal( Solver % Values, 'Adams Zeta', GotIt )
       IF ( .NOT. Gotit) zeta = 1.0_dp
       PredCorrOrder = ListGetInteger( Solver % Values, &
           'Predictor-Corrector Scheme Order', GotIt)
       IF (.NOT. GotIt) PredCorrOrder = 2
       PredCorrOrder = MIN(PredCorrOrder, Solver % DoneTime /2)       
       CALL AdamsBashforth( n*DOFs, dt, MassMatrix, StiffMatrix, Force, &
           PrevSol(:,1), zeta, PredCorrOrder)
       
     CASE('adams-moulton')
       PredCorrOrder = ListGetInteger( Solver % Values, &
           'Predictor-Corrector Scheme Order', GotIt)
       IF (.NOT. GotIt) PredCorrOrder = 2
       PredCorrOrder = MIN(PredCorrOrder, Solver % DoneTime /2)
       CALL AdamsMoulton( n*DOFs, dt, MassMatrix, StiffMatrix, Force, &
           PrevSol, PredCorrOrder )      
       
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
     TYPE(Variable_t), POINTER :: DtVar
     REAL(KIND=dp) :: Dts(Solver % Order)
     REAL(KIND=dp), POINTER :: PrevSol(:,:), ML(:), CurrSol(:)
     INTEGER, POINTER :: Rows(:), Cols(:)
     LOGICAL :: ConstantDt, Lumped, Found
     CHARACTER(:), ALLOCATABLE :: Method
!------------------------------------------------------------------------------

     CALL Info('Add1stOrderTime_CRS','Adding time discretization to CRS matrix',Level=20)

!------------------------------------------------------------------------------
     Order = MIN(Solver % DoneTime, Solver % Order)
     IF(Order == 0 ) THEN
       CALL Info('Add1stOrderTime_CRS','Zeroth order, nothing to do!',Level=20)
       RETURN
     END IF
     
     Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )
     CurrSol => Solver % Variable % Values
     PrevSol => Solver % Variable % PrevValues
     
     SELECT CASE( Method )
       
     CASE( 'fs' ) 
       CALL FractionalStep_CRS( dt, Matrix, Force, PrevSol(:,1), Solver )

     CASE('bdf')
       ConstantDt = .TRUE.
       IF(Order > 1) THEN
         Dts(1) = Dt
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

     CASE('adams-bashforth')
       CALL Fatal('Add1stOrderTime_CRS','Not implemented for method: '//TRIM(Method))
       
     CASE('adams-moulton')
       CALL Fatal('Add1stOrderTime_CRS','Not implemented for method: '//TRIM(Method))
       
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
     REAL(KIND=dp) :: s,t
!    CHARACTER(:), ALLOCATABLE :: Method
     REAL(KIND=dp) :: X(DOFs*n),V(DOFs*N),A(DOFs*N),A2(DOFs*N), LForce(n*DOFs)

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
!          SELECT CASE(Method)
!          CASE DEFAULT
             X(K) = Solver % Variable % PrevValues(L,3)
             V(K) = Solver % Variable % PrevValues(L,4)
             A(K) = Solver % Variable % PrevValues(L,5)
             A2(K) = Solver % Variable % PrevValues(L,7)
!          END SELECT
         END IF
       END DO
     END DO

     LForce(1:n*DOFs) = Force(1:n*DOFs)
     CALL UpdateGlobalForce( Solver % Matrix % Force(:,1), LForce, &
                  n, DOFs, NodeIndexes )
!------------------------------------------------------------------------------
!    Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )
!    SELECT CASE(Method)
     CALL Time2ndOrder( n*DOFs, dt, MassMatrix, DampMatrix, StiffMatrix, &
                Force, X, V, A, A2, Solver % Alpha, Solver % Beta )
!    END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE Add2ndOrderTime
!------------------------------------------------------------------------------

   
   SUBROUTINE Add2ndOrderTime_CRS( Matrix, Force, &
       dt, PrevValues, Solver )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: Matrix  !< Global matrix (including stiffness and mass)
     REAL(KIND=dp) :: Force(:)          !< Global right-hand-side vector.
     REAL(KIND=dp) :: dt                !< Simulation timestep size
     REAL(KIND=dp), POINTER :: PrevValues(:,:)
     TYPE(Solver_t), POINTER :: Solver
     
     IF( Matrix % Lumped ) THEN
       CALL Fatal('Add2ndOrderTime_CRS','Implement matrix lumping for this!')
     END IF

     CALL Time2ndOrder_CRS(dt, Matrix, Force, Solver % Alpha, Solver % Beta, PrevValues )
     
   END SUBROUTINE Add2ndOrderTime_CRS
     

   
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
      ForceVector, LocalForce, n, NDOFs, DofIndexes, RotateNT, UElement, &
              GlobalValues )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix  !< The global matrix
     REAL(KIND=dp) :: LocalStiffMatrix(:,:)  !< Local matrix to be added to the global matrix.
     REAL(KIND=dp) :: LocalForce(:)          !< Element local force vector.
     REAL(KIND=dp) :: ForceVector(:)         !< The global RHS vector.
     INTEGER :: n                            !< Number of degrees of freedom in element for each component
     INTEGER :: NDOFs                        !< Number of components for vector field.
     INTEGER :: DofIndexes(:)                !< Element node/edge/face to global node/edge/face numbering mapping.
     LOGICAL, OPTIONAL :: RotateNT           !< Should the global equation be done in local normal-tangential coordinates.
     TYPE(Element_t), OPTIONAL, TARGET :: UElement !< Element to be updated
     REAL(KIND=dp), OPTIONAL :: GlobalValues(:)
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,np,dim, NormalIndexes(n), pIndexes(64)
     LOGICAL :: Rotate
     TYPE(Element_t), POINTER :: Element
     TYPE(NormalTangential_t), POINTER :: NT
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
     IF( ndofs < dim ) Rotate = .FALSE.

     IF( Rotate ) THEN
       NT => CurrentModel % Solver % NormalTangential
       Rotate = ( NT % NormalTangentialNOFNodes > 0 )
     END IF

     IF ( Rotate ) THEN
       NormalIndexes = 0

       np = mGetElementDOFs(pIndexes,Element)
       np = MIN(np, n)
       DO i=1,np
         j = pIndexes(i)
         IF(j>0 .AND. j<= SIZE(NT % BoundaryReorder)) THEN
           NormalIndexes(i) = NT % BoundaryReorder(j)
         END IF
       END DO
       
       CALL RotateMatrix( LocalStiffMatrix, LocalForce, n, dim, NDOFs, &
          NormalIndexes, NT % BoundaryNormals, NT % BoundaryTangent1, NT % BoundaryTangent2 )
     END IF
!------------------------------------------------------------------------------
     IF ( ASSOCIATED( StiffMatrix ) ) THEN
       SELECT CASE( StiffMatrix % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_GlueLocalMatrix( StiffMatrix,n,NDOFs, &
                      DofIndexes, LocalStiffMatrix, GlobalValues )

       CASE( MATRIX_LIST )
         CALL List_GlueLocalMatrix( StiffMatrix % ListMatrix,n,NDOFs,DofIndexes, &
                          LocalStiffMatrix )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_GlueLocalMatrix( StiffMatrix,n,NDOFs,DofIndexes, &
                          LocalStiffMatrix )
       END SELECT
     END IF

     DO i=1,n
       IF ( DofIndexes(i) > 0 ) THEN
         DO j=1,NDOFs
           k = NDOFs * (DofIndexes(i)-1) + j
!$omp atomic
           ForceVector(k) = ForceVector(k) + LocalForce(NDOFs*(i-1)+j)
         END DO
       END IF
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE UpdateGlobalEquations
!------------------------------------------------------------------------------


!> Add element local matrices & vectors to global matrices and vectors.
!> Vectorized version, does not support normal or tangential boundary
!> conditions yet.
   SUBROUTINE UpdateGlobalEquationsVec( Gmtr, Lmtr, Gvec, Lvec, n, &
           NDOFs, DofIndexes, RotateNT, UElement, MCAssembly )
     TYPE(Matrix_t), POINTER :: Gmtr         !< The global matrix
     REAL(KIND=dp) CONTIG :: Lmtr(:,:)              !< Local matrix to be added to the global matrix.
     REAL(KIND=dp) CONTIG :: Gvec(:)                !< Element local force vector.
     REAL(KIND=dp) CONTIG :: Lvec(:)                !< The global RHS vector.
     INTEGER :: n                            !< Number of nodes.
     INTEGER :: NDOFs                        !< Number of degrees of free per node.
     INTEGER CONTIG :: DofIndexes(:)               !< Element node to global node numbering mapping.
     LOGICAL, OPTIONAL :: RotateNT           !< Should the global equation be done in local normal-tangential coordinates.
     TYPE(Element_t), OPTIONAL, TARGET :: UElement !< Element to be updated
     LOGICAL, OPTIONAL :: MCAssembly   !< Assembly process is multicoloured and guaranteed race condition free 

     ! Local variables
     INTEGER :: dim, i,j,k,np
     INTEGER :: NormalIndexes(n),pIndexes(64)
     REAL(KIND=dp) :: Vals(n*NDOFs)
!DIR$ ATTRIBUTES ALIGN:64::Vals

     TYPE(Element_t), POINTER :: Element
     LOGICAL :: Rotate
     LOGICAL :: ColouredAssembly, NeedMasking
     TYPE(NormalTangential_t), POINTER :: NT
     
     IF (PRESENT(UElement)) THEN
       Element => UElement
     ELSE
       Element => CurrentModel % CurrentElement
     END IF
     
     IF ( CheckPassiveElement(Element) )  RETURN
     Rotate = .TRUE.
     IF ( PRESENT(RotateNT) ) Rotate = RotateNT
     
     ColouredAssembly = .FALSE.
     IF ( PRESENT(MCAssembly) ) ColouredAssembly = MCAssembly

     dim = CoordinateSystemDimension()

     IF(ndofs < dim) Rotate = .FALSE.

     IF( Rotate ) THEN
       NT => CurrentModel % Solver % NormalTangential
       Rotate = ( NT % NormalTangentialNOFNodes > 0 )
     END IF

     IF ( Rotate ) THEN
       NormalIndexes = 0

       np = mGetElementDOFs(pIndexes,Element)
       np = MIN(np,n)

       DO i=1,np
         j = pIndexes(i)
         IF(j>0 .AND. j<= SIZE(NT % BoundaryReorder)) THEN
           NormalIndexes(i) = NT % BoundaryReorder(j)
         END IF
       END DO
       
       ! TODO: See that RotateMatrix is vectorized
       CALL RotateMatrix( Lmtr, Lvec, n, dim, NDOFs, NormalIndexes, NT % BoundaryNormals, &
                    NT % BoundaryTangent1, NT % BoundaryTangent2 )

       !IF ( Rotate .AND. NormalTangentialNOFNodes > 0 .AND. ndofs>=dim) THEN
       !  CALL Fatal('UpdateGlobalEquationsVec', &
       !          'Normal or tangential boundary conditions not supported yet!')
     END IF

     NeedMasking = .FALSE.
     DO i=1,n
       IF (DofIndexes(i)<=0) THEN
         NeedMasking = .TRUE.
         EXIT
       END IF
     END DO
     
     IF ( ASSOCIATED( Gmtr ) ) THEN
       SELECT CASE( Gmtr % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_GlueLocalMatrixVec(Gmtr, n, NDOFs, DofIndexes, Lmtr, ColouredAssembly, NeedMasking)
       CASE DEFAULT
         CALL Fatal('UpdateGlobalEquationsVec','Not implemented for given matrix type')
       END SELECT
     END IF
     
     ! Check for multicolored assembly
     IF (ColouredAssembly) THEN
       IF (NeedMasking) THEN
         ! Vector masking needed, no ATOMIC needed
         !_ELMER_OMP_SIMD PRIVATE(j,k)
         DO i=1,n
           IF (DofIndexes(i)>0) THEN
             DO j=1,NDOFs
               k = NDOFs*(DofIndexes(i)-1) + j
               Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
             END DO
           END IF
         END DO
       ELSE
         ! No vector masking needed, no ATOMIC needed
         IF (NDOFS>1) THEN
           !_ELMER_OMP_SIMD PRIVATE(j,k)
           DO i=1,n
             DO j=1,NDOFs
               k = NDOFs*(DofIndexes(i)-1) + j
               Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
             END DO
           END DO
         ELSE
           !_ELMER_OMP_SIMD
           DO i=1,n
             Gvec(DofIndexes(i)) = Gvec(DofIndexes(i)) + Lvec(i)
           END DO
         END IF
       END IF ! Vector masking
     ELSE
       IF (NeedMasking) THEN
         ! Vector masking needed, ATOMIC needed
         DO i=1,n
           IF (DofIndexes(i)>0) THEN
!DIR$ IVDEP
             DO j=1,NDOFs
               k = NDOFs*(DofIndexes(i)-1) + j
               !$OMP ATOMIC
               Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
             END DO
           END IF
         END DO
       ELSE
         ! No vector masking needed, ATOMIC needed
         DO i=1,n
!DIR$ IVDEP
           DO j=1,NDOFs
             k = NDOFs*(DofIndexes(i)-1) + j
             !$OMP ATOMIC
             Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
           END DO
         END DO
       END IF ! Vector masking
     END IF ! Coloured assembly
   END SUBROUTINE UpdateGlobalEquationsVec

!------------------------------------------------------------------------------
!> Update the global vector with the local vector entry.
!------------------------------------------------------------------------------
   SUBROUTINE UpdateGlobalForce(ForceVector, LocalForce, n, &
             NDOFs, DofIndexes, RotateNT, UElement )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: LocalForce(:)          !< Element local force vector.
     REAL(KIND=dp) :: ForceVector(:)         !< The global RHS vector.
     INTEGER :: n                            !< Number of nodes.
     INTEGER :: NDOFs                        !< Number of element nodes. 
     INTEGER :: DofIndexes(:)                !< Element node/edge/face to global node/edge/face numbering mapping.
     LOGICAL, OPTIONAL :: RotateNT           !< Should the global equation be done in local normal-tangential coordinates.
     TYPE(Element_t), OPTIONAL, TARGET :: UElement !< Element to be updated
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k,np,dim,NormalIndexes(n),pIndexes(64)
     LOGICAL :: Rotate
     REAL(KIND=dp) :: LocalStiffMatrix(n*NDOFs,n*NDOFs), LForce(n*NDOFs)
     TYPE(NormalTangential_t), POINTER :: NT
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
     IF( Rotate ) THEN
       NT => CurrentModel % Solver % NormalTangential
       Rotate = ( NT % NormalTangentialNOFNodes > 0 )
     END IF
            
     IF ( Rotate ) THEN
       dim = CoordinateSystemDimension()

       NormalIndexes = 0

       np = mGetElementDOFs(pIndexes,Element)
       np = MIN(np,n)

       DO i=1,np
         j = pIndexes(i)
         IF(j>0 .AND. j<= SIZE(NT % BoundaryReorder)) THEN
           NormalIndexes(i) = NT % BoundaryReorder(j)
         END IF
       END DO

       CALL RotateMatrix( LocalStiffMatrix, LocalForce, n, dim, NDOFs, &
          NormalIndexes, NT % BoundaryNormals, NT % BoundaryTangent1, NT % BoundaryTangent2 )
     END IF

     DO i=1,n
       IF ( DofIndexes(i) > 0 ) THEN
         DO j=1,NDOFs
           k = NDOFs * (DofIndexes(i)-1) + j
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
              n, NDOFs, DofIndexes, GlobalValues )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix  !< The global matrix structure
     REAL(KIND=dp) :: LocalMassMatrix(:,:)   !< Local matrix to be added to the global matrix
     INTEGER :: n                            !<  number of nodes in element
     INTEGER :: NDOFs                        !< number of DOFs per node
     INTEGER :: DofIndexes(:)               !< Element node to global node numbering mapping
     REAL(KIND=dp), OPTIONAL, TARGET :: GlobalValues(:)
!------------------------------------------------------------------------------
     INTEGER :: i,j,k
     REAL(KIND=dp) :: s,t
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


     SELECT CASE( StiffMatrix % Format )
        CASE( MATRIX_CRS )
           CALL CRS_GlueLocalMatrix( StiffMatrix, &
                n, NDOFs, DofIndexes, LocalMassMatrix, GlobalValues )

!       CASE( MATRIX_LIST )
!          CALL List_GlueLocalMatrix( StiffMatrix % ListMatrix, &
!               n, NDOFs, DofIndexes, LocalMassMatrix )

!      CASE( MATRIX_BAND,MATRIX_SBAND )
!          CALL Band_GlueLocalMatrix( StiffMatrix, &
!               n, NDOFs, DofIndexes, LocalMassMatrix )

        CASE DEFAULT
          CALL FATAL( 'UpdateMassMatrix', 'Unexpected matrix format')
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE UpdateMassMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Get normal vector for solver using precomputed data.  
!------------------------------------------------------------------------------
   FUNCTION ConsistentNormalVector( Solver, NrmVar, Element, Found, Basis, node ) RESULT ( Normal ) 

     TYPE(Solver_t), POINTER :: Solver
     TYPE(Variable_t), POINTER :: NrmVar
     TYPE(Element_t), POINTER :: Element
     LOGICAL, OPTIONAL :: Found
     REAL(KIND=dp), OPTIONAL :: Basis(:)
     INTEGER, OPTIONAL :: node
     REAL(KIND=dp) :: Normal(3)
     
     TYPE(Solver_t), POINTER :: PrevSolver => NULL()
     TYPE(NormalTangential_t), POINTER :: NT => NULL()
     INTEGER :: NrmPerm(27),dofs,dim,i,j,m,n
     LOGICAL :: GotIt, uFound
     REAL(KIND=dp) :: NrmLen, NormalTest(3)
     
     SAVE PrevSolver, NT, dofs, dim
     
     uFound = .FALSE.
     Normal = 0.0_dp
     
     IF(.NOT. ASSOCIATED(Solver,PrevSolver) ) THEN
       dofs = 0
       
       !IF( .NOT. ASSOCIATED(NT) ) THEN
       !str = ListGetString( Solver % Values,'Normal Vector Name',GotIt )
       !IF(.NOT. GotIt) str = 'Normal Vector'         
       !NrmVar => VariableGet( Solver % Mesh % Variables, str, ThisOnly=.TRUE. ) 

       IF( ASSOCIATED(NrmVar) ) THEN
         ! If we have given Normal variable use that!
         dofs = NrmVar % Dofs
       END IF       

       ! If we have precomputed normal-tangential vector use that!
       NULLIFY(NT)       
       IF( ASSOCIATED(NT) ) THEN
         IF( NT % NormalTangentialNOFNodes == 0 ) NULLIFY(NT)
       END IF

       dim = CoordinateSystemDimension()
       PrevSolver => Solver
     END IF

     ! Note that we need to have full hit, otherwise return .FALSE.
     ! and use elemental normal vector in the code. 
     uFound = .FALSE.
     IF(dofs > 0) THEN
       IF(PRESENT( Node ) ) THEN
         j = NrmVar % Perm(node)
         IF( j>0 ) THEN
           Normal(1:dim) = NrmVar % Values(dofs*(j-1)+1:dofs*(j-1)+dim) 
           ! Some legacy ways to compute normal vector do not compute the vector even though
           ! the permutation is positive. This tries to deal with that kind of issue.
           NrmLen = SQRT(SUM(Normal(1:dim)**2))
           uFound = (NrmLen > 0.5_dp)
         END IF         
       ELSE IF( PRESENT( Basis ) ) THEN
         n = Element % TYPE % NumberOfNodes
         m = 0
         DO i=1,n
           j = NrmVar % Perm(Element % NodeIndexes(i))
           IF(j>0) THEN
             NormalTest(1:dim) = NrmVar % Values(dofs*(j-1)+1:dofs*(j-1)+dim) 
             NrmLen = SQRT(SUM(NormalTest(1:dim)**2))
             IF( NrmLen > 0.5_dp ) m = m+1
           END IF
         END DO
         IF( m == n ) THEN
           NrmPerm(1:n) = NrmVar % Perm( Element % NodeIndexes) - 1
           Normal(1) = SUM( Basis(1:n) * NrmVar % Values(dofs*NrmPerm(1:n)+1) ) 
           Normal(2) = SUM( Basis(1:n) * NrmVar % Values(dofs*NrmPerm(1:n)+2) )
           IF( dim == 3 ) THEN
             Normal(3) = SUM( Basis(1:n) * NrmVar % Values(dofs*NrmPerm(1:n)+3) ) 
           END IF
           uFound = .TRUE.
         END IF
       ELSE
         CALL Fatal('ConsistentNormalvector','Either Basis of Node is required!')
       END IF
     END IF
       
     ! We can also try to use the existing NT coordinate system associated normals.
     IF(.NOT. uFound .AND. ASSOCIATED(NT) ) THEN
       IF(PRESENT( Node ) ) THEN
         j = NT % BoundaryReorder(node)
         IF( j>0 ) THEN
           Normal(1:dim) = NT % BoundaryNormals(j,1:dim)
           uFound = .TRUE.
         END IF
       ELSE IF( PRESENT( Basis ) ) THEN
         n = Element % TYPE % NumberOfNodes       
         m = COUNT( NT % BoundaryReorder(Element % NodeIndexes) > 0 )
         IF( m == n ) THEN
           DO i=1,n
             j = NT % BoundaryReorder(Element % NodeIndexes(i))
             Normal(1:dim) = Normal(1:dim) + Basis(i) * NT % BoundaryNormals(j,1:dim)
           END DO
           uFound = .TRUE.
         END IF
       ELSE
         CALL Fatal('ConsistentNormalvector','Either Basis of Node is required!')
       END IF         
     END IF

     IF( uFound ) THEN
       NrmLen = SQRT(SUM(Normal**2))
       IF( ABS(1.0_dp-NrmLen) > 0.5_dp ) THEN         
         PRINT *,'NormalVector:',dofs,Element % ElementIndex, Element % NodeIndexes, Normal(1:dim)
         PRINT *,'Called by solver:',Solver % SolverId, ASSOCIATED(NT), ASSOCIATED(NrmVar), &
             PRESENT(Node), PRESENT(Basis)
         IF(PRESENT(Node)) THEN
           PRINT *,'Node:',Solver % Mesh % Nodes % x(node), Solver % Mesh % Nodes % y(node),&
               Solver % Mesh % Nodes % z(node), Node, nrmVar % Perm(node)
           PRINT *,'nrmVar:',SIZE(NrmVar % Values), Solver % Mesh % NumberOfNodes, &
               COUNT(NrmVar % Perm > 0)           
         END IF
         CALL Warn('ConsistentNormalVector','NormalVector should have a norm close to one!')
       END IF
     END IF       
     
     IF(PRESENT(Found)) Found = uFound

     
   END FUNCTION ConsistentNormalVector


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
     REAL(KIND=dp) :: LimitSign, EqSign, ValEps, LoadEps, ValEps0, LoadEps0, &
         MaxValue, MaxLoad, val
     INTEGER, POINTER :: FieldPerm(:), NodeIndexes(:)
     LOGICAL :: Found,AnyLimitBC, AnyLimitBF, GotInit, GotActive
     LOGICAL, ALLOCATABLE :: LimitDone(:)
     LOGICAL, POINTER :: LimitActive(:)
     TYPE(ValueList_t), POINTER :: Params, Entity
     LOGICAL, ALLOCATABLE :: InterfaceDof(:)
     INTEGER :: ConservativeAfterIters, NonlinIter, CoupledIter, DownStreamDirection
     LOGICAL :: Conservative, ConservativeAdd, ConservativeRemove, RelativeEps, &
         DoAdd, DoRemove, DirectionActive, FirstTime, DownStreamRemove
     TYPE(Mesh_t), POINTER :: Mesh

     CHARACTER(:), ALLOCATABLE :: Name,LimitName, InitName, ActiveName
     CHARACTER(*), PARAMETER :: Caller = 'DetermineSoftLimiter'

     Model => CurrentModel
     Var => Solver % Variable
     Mesh => Solver % Mesh
     
     ! Check the iterations counts and determine whether this is the first 
     ! time with this solver. 
     !------------------------------------------------------------------------
     FirstTime = .TRUE.
     iterV => VariableGet( Mesh % Variables,'nonlin iter')
     IF( ASSOCIATED( iterV ) ) THEN
       NonlinIter =  NINT( iterV % Values(1) ) 
       IF( NonlinIter > 1 ) FirstTime = .FALSE.
     END IF

     iterV => VariableGet( Mesh % Variables,'coupled iter')
     IF( ASSOCIATED( iterV ) ) THEN
       CoupledIter = NINT( iterV % Values(1) )
       IF( CoupledIter > 1 ) FirstTime = .FALSE.
     END IF
     IF( FirstTime ) THEN
       CALL Info(Caller,'Initializing soft limiter for solver',Level=10)
     END IF

     
     ! Determine variable for computing the contact load used to determine the 
     ! soft limit set.
     !------------------------------------------------------------------------
     CALL Info(Caller,'Determining soft limiter problems',Level=8)
     LoadVar => VariableGet( Model % Variables, &
         GetVarName(Var) // ' Contact Load',ThisOnly = .TRUE. )
     CALL CalculateLoads( Solver, Solver % Matrix, Var % Values, Var % DOFs, .FALSE., LoadVar ) 

     IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
       CALL Fatal(Caller, &
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
         CALL Info(Caller,'Adding dofs in conservative fashion',Level=8)
       END IF
     END IF

     ConservativeRemove = .FALSE.
     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Remove After Iterations',Found )      
     IF( Found ) THEN
       Conservative = .TRUE.  
       ConservativeRemove = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeRemove ) THEN
         CALL Info(Caller,'Removing dofs in conservative fashion',Level=8)
       END IF
     END IF

     DownStreamRemove = ListGetLogical( Params,'Apply Limiter Remove Downstream',Found)
     IF( DownStreamRemove ) THEN
       CALL Info(Caller,'Removing contact dofs only in downstream',Level=8)      
       ConservativeRemove = .TRUE.
       Conservative = .TRUE.
       DownStreamDirection = ListGetInteger( Params,'Apply Limiter Downstream Direction',Found)
       IF(.NOT. Found ) DownStreamDirection = 1
     END IF
       
     LoadEps0 = ListGetConstReal(Params,'Limiter Load Tolerance',Found ) 
     IF(.NOT. Found ) LoadEps0 = 1.0e-8_dp
     LoadEps = LoadEps0
     
     ValEps0 = ListGetConstReal(Params,'Limiter Value Tolerance',Found ) 
     IF(.NOT. Found ) ValEps0 = 1.0e-8_dp
     ValEps = ValEps0

     RelativeEps = ListGetLogical(Params,'Limiter Relative Tolerance', Found ) 
     
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

         Name = TRIM(Var % name)
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
         
         CALL Info(Caller,'Applying limit: '//TRIM(LimitName),Level=8)

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
           IF(ASSOCIATED(Var % LowerLimitActive)) THEN
             IF(SIZE(Var % LowerLimitActive) /= totsize) THEN
               DEALLOCATE(Var % LowerLimitActive)
             END IF
           END IF
           IF( .NOT. ASSOCIATED(Var % LowerLimitActive ) ) THEN
             ALLOCATE( Var % LowerLimitActive( totsize ) )
             Var % LowerLimitActive = .FALSE.
           END IF
           LimitActive => Var % LowerLimitActive
         ELSE
           IF(ASSOCIATED(Var % UpperLimitActive)) THEN
             IF(SIZE(Var % UpperLimitActive) /= totsize) THEN
               DEALLOCATE(Var % UpperLimitActive)
             END IF
           END IF
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

         IF( RelativeEps ) THEN         
           IF( dofs == 1 ) THEN
             MaxLoad = MAXVAL( ABS( LoadValues ) )
             MaxValue = MAXVAL( ABS( Var % Values ) ) 
           ELSE
             MaxLoad = MAXVAL( ABS( LoadValues(dof::dofs ) ) )
             MaxValue = MAXVAL( ABS( Var % Values(dof::dofs) ) ) 
           END IF
           MaxLoad = ParallelReduction(MaxLoad)
           MaxValue = ParallelReduction(MaxValue)
           IF( InfoActive(20) ) THEN
             WRITE(Message,'(A,ES12.3)') 'Using load maximum for scaling: ',MaxLoad
             CALL Info(Caller,Message)
             WRITE(Message,'(A,ES12.3)') 'Using value maximum for scaling: ',MaxValue
             CALL Info(Caller,Message)
           END IF
           LoadEps = LoadEps0 * MaxLoad
           ValEps = ValEps0 * MaxValue 
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
               IF(Element % BodyId == 0) CYCLE
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
                 IF(.NOT. LimitActive(ind)) THEN
                   added = added + 1
                   LimitActive(ind) = .TRUE.
                 END IF
               ELSE IF( GotActive .AND. ElemActive(i) > 0.0_dp ) THEN
                 IF(.NOT. LimitActive(ind)) THEN
                   added = added + 1
                   LimitActive(ind) = .TRUE.
                 END IF
               ELSE
                 LimitActive(ind) = .FALSE.
               END IF

               ! Enforce the values to limits because nonlinear material models
               ! may otherwise lead to divergence of the iteration
               !--------------------------------------------------------------
               IF( LimitActive(ind) ) THEN
                 IF( Upper == 0 ) THEN
                   Var % Values(ind) = MAX( Var % Values(ind), ElemLimit(i) )
                 ELSE
                   Var % Values(ind) = MIN( Var % Values(ind), ElemLimit(i) )
                 END IF
               END IF

               LimitDone(ind) = .TRUE.             
             END DO
           END DO

           CYCLE
         END IF ! FirstTime


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
             
             IF( DownStreamRemove ) THEN
               ! This includes only interface dofs donwstream from
               ! non-contact zone.
               BLOCK
                 REAL(kind=DP) :: r1(3),r2(3),dr(3),reps=1.0d-6
                 
                 DO i=1,n
                   j = FieldPerm( NodeIndexes(i) )
                   IF( j == 0 ) CYCLE
                   ind = Dofs * ( j - 1) + Dof
                   
                   ! Downstream of non-contact zone
                   IF(LimitActive(ind)) CYCLE
                                      
                   DO i2 = i,n
                     IF( i2 == i ) CYCLE                   
                     j2 = FieldPerm( NodeIndexes(i2) )
                     IF( j2 == 0 ) CYCLE
                     ind2 = Dofs * ( j2 - 1) + Dof
                     
                     IF( LimitActive(ind2) ) THEN
                       r2(1) =  Mesh % Nodes % x(NodeIndexes(i2))
                       r2(2) =  Mesh % Nodes % y(NodeIndexes(i2))
                       r2(3) =  Mesh % Nodes % z(NodeIndexes(i2))
                       
                       r1(1) = Mesh % Nodes % x(NodeIndexes(i))
                       r1(2) = Mesh % Nodes % y(NodeIndexes(i))
                       r1(3) = Mesh % Nodes % z(NodeIndexes(i))

                       k = DownStreamDirection 
                       IF( k > 0 ) THEN
                         dr = r2 - r1
                       ELSE
                         dr = r1 - r2
                         k = -k
                       END IF
                       
                       IF( dr(k) < reps ) CYCLE
                       
                       IF( dr(k) > 0.5*SQRT(SUM(dr*dr)) ) THEN
                         InterfaceDof(ind2) = .TRUE.
                         !PRINT *,'downstream coord:',dr
                       END IF
                     END IF
                   END DO
                 END DO
               END BLOCK
             ELSE
               ! This includes all interface dofs
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
             END IF
           END DO

           CALL Info(Caller,&
               'Number of interface dofs: '//I2S(COUNT(InterfaceDof)),Level=8)
         END IF

         IF( DownStreamRemove ) THEN
           t = COUNT(InterfaceDof)
           CALL Info(Caller,'Downstream contact set dofs:'//I2S(t),Level=8)
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
             IF(Element % BodyId == 0) CYCLE
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
                   IF(LimitActive(ind)) THEN                     
                     removed = removed + 1
                     LimitActive(ind) = .FALSE.
                   END IF
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

             IF( LimitActive(ind) ) THEN
               ! Set the Dirichlet conditions already here!
               Solver % Matrix % DValues(ind) = ElemLimit(i)
               Solver % Matrix % ConstrainedDOF(ind) = .TRUE.

               ! Enforce the values to limits because nonlinear material models
               ! may otherwise lead to divergence of the iteration
               IF( Upper == 0 ) THEN
                 Var % Values(ind) = MAX( Var % Values(ind), ElemLimit(i) )
               ELSE
                 Var % Values(ind) = MIN( Var % Values(ind), ElemLimit(i) )
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
           WRITE(Message,'(A,I0)') 'Number of lower limited dofs for '&
               //TRIM(GetVarName(Var))//': ',COUNT( LimitActive )
         ELSE
           WRITE(Message,'(A,I0)') 'Number of upper limited dofs for '&
               //TRIM(GetVarName(Var))//': ',COUNT( LimitActive )
         END IF
         CALL Info(Caller,Message,Level=5)
         
         IF(added + removed >= 0) THEN
           CALL Info(Caller,'Added '//I2S(added)//' and removed '&
               //I2S(removed)//' dofs in contact set',Level=6)
         END IF
       END IF
         
     END DO
                
     ! Optionally save the limiters as a field variable so that 
     ! lower limit is given value -1.0 and upper limit value +1.0.
     IF( ListGetLogical( Params,'Save Limiter',Found ) ) THEN
       
       LimitVar => VariableGet( Model % Variables, &
           GetVarName(Var) // ' Contact Active',ThisOnly = .TRUE. )
       IF(.NOT. ASSOCIATED( LimitVar ) ) THEN
         CALL Info(Caller,'Creating field for contact: '//TRIM(GetVarName(Var)),Level=7)
         CALL VariableAddVector( Model % Variables, Solver % Mesh, Solver,&
             GetVarName(Var) //' Contact Active', Var % Dofs, Perm = FieldPerm )
         LimitVar => VariableGet( Model % Variables, &
             GetVarName(Var) // ' Contact Active',ThisOnly = .TRUE. )
       END IF
              
       LimitVar % Values = 0.0_dp
       IF( ASSOCIATED( Var % LowerLimitActive ) ) THEN
         IF(SIZE(LimitVar % Values) /= SIZE(Var % LowerLimitActive)) THEN
           CALL Fatal(Caller,'Mismatch in size for LimitVar values!')
         END IF
         WHERE( Var % LowerLimitActive ) 
           LimitVar % Values = -1.0_dp
         END WHERE           
       END IF
       IF( ASSOCIATED( Var % UpperLimitActive ) ) THEN
         IF(SIZE(LimitVar % Values) /= SIZE(Var % UpperLimitActive)) THEN
           CALL Fatal(Caller,'Mismatch in size for LimitVar values!')
         END IF
         WHERE( Var % UpperLimitActive ) 
           LimitVar % Values = 1.0_dp
         END WHERE           
       END IF
     END IF

     IF( ALLOCATED( LimitDone ) ) THEN
       DEALLOCATE( LimitDone, ElemLimit, ElemInit, ElemActive ) 
     END IF
     
     IF( ALLOCATED( InterfaceDof ) ) THEN
       DEALLOCATE( InterfaceDof )
     END IF

     CALL Info(Caller,'All done',Level=12)

   END SUBROUTINE DetermineSoftLimiter
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Gets a string from the list by its name, if not found return empty string.
!------------------------------------------------------------------------------
   FUNCTION LagrangeMultiplierName( Solver, SetUnfound ) RESULT( Name )
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver
     LOGICAL, OPTIONAL :: SetUnfound
     CHARACTER(:), ALLOCATABLE :: Name
!------------------------------------------------------------------------------
     LOGICAL :: Found     
     Name = ListGetString( Solver % Values,'Lagrange Multiplier Name', Found )
     IF(.NOT. Found ) THEN
       Name = 'Lagrange Multiplier '//TRIM(ListGetString(Solver % Values,'equation'))
       IF(PRESENT(SetUnfound)) THEN
         IF(SetUnfound ) THEN 
           CALL Info('LagrangeMultiplierName','Defaulting name to: '//TRIM(Name))
           CALL ListAddString( Solver % Values,'Lagrange Multiplier Name', TRIM(Name) ) 
         END IF
       END IF
     END IF
     
   END FUNCTION LagrangeMultiplierName
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
         WeightVar, NormalActiveVar, StickActiveVar, GapVar, ContactLagrangeVar
     TYPE(Element_t), POINTER :: Element
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: i,j,k,l,n,m,t,ind,dofs,cdofs,dim, bf, Upper, &
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
     CHARACTER(:), ALLOCATABLE ::  Str,VarName, ContactType
     INTEGER :: ConservativeAfterIters, ActiveDirection, NonlinIter, CoupledIter, EasyIters
     LOGICAL :: ConservativeAdd, ConservativeRemove, &
         DoAdd, DoRemove, DirectionActive, Rotated, FlatProjector, PlaneProjector, &
         RotationalProjector, NormalProjector, FirstTime = .TRUE., &
         AnyRotatedContact, ThisRotatedContact, StickContact, TieContact, FrictionContact, SlipContact, &
         CalculateVelocity, NodalNormal, ResidualMode, AddDiag, SkipFriction, DoIt
     TYPE(MortarBC_t), POINTER :: MortarBC
     TYPE(Matrix_t), POINTER :: Projector, DualProjector
     TYPE(ValueList_t), POINTER :: BC, MasterBC
     REAL(KIND=dp), POINTER :: nWrk(:,:)
     LOGICAL :: CreateDual, pContact
     INTEGER, TARGET :: pIndexes(12)
     TYPE(Variable_t), POINTER :: UseLoadVar
     LOGICAL :: UseLagrange
     TYPE(NormalTangential_t), POINTER :: NT
     CHARACTER(*), PARAMETER :: Caller = 'DetermineContact'

     
     SAVE FirstTime

     CALL Info(Caller,'Setting up contact conditions',Level=8)
     
     Model => CurrentModel
     Var => Solver % Variable
     VarName = GetVarName( Var ) 
     Mesh => Solver % Mesh
     NT => Model % Solver % NormalTangential
     
     ! Is any boundary rotated or not
     AnyRotatedContact = ( NT % NormalTangentialNOFNodes > 0 ) 

     ! The variable to be constrained by the contact algorithm
     ! Here it is assumed to be some "displacement" i.e. a vector quantity
     FieldValues => Var % Values
     FieldPerm => Var % Perm
     totsize = SIZE( FieldValues )
     dofs = Var % Dofs
     dim = Mesh % MeshDim 
     Params => Solver % Values

     IF(dofs == dim) THEN
       cdofs = dim
     ELSE IF(dofs == dim+1) THEN
       CALL Info(Caller,'We seem to have mixed formulation, ignoring pressure!',Level=7)
       cdofs = dim
     ELSE 
       CALL Fatal(Caller,'Invalid number of dofs for contact problem: '//I2S(dofs))
     END IF     

     pContact = IsPelement(Mesh % Elements(1) )
     IF( pContact ) THEN
       ! We only have to deal with the middle dofs if they are not condensated away!
       IF( .NOT. ListGetLogical( Params,'Bubbles in Global System',Found ) ) THEN
         IF(Found) pContact = .FALSE.   
       END IF
     END IF
     IF( ListGetLogical( Params,'Contact Linear Basis',Found ) ) THEN
       pContact = .FALSE.
     END IF
    
     IF( pContact ) THEN
       CALL Info(Caller,'Using p-elements for contact, if available in projector!',Level=8)
     END IF
     
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
         CALL Info(Caller,'Adding dofs in conservative fashion',Level=8)
       END IF
     END IF

     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Remove After Iterations',ConservativeRemove ) 
     IF( ConservativeRemove ) THEN
       IF( CoupledIter == 1 ) ConservativeRemove = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeRemove ) THEN
         CALL Info(Caller,'Removing dofs in conservative fashion',Level=8)
       END IF
     END IF
         
     EasyIters = ListGetInteger(Params,'Nonlinear System Initial Friction Iterations',Found )
     IF(.NOT. Found ) EasyIters = 1          

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
       CALL Fatal(Caller,'Cannot apply contact without projectors!')
     END IF

     ! a) Create rotateted contact if needed
     CALL RotatedDisplacementField() 

     CALL PickLagrangeMultiplier()

     ! b) Create and/or obtain pointers to boundary variables 
     CALL GetContactFields( FirstTime )

     ! c) Calculate the contact loads to the normal direction
     LoadVar => CalculateContactLoad() 
     LoadValues => LoadVar % Values

     UseLagrange = ListGetLogical( Params,'Use Lagrange Multiplier for Contact',Found )
     IF( UseLagrange ) THEN
       CALL Info(Caller,'Using Lagrange multiplier to determine contact condition!')
       UseLoadVar => ContactLagrangeVar 
     ELSE
       UseLoadVar => NormalLoadVar 
     END IF     
     
     ! Loop over each contact pair
     !--------------------------------------------------------------
     DO bc_ind = 1, Model % NumberOfBCs
       
       MortarBC => Model % Solver % MortarBCs(bc_ind)  
       IF( .NOT. ASSOCIATED( MortarBC ) ) CYCLE

       Projector => MortarBC % Projector
       IF(.NOT. ASSOCIATED(Projector) ) CYCLE

       BC => Model % BCs(bc_ind) % Values

       CALL Info(Caller,'Set contact for boundary: '&
           //I2S(bc_ind),Level=8)
       Model % Solver % MortarBCsChanged = .TRUE.
       
       FlatProjector = ListGetLogical( BC, 'Flat Projector',Found ) 
       PlaneProjector = ListGetLogical( BC, 'Plane Projector',Found )
       RotationalProjector = ListGetLogical( BC, 'Rotational Projector',Found ) .OR. &
           ListGetLogical( BC, 'Cylindrical Projector',Found )
       NormalProjector = ListGetLogical( BC, 'Normal Projector',Found )
       
       ! Is the current boundary rotated or not
       ThisRotatedContact = ListGetLogical( BC,'Normal-Tangential '//TRIM(VarName),Found)

       IF( FlatProjector ) THEN
         ActiveDirection = ListGetInteger( BC, 'Flat Projector Coordinate',Found )
         IF( .NOT. Found ) ActiveDirection = cdofs       
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
           CALL Info(Caller,'Active direction set to: '//I2S(ActiveDirection),Level=6)
         END IF
       ELSE IF( RotationalProjector .OR. NormalProjector ) THEN
         ActiveDirection = 1
         IF( .NOT. ThisRotatedContact ) THEN
           CALL Warn(Caller,'Rotational and normal projectors should only work with N-T coordinates!')
         END IF
       ELSE
         CALL Fatal(Caller,'Projector must be current either flat, plane, cylindrical or rotational!')
       END IF
      

       ! Get the pointer to the other side i.e. master boundary  
       master_ind = ListGetInteger( BC,'Mortar BC',Found )
       IF( .NOT. Found ) master_ind = ListGetInteger( BC,'Contact BC',Found )
       MasterBC => Model % BCs(master_ind) % Values
       
       ! If we have dual projector we may use it to map certain quantities directly to master nodes
       DualProjector => Projector % Ematrix
       CreateDual = ASSOCIATED( DualProjector )
       IF( CreateDual ) THEN
         CALL Info(Caller,'Using also the dual projector',Level=8)
       END IF
      
       ! If we have N-T system then the mortar condition for the master side
       ! should have reverse sign as both normal displacement diminish the gap.
       IF( ThisRotatedContact ) THEN
         IF( master_ind > 0 ) THEN
           IF( .NOT. ListGetLogical( MasterBC, &
               'Normal-Tangential '//TRIM(VarName),Found) ) THEN
             CALL Fatal(Caller,'Master boundary '//I2S(master_ind)//&
                 ' should also have N-T coordinates!')
           END IF
         END IF

         CALL Info(Caller,'We have a normal-tangential system',Level=6)
         MortarBC % MasterScale = -1.0_dp
         DofN = 1
       ELSE                 
         DofN = ActiveDirection 
       END IF

       ! Get the degrees of freedom related to the normal and tangential directions
       DofT1 = 0; DofT2 = 0
       DO i=1,cdofs
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
           CALL Fatal(Caller,'Unknown contact type: '//TRIM(ContactType))
         END SELECT
       ELSE
         StickContact = ListGetLogical( BC,'Stick Contact',Found )
         IF(.NOT. Found ) TieContact = ListGetLogical( BC,'Tie Contact',Found )
         IF(.NOT. Found ) FrictionContact = ListGetLogical( BC,'Friction Contact',Found )
         IF(.NOT. Found ) SlipContact = ListGetLogical( BC,'Slip Contact',Found )
         IF(.NOT. Found ) SlipContact = ListGetLogical( BC,'Slide Contact',Found )
         IF(.NOT. Found ) THEN 
           CALL Warn(Caller,'No contact type given, assuming > Slip Contact <')
           SlipContact = .TRUE.
         END IF
       END IF

       IF( StickContact ) CALL Info(Caller,'Using "Stick Contact" for displacement',Level=10)
       IF( TieContact ) CALL Info(Caller,'Using "Tie Contact" for displacement',Level=10)
       IF( FrictionContact ) CALL Info(Caller,'Using "Friction Contact" for displacement',Level=10)
       IF( SlipContact ) CALL Info(Caller,'Using "Slip Contact" for displacement',Level=10)
       

       ! At the start of nonlinear or transient iteration it may be beneficial to assume
       ! less challenging initial contat model
       ! friction models -> stick contact /slip contact -> tie contact
       BLOCK
         LOGICAL :: InitialTie, InitialStick, InitialSlip 

         InitialTie = .FALSE.
         InitialStick = .FALSE.
         InitialSlip = .FALSE.
         
         DO i=1,3
           IF(i==1) THEN
             Str = 'Nonlinear System'
             IF( NonlinIter > EasyIters ) CYCLE
           ELSE IF(i==2) THEN
             Str = 'Transient System'
             IF( TimeStep > 1 ) CYCLE 
           ELSE
             Str = ''
             IF( NonlinIter > 1 .OR. Timestep > 1 ) CYCLE
           END IF

           IF( ListGetLogical( BC,TRIM(Str)//' Initial Tie Contact', Found ) ) InitialTie = .TRUE.
           IF( ListGetLogical( BC,TRIM(Str)//' Initial Stick Contact', Found ) ) InitialStick = .TRUE.
           IF( ListGetLogical( BC,TRIM(Str)//' Initial Slip Contact', Found ) ) InitialSlip = .TRUE.
         END DO
         
         SkipFriction = .FALSE.
         IF( InitialTie .OR. InitialStick .OR. InitialSlip ) THEN
           FrictionContact = .FALSE.
           StickContact = .FALSE.
           SlipContact = .FALSE.
           TieContact = .FALSE.
           SkipFriction = .TRUE.
         END IF

         IF( InitialTie ) THEN
           TieContact = .TRUE.
           CALL Info(Caller,'Reverting to "Tie Contact" for the start',Level=10)                      
         ELSE IF(InitialStick) THEN
           StickContact = .TRUE.
           CALL Info(Caller,'Reverting to "Stick Contact" for the start',Level=10)                      
         ELSE IF(InitialSlip ) THEN
           SlipContact = .TRUE.
           CALL Info(Caller,'Reverting to "Slip Contact" for the start',Level=10)                      
         END IF
       END BLOCK

       ! If we have stick contact then create a diagonal entry to the projection matrix.
       IF( StickContact .OR. FrictionContact ) THEN
         AddDiag = ListCheckPresent( BC,'Stick Contact Coefficient')      
       ELSE
         AddDiag = .FALSE.
       END IF

       IF(InfoActive(30)) THEN
         PRINT *,'Contact Flags:',TieContact, FrictionContact, StickContact, SlipContact, SkipFriction, AddDiag 
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

       ! i) If requested ensure a minimum number of contact nodes
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
     CALL Info(Caller,'All done',Level=10)

   CONTAINS


     ! Given the cartesian solution compute the rotated solution.
     !-------------------------------------------------------------------------
     SUBROUTINE RotatedDisplacementField( ) 

       REAL(KIND=dp) :: RotVec(3)
       INTEGER :: i,j,k,n,m

       IF( .NOT. AnyRotatedContact ) RETURN

       CALL Info(Caller,'Rotating displacement field',Level=8)
       ALLOCATE( RotatedField(Solver % Matrix % NumberOfRows ) )
       RotatedField = Var % Values

       n = SIZE( FieldPerm ) 
       m = SIZE( NT % BoundaryReorder )
       IF( n > m ) THEN
         i = COUNT(FieldPerm(m+1:n) > 0 )
         IF( i > 0 ) THEN
           CALL Fatal(Caller,'Number of potential untreated rotations: '//I2S(i))
         END IF
       END IF
       
       DO i=1,SIZE(FieldPerm)
         j = FieldPerm(i)
         IF( j == 0 ) CYCLE
         m = NT % BoundaryReorder(i)
         IF( m == 0 ) CYCLE
         
         RotVec = 0._dp
         DO k=1,cdofs
           RotVec(k) = RotatedField(dofs*(j-1)+k)
         END DO
         CALL RotateNTSystem( RotVec, i )
         DO k=1,cdofs
           RotatedField(dofs*(j-1)+k) = RotVec( k )
         END DO
       END DO

       IF( InfoActive(30) ) THEN
         CALL VectorValuesRange(RotatedField, SIZE(RotatedField),'RotatedField')
       END IF
       
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


       CALL Info(Caller,'Determining reaction forces for contact problem',Level=10)

       LoadVar => VariableGet( Model % Variables, &
           TRIM(VarName) // ' Contact Load',ThisOnly = .TRUE. )
       IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
         CALL Fatal(Caller, &
             'No Loads associated with variable: '//GetVarName(Var) )
       END IF

       IF( AnyRotatedContact ) THEN
         TempX => RotatedField 
       ELSE
         TempX => FieldValues
       END IF

       CALL CalculateLoads( Solver, Solver % Matrix, TempX, dofs, .FALSE., LoadVar ) 

       IF( InfoActive(30) ) THEN
         CALL VectorValuesRange(LoadVar % Values, SIZE(LoadVar % Values),'ContactLoad')
       END IF
       
     END FUNCTION CalculateContactLoad

     
     ! Given the previous solution and the related Lagrange multiplier pick the
     ! new multiplier such that it may be visualized as a field.
     !-------------------------------------------------------------------------
     SUBROUTINE PickLagrangeMultiplier( ) 

       TYPE(Variable_t), POINTER :: LinSysVar, ContactSysVar, ActiveVar
       INTEGER :: i,j,k,l,n
       INTEGER, POINTER :: InvPerm(:)
       
       CALL Info(Caller,'Pick lagrange coefficient from the active set to whole set',Level=10)

       LinSysVar => VariableGet( Model % Variables, &
           LagrangeMultiplierName(Solver),ThisOnly = .TRUE. )
       IF( .NOT. ASSOCIATED( LinSysVar ) ) THEN
         CALL Warn(Caller, &
             'No Lagrange multiplier field associated with linear system: '//GetVarName(Var) )
         RETURN
       END IF
       
       ContactSysVar => VariableGet( Model % Variables, &
           TRIM(VarName) // ' Lagrange Multiplier',ThisOnly = .TRUE. )
       IF( .NOT. ASSOCIATED( ContactSysVar ) ) THEN
         CALL Fatal(Caller, &
             'No Lagrange multiplier field associated with: '//GetVarName(Var) )
       END IF
       ContactSysVar % Values = 0.0_dp

       IF(.NOT. ASSOCIATED( Solver % Matrix % ConstraintMatrix ) ) THEN
          CALL Fatal(Caller, &
             'No constraint matrix associated with: '//GetVarName(Var) )          
       END IF            
       
       InvPerm => Solver % Matrix % ConstraintMatrix % InvPerm
       n = Solver % Matrix % ConstraintMatrix % NumberOfRows
       
       DO i=1,SIZE(InvPerm)
         ! This is related to the full matrix equation
         j = InvPerm(i)

         IF( MODULO(j,dofs) /= 1 ) CYCLE
         l = (j-1)/dofs+1
         
         IF( l > 0 .AND. l <= SIZE( ContactSysVar % Perm ) ) THEN
           k = ContactSysVar % Perm(l)
           IF( k > 0 ) THEN
             ContactSysVar % Values(k) = LinSysVar % Values(i)
           END IF
         END IF
       END DO

       IF(InfoActive(30)) THEN
         CALL VectorValuesRange( LinsysVar % Values,SIZE(LinsysVar % Values),'LinsysValues')
         CALL VectorValuesRange( ContactSysVar % Values,SIZE(ContactSysVar % Values),'ContactSysValues')
       END IF
                     
     END SUBROUTINE PickLagrangeMultiplier

     
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
         CALL Info(Caller,'Creating contact fields',Level=8)

         n = SIZE( FieldPerm ) 
         ALLOCATE( BoundaryPerm(n) )
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
                 IF( pContact ) THEN
                   n = mGetElementDOFs(pIndexes,Element)                   
                   BoundaryPerm(pIndexes(1:n)) = 1
                 ELSE                 
                   BoundaryPerm( Element % NodeIndexes ) = 1
                 END IF
               END IF
             END IF
           END DO
         END DO

         DEALLOCATE( ActiveBCs )

         j = 0
         DO i=1,SIZE(BoundaryPerm)
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
               TRIM(VarName)//' Contact Velocity',cDofs,Perm = BoundaryPerm )
         END IF
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Lagrange Multiplier',1,Perm = BoundaryPerm )
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
       
       ContactLagrangeVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Lagrange Multiplier')       
       
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
       totsize = cDofs * onesize

       IF( .NOT. AddDiag .AND. ASSOCIATED(MortarBC % Diag) ) THEN
         DEALLOCATE( MortarBC % Diag ) 
       END IF

       ! Create the permutation that is later needed in putting the diag and rhs to correct position.
       ALLOCATE( Perm( SIZE( FieldPerm ) ) )
       Perm = 0
       DO i=1,SIZE( Projector % InvPerm )
         j = Projector % InvPerm(i) 
         IF( j == 0 ) CYCLE
         Perm( j ) = i
       END DO

       ! First time nothing is allocated
       IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN
         CALL Info(Caller,'Allocating projector mortar vectors of size: '//I2S(totsize),Level=10)
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
       SamePerm = ALL( Perm == MortarBC % Perm )

       ! Permutation unchanged, just return
       IF( SamePerm ) THEN
         DEALLOCATE( Perm ) 
         RETURN
       END IF
       
       ! Permutation changes, and also sizes changed?
       SameSize = ( SIZE(MortarBC % Rhs) == totsize )
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

         DO l=1,cDofs
           Active(cDofs*(j-1)+l) = MortarBC % Active(cDofs*(k-1)+l)
         END DO
       END DO

       IF( ASSOCIATED( MortarBC % Active ) ) DEALLOCATE( MortarBC % Active )
       IF( ASSOCIATED( MortarBC % Perm ) ) DEALLOCATE( MortarBC % Perm )
       MortarBC % Active => Active 
       MortarBC % Perm => Perm 

       IF( AddDiag ) THEN
         IF( ASSOCIATED( MortarBC % Diag ) ) DEALLOCATE( MortarBC % Diag ) 
         MortarBC % Diag => Diag 
       END IF

       CALL Info(Caller,'Copied > Active < flag to changed projector',Level=8)

     END SUBROUTINE InitializeMortarVectors

     

     ! Make a list of interface dofs to allow conservative algorithms. 
     ! There only nodes that are at the interface are added or removed from the set.
     !------------------------------------------------------------------------------
     SUBROUTINE MarkInterfaceDofs()
       
       INTEGER :: i,j,i2,j2,k,k2,l,n,ind,ind2,elem
       INTEGER, POINTER :: Indexes(:)
       TYPE(Element_t), POINTER :: Element
       
       CALL Info(Caller,'Marking interface dofs for conservative adding/removal',Level=8)

       IF(.NOT. ALLOCATED( InterfaceDof ) ) THEN
         ALLOCATE( InterfaceDof( SIZE(MortarBC % Active) ) )
       END IF
       InterfaceDof = .FALSE. 


       DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         
         Element => Mesh % Elements( elem )         
         IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE

         IF( pContact ) THEN
           n = mGetElementDOFs(pIndexes,Element)                   
           Indexes => pIndexes
         ELSE         
           n = Element % TYPE % NumberOfNodes         
           Indexes => Element % NodeIndexes
         END IF
           
         DO i=1,n
           j = FieldPerm( Indexes(i) )
           IF( j == 0 ) CYCLE
           k = MortarBC % Perm( Indexes(i) )
           
           DO i2 = i+1,n
             j2 = FieldPerm( Indexes(i2) )
             IF( j2 == 0 ) CYCLE
             k2 = MortarBC % perm( Indexes(i2) )
             
             DO l=1,cDofs             
               ind = cDofs * ( k - 1 ) + l
               ind2 = cDofs * ( k2 - 1) + l
               
               IF( MortarBC % Active(ind) .NEQV. MortarBC % Active(ind2) ) THEN
                 InterfaceDof(ind) = .TRUE.
                 InterfaceDof(ind2) = .TRUE.
               END IF
             END DO
           END DO
         END DO
       END DO

       n = COUNT(InterfaceDof)
       CALL Info(Caller,'Number of interface dofs: '//I2S(n),Level=8)
       
     END SUBROUTINE MarkInterfaceDofs
     

     ! Calculates the signed distance that is used to define whether we have contact or not.
     ! If distance is negative then we can later add the corresponding node to the contact set
     ! Also computes the right-hand-side of the mortar equality constrained which is the 
     ! desired distance in the active direction. Works also for residual mode which greatly 
     ! improves the convergence for large displacements.  
     !----------------------------------------------------------------------------------------
     SUBROUTINE CalculateMortarDistance()

       REAL(KIND=dp) :: Disp(3), Coord(3), PrevDisp(3), Velo(3), ContactVec(3), ContactVelo(3), &
           LocalNormal0(3), SlipCoord(3), CartVec(3), ContactDist
       REAL(KIND=dp), POINTER :: DispVals(:), PrevDispVals(:) 
       REAL(KIND=dp) :: MinDist, MaxDist, wsum, wsumM, mult
       TYPE(Matrix_t), POINTER :: ActiveProjector
       LOGICAL :: IsSlave, IsMaster, DistanceSet
       LOGICAL, ALLOCATABLE :: SlaveNode(:), MasterNode(:), NodeDone(:)
       INTEGER, POINTER :: Indexes(:)
       INTEGER :: elemcode, CoeffSign
       REAL(KIND=dp), ALLOCATABLE :: CoeffTable(:)
       INTEGER :: l2,elem,i1,i2,j1,j2,n
       LOGICAL :: LinearContactGap, DebugNormals
       
       
       CALL Info('CalculateMortarDistance','Computing distance between mortar boundaries',Level=14)

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
         IF(.NOT. ASSOCIATED( PrevDispVals ) ) CALL Fatal('CalculateMortarDistance',&
             'Previous displacement field required!')
       END IF

       LinearContactGap = ListGetLogical( Model % Simulation,&
           'Contact BCs linear gap', Found )      

       ALLOCATE( SlaveNode( SIZE( FieldPerm ) ) )
       SlaveNode = .FALSE.

       IF( CreateDual ) THEN
         ALLOCATE( MasterNode( SIZE( FieldPerm ) ) )
         MasterNode = .FALSE.
       END IF

       DO i=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( i )                  
         IF( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) THEN
           IF( pContact ) THEN
             n = mGetElementDOFs(pIndexes,Element)                   
             SlaveNode(pIndexes(1:n)) = .TRUE.
           ELSE
             SlaveNode( Element % NodeIndexes ) = .TRUE.
           END IF
         END IF
         IF( CreateDual ) THEN
           IF ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) THEN
             IF( pContact ) THEN
               n = mGetElementDOFs(pIndexes,Element)                   
               MasterNode(pIndexes(1:n)) = .TRUE.
             ELSE               
               MasterNode( Element % NodeIndexes ) = .TRUE.
             END IF
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

       IF( DebugNormals ) THEN
         PRINT *,'Flags:',TieContact,ResidualMode,ThisRotatedContact,NodalNormal,StickContact,RotationalProjector
       END IF

       
100    CONTINUE

       DO i = 1,ActiveProjector % NumberOfRows

         j = ActiveProjector % InvPerm(i)

         IF( j == 0 ) CYCLE
         
         wsum = 0.0_dp
         wsumM = 0.0_dp
         Dist = 0.0_dp
         DistN = 0.0_dp
         DistT1 = 0.0_dp
         DistT2 = 0.0_dp
         ContactVelo = 0.0_dp
         ContactVec = 0.0_dp
         DistanceSet = .FALSE.
         ContactDist = 0.0_dp
         CartVec = 0.0_dp
         
         ! This is the most simple contact condition. We just want no slip on the contact.
         IF( TieContact .AND. .NOT. ResidualMode ) GOTO 200

         ! Get the normal of the slave surface.
         IF( ThisRotatedContact ) THEN
           Rotated = GetSolutionRotation(NTT, j )
           LocalNormal = NTT(:,1)
           LocalNormal0 = LocalNormal
           LocalT1 = NTT(:,2)
           IF( cDofs == 3 ) LocalT2 = NTT(:,3)
         ELSE
           LocalNormal = ContactNormal
           LocalT1 = ContactT1
           IF( cDofs == 3 ) LocalT2 = ContactT2 
         END IF

         ! Compute normal of the master surface from the average sum of normals
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
             IF( cDofs == 3 ) LocalT2 = LocalT2 + coeff * NTT(:,3)
           END DO

           ! Normalize the unit vector length to one
           LocalNormal = LocalNormal / SQRT( SUM( LocalNormal**2 ) )
           LocalT1 = LocalT1 / SQRT( SUM( LocalT1**2 ) )
           IF( cDofs == 3 ) LocalT2 = LocalT2 / SQRT( SUM( LocalT1**2 ) )

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

         ! Compute the weigted distance in the normal direction.
         DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
           k = ActiveProjector % Cols(j)

           l = FieldPerm( k ) 
           IF( l == 0 ) CYCLE

           coeff = ActiveProjector % Values(j)
                  
           ! Only compute the sum related to the active projector
           IF( SlaveNode(k) ) THEN
             wsum = wsum + coeff
           ELSE
             wsumM = wsumM + coeff
           END IF           
         END DO

         IF( ABS( wsum ) <= TINY( wsum ) ) THEN
           CALL Fatal('CalculateMortarDistance','wsum seems to be almost zero!')
         END IF
         IF( ABS( wsumM ) <= TINY( wsumM ) ) THEN
           CALL Fatal('CalculateMortarDistance','wsumM seems to be almost zero!')
         END IF

         ! Slave and master multipliers should sum up to same value
         mult = ABS( wsum / wsumM ) 
         
         ! Compute the weigted distance in the normal direction.
         DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
           k = ActiveProjector % Cols(j)

           IF( k > SIZE( FieldPerm ) ) THEN
             PRINT *,'k:',K,SIZE(FieldPerm)
             CALL Fatal('','Index too large')
           END IF
           
           l = FieldPerm( k ) 
           IF( l == 0 ) CYCLE

           ! This includes only the coordinate since the displacement
           ! is added to the coordinate!
           coeff = ActiveProjector % Values(j)

           CoeffSign = 1
           
           ! Only compute the sum related to the active projector
           IF( .NOT. SlaveNode(k) ) THEN
             coeff = mult * coeff
             IF( ThisRotatedContact ) CoeffSign = -1
           END IF
           
           disp(1) = DispVals( dofs * (l-1) + 1)
           disp(2) = DispVals( dofs * (l-1) + 2 )
           IF( cdofs == 2 ) THEN
             disp(3) = 0.0_dp
           ELSE
             disp(3) = DispVals( dofs * (l-1) + 3 )
           END IF

           ! If nonlinear analysis is used we may need to cancel the introduced gap due to numerical errors 
           IF( TieContact .AND. ResidualMode ) THEN !.AND. k <= dofs * Mesh % NumberOfNodes ) THEN
             IF( ThisRotatedContact ) THEN
               ContactVec(1) = ContactVec(1) + coeff * SUM( LocalNormal * Disp )
               ContactVec(2) = ContactVec(2) + coeff * SUM( LocalT1 * Disp )
               IF( cDofs == 3) ContactVec(3) = ContactVec(3) + coeff * SUM( LocalT2 * Disp )
             ELSE
               ContactVec(1) = ContactVec(1) + coeff * SUM( ContactNormal * Disp )
               ContactVec(2) = ContactVec(2) + coeff * SUM( ContactT1 * Disp )
               IF( cDofs == 3 ) ContactVec(3) = ContactVec(3) + coeff * SUM( ContactT2 * Disp ) 
             END IF
             CYCLE
           END IF

           coord(1) = Mesh % Nodes % x( k ) 
           coord(2) = Mesh % Nodes % y( k ) 
           coord(3) = Mesh % Nodes % z( k ) 

           PrevDisp = 0._dp
           IF( CalculateVelocity ) THEN
             PrevDisp(1) = PrevDispVals( dofs * (l-1) + 1)
             PrevDisp(2) = PrevDispVals( dofs * (l-1) + 2 )
             IF( cdofs == 2 ) THEN
               PrevDisp(3) = 0.0_dp
             ELSE
               PrevDisp(3) = PrevDispVals( dofs * (l-1) + 3 )
             END IF
           END IF

           ! If the linear system is in residual mode also set the current coordinate in residual mode too!
           ! Note that displacement field is given always in cartesian coordinates!
           IF( ResidualMode ) THEN
             Coord = Coord + Disp
           END IF

           ! DistN is used to give the distance that we need to move the original coordinates
           ! in the wanted direction in order to have contact.
           IF( ThisRotatedContact ) THEN
             ContactVec(1) = ContactVec(1) + coeff * SUM( LocalNormal * Coord )
           ELSE
             ContactVec(1) = ContactVec(1) + coeff * SUM( ContactNormal * Coord )
           END IF

           ! Tangential distances needed to move the original coordinates to the contact position
           ! If stick is required then we want to keep the tangential slip zero. 
           IF( StickContact ) THEN             
             SlipCoord = -PrevDisp 
             IF( ResidualMode ) SlipCoord = SlipCoord + Disp 

             IF( ThisRotatedContact ) THEN
               ContactVec(2) = ContactVec(2) + coeff * SUM( LocalT1 * SlipCoord )
               IF( cDofs == 3) ContactVec(3) = ContactVec(3) + coeff * SUM( LocalT2 * SlipCoord )
             ELSE
               ContactVec(2) = ContactVec(2) + coeff * SUM( ContactT1 * SlipCoord )
               IF( cDofs == 3 ) ContactVec(3) = ContactVec(3) + coeff * SUM( ContactT2 * SlipCoord )
             END IF
           END IF

           ! If not in the residual mode still take into account the displacement for the condition
           IF( .NOT. ResidualMode ) Coord = Coord + Disp

           ! Dist is used to compute the current signed distance that is used to determine
           ! whether we have contact or not. 
           IF( RotationalProjector ) THEN
             Dist = Dist + coeff * SQRT( SUM( Coord**2 ) )
           ELSE IF( NormalProjector ) THEN
             Dist = Dist + coeff * SUM( LocalNormal * Coord )
           ELSE             
             Dist = Dist + coeff * SUM( ContactNormal * Coord )
           END IF

           CartVec = CartVec + coeff * Coord
           
           IF( CalculateVelocity ) THEN
             Velo = ( Disp - PrevDisp ) !/ dt
             ContactVelo(1) = ContactVelo(1) + coeff * SUM( Velo * LocalNormal ) 
             ContactVelo(2) = ContactVelo(2) + coeff * SUM( Velo * LocalT1 )
             ContactVelo(3) = ContactVelo(3) + coeff * SUM( Velo * LocalT2 ) 
           END IF
           DistanceSet = .TRUE.
         END DO

         ! Divide by weight to get back to real distance in the direction of the normal
         ContactVec = ContactVec / wsum 
         Dist = DistSign * Dist / wsum
         IF( CalculateVelocity ) THEN
           ContactVelo = ContactVelo / wsum
         END IF
         CartVec = CartVec / wsum
         
200      IF( IsSlave ) THEN

           MortarBC % Rhs(cDofs*(i-1)+DofN) = -ContactVec(1)
           IF( StickContact .OR. TieContact ) THEN
             MortarBC % Rhs(cDofs*(i-1)+DofT1) = -ContactVec(2) 
             IF( cDofs == 3 ) THEN
               MortarBC % Rhs(cDofs*(i-1)+DofT2) = -ContactVec(3)
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

         GapVar % Values( j ) = ContactVec(1)

         IF( CalculateVelocity ) THEN
           DO k=1,cDofs             
             VeloVar % Values( cDofs*(j-1)+k ) = ContactVelo(k) 
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

       
       IF( LinearContactGap .OR. pContact ) THEN       
         DO elem=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
           
           Element => Mesh % Elements( elem )         
           
           IsSlave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) 
           IsMaster = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) 
           IF( .NOT. ( IsSlave .OR. ( CreateDual .AND. IsMaster ) ) ) CYCLE
           
           ElemCode = Element % TYPE % ElementCode           
           IF( pContact ) THEN
             n = mGetElementDOFs(pIndexes,Element)                   
             Indexes => pIndexes
           ELSE            
             n = Element % TYPE % NumberOfNodes
             Indexes => Element % NodeIndexes         
           END IF
           
           SELECT CASE ( ElemCode )
           CASE( 202, 203 )
             i=3
             
             j = DistVar % Perm(Indexes(i))
             IF( j > 0 ) THEN
               IF( pContact ) THEN
                 DistVar % Values(j) = 0.0_dp
                 GapVar % Values(j) = 0.0_dp
                 IF( CalculateVelocity ) THEN
                   DO k=1,cDofs
                     VeloVar % Values( cDofs*(j-1)+k ) = 0.0_dp
                   END DO
                 END IF
               ELSE
                 i1=1
                 i2=2
                 j1 = DistVar % Perm(Indexes(i1))
                 j2 = DistVar % Perm(Indexes(i2))
                 
                 DistVar % Values(j) = 0.5_dp * &
                     ( DistVar % Values(j1) + DistVar % Values(j2))
                 GapVar % Values(j) = 0.5_dp * &
                     ( GapVar % Values(j1) + GapVar % Values(j2))
                 
                 IF( CalculateVelocity ) THEN
                   DO k=1,cDofs
                     VeloVar % Values( cDofs*(j-1)+k ) = 0.5_dp * &
                         ( VeloVar % Values(cDofs*(j1-1)+k) + VeloVar % Values(cDofs*(j2-1)+k))  
                   END DO
                 END IF
               END IF
             END IF

           CASE( 404, 408 )
             DO i=5,8
               
               j = DistVar % Perm(Indexes(i))
               IF( j == 0 ) CYCLE
               
               IF( pContact ) THEN
                 DistVar % Values(j) = 0.0_dp
                 GapVar % Values(j) = 0.0_dp
                 IF( CalculateVelocity ) THEN
                   DO k=1,cDofs
                     VeloVar % Values( cDofs*(j-1)+k ) = 0.0_dp
                   END DO
                 END IF
               ELSE
                 i1=i-4
                 i2=i1+1
                 IF(i2==5) i2=1
                 j1 = DistVar % Perm(Indexes(i1))
                 j2 = DistVar % Perm(Indexes(i2))
                 
                 DistVar % Values(j) = 0.5_dp * &
                     ( DistVar % Values(j1) + DistVar % Values(j2))
                 GapVar % Values(j) = 0.5_dp * &
                     ( GapVar % Values(j1) + GapVar % Values(j2))
                 
                 IF( CalculateVelocity ) THEN
                   DO k=1,cDofs
                     VeloVar % Values( cDofs*(j-1)+k ) = 0.5_dp * &
                         ( VeloVar % Values(cDofs*(j1-1)+k) + VeloVar % Values(cDofs*(j2-1)+k))  
                   END DO
                 END IF
               END IF
             END DO
               
           CASE DEFAULT
             CALL Fatal('CalculateMortarDistance','Implement linear gaps for: '//I2S(ElemCode))
           END SELECT

         END DO
       END IF
       
       DEALLOCATE( SlaveNode )
       IF( CreateDual ) DEALLOCATE( MasterNode )
       
       IF( InfoActive(30) ) THEN
         ! We don't know if other partitions are here, so let us not make parallel reductions!
         CALL VectorValuesRange(DistVar % Values,SIZE(DistVar % Values),'Dist',.TRUE.)
         CALL VectorValuesRange(GapVar % Values,SIZE(GapVar % Values),'Gap',.TRUE.)
         CALL VectorValuesRange(MortarBC % rhs,SIZE(MortarBC % rhs),'Mortar Rhs',.TRUE.)       
       END IF
         
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
       LOGICAL :: LinearContactLoads
       INTEGER :: i1,i2,j1,j2,ElemCode,m
       
       n = Mesh % MaxElementNodes
       ALLOCATE(Basis(2*n), Nodes % x(2*n), Nodes % y(2*n), Nodes % z(2*n) )
       Nodes % x = 0.0_dp; Nodes % y = 0.0_dp; Nodes % z = 0.0_dp

       CALL Info(Caller,'Computing pressure for contact problems',Level=20)
       
       CoordSys = CurrentCoordinateSystem()
       NodalForce = 0.0_dp

       NormalSign0 = 0
       NormalCount = 0
       
       ALLOCATE( NodeDone( SIZE( FieldPerm ) ) )
       NodeDone = .FALSE.

       LinearContactLoads = ListGetLogical( Model % Simulation,&
           'Contact BCs linear loads', Found )

       
100    DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         
         Element => Mesh % Elements( elem )         

         IsSlave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) 
         IsMaster = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) 

         IF( .NOT. ( IsSlave .OR. ( CreateDual .AND. IsMaster ) ) ) CYCLE
                  
         IF( pContact ) THEN
           n = mGetElementDOFs(pIndexes,Element)                   
           Indexes => pIndexes
         ELSE         
           n = Element % TYPE % NumberOfNodes         
           Indexes => Element % NodeIndexes
         END IF
         
         Nodes % x(1:n) = Mesh % Nodes % x(Indexes(1:n))
         Nodes % y(1:n) = Mesh % Nodes % y(Indexes(1:n))
         Nodes % z(1:n) = Mesh % Nodes % z(Indexes(1:n))
         
         IntegStuff = GaussPoints( Element )

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
           IF( IsSlave .AND. ( FlatProjector .OR. PlaneProjector .OR. NormalProjector ) ) THEN
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
             IF( j == 0 ) CYCLE             
             
             IF( .NOT. NodeDone( Indexes(i) ) ) THEN             
               NodeDone( Indexes(i) ) = .TRUE.
               WeightVar % Values(j) = 0.0_dp
               NormalLoadVar % Values(j) = 0.0_dp
               SlipLoadVar % Values(j) = 0.0_dp
             END IF

             k = FieldPerm( Indexes(i) )
             IF( k == 0 ) CYCLE
             
             DO l=1,cdofs
               NodalForce(l) = LoadValues(dofs*(k-1)+l)
             END DO

             IF( ThisRotatedContact ) THEN
               NormalForce = NodalForce(1)
             ELSE
               NormalForce = SUM( NodalForce * Normal ) 
             END IF
             ! By construction the expression should be positive but due to rounding errors we have to ensure it!
             SlipForce = SQRT( MAX(0.0_dp, SUM( NodalForce**2 ) - NormalForce**2 ) )
             
             NormalLoadVar % Values(j) = NormalLoadVar % Values(j) - &
                 s * Basis(i) * NormalForce
             SlipLoadVar % Values(j) = SlipLoadVar % Values(j) + &
                 s * Basis(i) * SlipForce
             
             WeightVar % Values(j) = WeightVar % Values(j) + s * Basis(i)
           END DO
           
         END DO
       END DO
       
       ! Normalize the computed normal loads such that the unit will be that of pressure
       DO i=1,SIZE(FieldPerm)
         IF( NodeDone( i ) ) THEN             
           j = WeightVar % Perm(i)
           IF(j==0) CYCLE
           s = WeightVar % Values(j)
           IF( s /= s ) CYCLE
           IF( ABS(s) > EPSILON(s) ) THEN
             SlipLoadVar % Values(j) = SlipLoadVar % Values(j) / s**2
             NormalLoadVar % Values(j) = NormalLoadVar % Values(j) / s**2
           END IF
         END IF
       END DO

       IF( LinearContactLoads ) THEN       
         DO elem=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
           
           Element => Mesh % Elements( elem )         

           IsSlave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) 
           IsMaster = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) 
           IF( .NOT. ( IsSlave .OR. ( CreateDual .AND. IsMaster ) ) ) CYCLE
           
           Indexes => Element % NodeIndexes         
           n = Element % TYPE % NumberOfNodes
           ElemCode = Element % TYPE % ElementCode           

           SELECT CASE ( ElemCode )
             
           CASE( 203, 306, 408 )
             n = MODULO( ElemCode, 100 )
             m = ElemCode / 100             

             DO i=m+1,n
               i1=i-m
               IF(i1==m) THEN                 
                 i2=m+1
               ELSE
                 i2=1
               END IF
               j = SlipLoadVar % Perm(Indexes(i))
               j1 = SlipLoadVar % Perm(Indexes(i1))
               j2 = SlipLoadVar % Perm(Indexes(i2))
               SlipLoadVar % Values(j) = 0.5_dp * &
                   ( SlipLoadVar % Values(j1) + SlipLoadVar % Values(j2))
               NormalLoadVar % Values(j) = 0.5_dp * &
                   ( NormalLoadVar % Values(j1) + NormalLoadVar % Values(j2))
             END DO
               
           CASE DEFAULT
             CALL Fatal(Caller,'Implement linear loads for: '//I2S(ElemCode))
           END SELECT
         END DO
       END IF

       IF( InfoActive(20) ) THEN
         CALL VariableValuesRange(SlipLoadVar,'SlipLoad',AlwaysSerial=.TRUE.)       
         CALL VariableValuesRange(NormalLoadVar,'NormalLoad',AlwaysSerial=.TRUE.)       
       END IF
       
       IF( FlatProjector .OR. PlaneProjector .OR. NormalProjector ) THEN
         IF( NormalCount == 0 ) THEN
           CALL Info(Caller,'All normals are consistently signed',Level=10)
         ELSE
           CALL Warn(Caller,'There are normals with conflicting signs: '&
               //I2S(NormalCount)  )
           NormalSign = 1
         END IF
         CALL Info(Caller,'Normal direction for distance measure: '&
             //I2S(NormalSign),Level=8)
         DistSign = NormalSign 
       END IF

       ! Check whether the normal sign has been enforced
       IF( ListGetLogical( BC,'Normal Sign Negative',Found ) ) DistSign = -1
       IF( ListGetLogical( BC,'Normal Sign Positive',Found ) ) DistSign = 1

       DEALLOCATE( Basis, Nodes % x, Nodes % y, Nodes % z, NodeDone )

       CALL Info(Caller,'Finished computing contact pressure',Level=30)
       
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

       IF( InfoActive(30) ) THEN      
         PRINT *,'InitialActiveSet:',COUNT( MortarBC % Active ), SIZE( MortarBC % Active )
       END IF
       
       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = UseLoadVar % Perm(j)

         ind = cDofs * (i-1) + DofN

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
           NodeLoad = UseLoadVar % Values(k)
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

       IF( InfoActive(20) ) THEN
         IF ( -HUGE(MaxDist) /= MaxDist ) THEN
           IF( MaxDist - MinDist >= 0.0_dp ) THEN
             PRINT *,'NormalContactSet Dist:',MinDist,MaxDist
           END IF
         END IF
         IF ( -HUGE(MaxLoad) /= MaxLoad) THEN
           IF( MaxLoad - MinLoad >= 0.0_dp ) THEN
             PRINT *,'NormalContactSet Load:',MinLoad,MaxLoad
           END IF
         END IF
         PRINT *,'NormalContactSet active count:',COUNT(MortarBC % Active)
         PRINT *,'NormalContactSet passive count:',COUNT(.NOT. MortarBC % Active)
       END IF

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' nodes to the set'
         CALL Info(Caller,Message,Level=6)
       END IF
       
       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' nodes from the set'
         CALL Info(Caller,Message,Level=6)
       END IF

       IF( InfoActive(30) ) THEN      
         PRINT *,'ModifiedActiveSet:',COUNT( MortarBC % Active ), SIZE( MortarBC % Active )
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

       LimitedNow = COUNT( MortarBC % active(DofN::cDofs) )      
       NewNodes = LimitedMin - LimitedNow
       IF( NewNodes <= 0 ) RETURN

       WRITE(Message,'(A,I0)') 'Initial number of contact nodes for '&
           //TRIM(VarName)//': ',LimitedNow 
       CALL Info(Caller,Message,Level=5)

       CALL Info(Caller,&
           'Setting '//I2S(NewNodes)//' additional contact nodes',Level=5)

       ALLOCATE( DistArray( NewNodes ), IndArray( NewNodes ) ) 
       DistArray = HUGE( DistArray ) 
       IndArray = 0

       ! Find additional contact nodes from the closest non-contact nodes
       DO i = 1,Projector % NumberOfRows
         ind = cDofs * (i-1) + DofN
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
         CALL Fatal(Caller,'Could not define sufficient number of new nodes!')
       END IF

       WRITE(Message,'(A,ES12.4)') 'Maximum distance needed for new nodes:',DistArray(NewNodes)
       CALL Info(Caller,Message,Level=8)

       MortarBC % Active( cDofs*(IndArray-1)+DofN ) = .TRUE.

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

       
       CALL Info(Caller,'Setting Tangent contact set',Level=20)
       
       IF( FrictionContact .AND. &
           ListGetLogical( BC,'Stick Contact Global',Found ) ) THEN
        
         ! Sum up global normal and slide forces
         DO i = 1,Projector % NumberOfRows
           j = Projector % InvPerm( i ) 
           IF( j == 0 ) CYCLE
           k = FieldPerm( j ) 
           IF( k == 0 ) CYCLE
           k = UseLoadVar % Perm(j)
                      
           ! If there is no contact there can be no stick either
           indN = cDofs * (i-1) + DofN
           IF( .NOT. MortarBC % Active(indN) ) CYCLE

           NodeLoad = UseLoadVar % Values(k)
           TangentLoad = SlipLoadVar % Values(k)
         
           mustatic = ListGetRealAtNode( BC,'Static Friction Coefficient', j )
           mudynamic = ListGetRealAtNode( BC,'Dynamic Friction Coefficient', j )
           IF( ( mustatic - mudynamic ) < -EPSILON(mustatic) ) THEN
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
         MortarBC % Active( DofT1 :: cDofs ) = .FALSE.
         IF( cDofs == 3 ) THEN
            MortarBC % Active( DofT2 :: cDofs ) = .FALSE.
          END IF
          GOTO 100 
       ELSE IF( StickContact .OR. TieContact ) THEN
         MortarBC % Active( DofT1 :: cDofs ) = MortarBC % Active( DofN :: cDofs )
         IF( cDofs == 3 ) THEN
           MortarBC % Active( DofT2 :: cDofs ) = MortarBC % Active( DofN :: cDofs ) 
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
         k = UseLoadVar % Perm(j)

         indN = cDofs * (i-1) + DofN
         indT1 = ind - DofN + DofT1
         IF(cDofs == 3 ) indT2 = ind - DofN + DofT2

         ! If there is no contact there can be no stick either
         IF( .NOT. MortarBC % Active(indN) ) THEN
           IF( MortarBC % Active(indT1) ) THEN
             removed0 = removed0 + 1
             MortarBC % Active(indT1) = .FALSE.
             IF( cDofs == 3 ) MortarBC % Active(indT2) = .FALSE.
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
           IF( cDofs == 3 ) MortarBC % Active(indT2) = .TRUE.
           CYCLE
         END IF

         ! Enforce no-stick condition (=slip)
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Stick Passive Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           IF( MortarBC % Active(IndT1) ) removed = removed + 1
           MortarBC % Active(indT1) = .FALSE.
           IF( cDofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           CYCLE
         END IF

         ! Remove nodes with too large tangent force
         !--------------------------------------------------------------------------       

         NodeLoad = UseLoadVar % Values(k)
         TangentLoad = SlipLoadVar % Values(k)
         
         mustatic = ListGetRealAtNode( BC,'Static Friction Coefficient', j )
         mudynamic = ListGetRealAtNode( BC,'Dynamic Friction Coefficient', j )

         IF( ( mustatic - mudynamic ) < -EPSILON(mustatic) ) THEN
           CALL Warn('TangentContactSet','Static friction coefficient should be larger than dynamic!')
         END IF

         IF( MortarBC % Active(IndT1) ) THEN
           IF( TangentLoad > mustatic * ABS( NodeLoad ) ) THEN
             removed = removed + 1
             MortarBC % Active(indT1) = .FALSE.
             IF( cDofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           END IF
         ELSE              
           stickcoeff = ListGetRealAtNode( BC,'Stick Contact Coefficient', j, Found )
           IF( Found ) THEN
             DO l=1,cDofs             
               du(l) = VeloVar % Values( cDofs*(k-1)+l ) 
             END DO
             IF( cDofs == 3 ) THEN
               Slip = SQRT(du(dofT1)**2 + du(DofT2)**2)
             ELSE
               Slip = ABS( du(dofT1) )
             END IF
             IF( stickcoeff * slip  < mudynamic * ABS( NodeLoad ) ) THEN
               added = added + 1
               MortarBC % Active(indT1) = .TRUE.
               IF( cDofs == 3 ) MortarBC % Active(indT2) = .TRUE.
             END IF
           END IF
         END IF
       END DO

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' nodes to the stick set'
         CALL Info(Caller,Message,Level=6)
       END IF
       
       IF(removed0 > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed0,' non-contact nodes from the stick set'
         CALL Info(Caller,Message,Level=6)
       END IF

       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' sliding nodes from the stick set'
         CALL Info(Caller,Message,Level=6)
       END IF


100    CALL Info(Caller,'Creating fields out of normal and stick contact sets',Level=10)

       DO i = 1, Projector % NumberOfRows
         j = Projector % InvPerm(i)
         IF( j == 0 ) CYCLE
         k = NormalActiveVar % Perm(j)

         IF( MortarBC % Active(cDofs*(i-1)+DofN) ) THEN
           NormalActiveVar % Values(k) = 1.0_dp
         ELSE
           NormalActiveVar % Values(k) = -1.0_dp
         END IF

         IF( MortarBC % Active(cDofs*(i-1)+DofT1) ) THEN
           StickActiveVar % Values(k) = 1.0_dp
         ELSE
           StickActiveVar % Values(k) = -1.0_dp
         END IF
       END DO

       IF( InfoActive(30) ) THEN
         PRINT *,'Active Tangent set:',COUNT( MortarBC % Active ) 
         CALL VariableValuesRange(NormalActiveVar,'NormalActive',AlwaysSerial=.TRUE.)
         CALL VariableValuesRange(StickActiveVar,'StickActive',AlwaysSerial=.TRUE.)
       END IF
         
     END SUBROUTINE TangentContactSet



     ! Sets the diagonal entry for slip in the tangent direction(s).
     ! This coefficient may be used to relax the stick condition, and also to
     ! revert back nodes from slip to stick set. 
     !----------------------------------------------------------------------------------
     SUBROUTINE StickCoefficientSet() 
       
       REAL(KIND=dp) :: NodeLoad, TangentLoad
       INTEGER :: i,j,k,ind,IndN, IndT1, IndT2
       LOGICAL :: Found

       CALL Info('StickCoefficientSet','Setting the stick coefficient entry for tangent components at stick',Level=10)

       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = UseLoadVar % Perm(j)

         indN = cDofs * (i-1) + DofN
         indT1 = cDofs * (i-1) + DofT1
         IF(cDofs == 3 ) indT2 = cDofs * (i-1) + DofT2

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
         IF( cDofs == 3 ) MortarBC % Diag(indT2) = coeff
       END DO

       IF(InfoActive(30)) THEN
         CALL VectorValuesRange( MortarBC % Diag, SIZE( MortarBC % Diag),'MortarBC Diag')
       END IF
       
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
             ElemInds(i) = cDofs * ( j - 1) + DofN
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
           CALL Fatal(Caller,'Cannot deal with element: '//I2S(elemcode))

         END SELECT
       END DO

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' quadratic nodes to contact set'
         CALL Info(Caller,Message,Level=6)
       END IF

       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' quadratic nodes from contact set'
         CALL Info(Caller,Message,Level=6)
       END IF
         
     END SUBROUTINE QuadraticContactSet


     ! Project contact fields from slave to master
     !----------------------------------------------------------------------------------------
     SUBROUTINE ProjectFromSlaveToMaster()

       REAL(KIND=dp) :: Disp(3), Coord(3), PrevDisp(3), Velo(3), ContactVelo(3), &
           LocalNormal0(3), SlipCoord(3)
       REAL(KIND=dp), POINTER :: DispVals(:), PrevDispVals(:) 
       REAL(KIND=dp) :: MinDist, MaxDist, CoeffEps
       LOGICAL, ALLOCATABLE :: SlaveNode(:), NodeDone(:)
       REAL(KIND=dp), ALLOCATABLE :: CoeffTable(:), RealActive(:)
       INTEGER :: i,j,k,l,l2,Indexes(12)

       CALL Info(Caller,'Mapping entities from slave to master',Level=10)

       n = SIZE( FieldPerm )
       ALLOCATE( SlaveNode( n ) )
       SlaveNode = .FALSE.
       
       DO i=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( i )                  
         IF( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE

         CurrentModel % CurrentElement => Element
         
         IF( pContact ) THEN
           n = mGetElementDOFs(pIndexes,Element)
           SlaveNode(pIndexes(1:n)) = .TRUE.
         ELSE
           SlaveNode( Element % NodeIndexes ) = .TRUE.
         END IF
       END DO

       IF( InfoActive(20) ) THEN
         n = COUNT( SlaveNode )
         CALL Info(Caller,'Number of dofs on slave side: '//I2S(n))
       END IF
         
       n = SIZE( DistVar % Values )
       ALLOCATE( CoeffTable( n ), NodeDone( n ) )
           
       CoeffTable = 0.0_dp
       NodeDone = .FALSE.
       

       DO i = 1,Projector % NumberOfRows             
         
         IF( Projector % InvPerm(i) == 0 ) CYCLE          
         l = DistVar % Perm( Projector % InvPerm(i) )
         
         IF(.NOT. pContact ) THEN
           IF( l > Mesh % NumberOfNodes ) CYCLE
         END IF

         DO j = Projector % Rows(i),Projector % Rows(i+1)-1
           k = Projector % Cols(j)

           IF(.NOT. pContact ) THEN
             IF( k > Mesh % NumberOfNodes ) CYCLE
           END IF                  
           
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
               DO k=1,cDofs             
                 VeloVar % Values( cDofs*(l2-1)+k ) = 0.0_dp
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
             DO k=1,cDofs             
               VeloVar % Values( cDofs*(l2-1)+k ) = VeloVar % Values( cDofs*(l2-1)+k ) + &
                   coeff * VeloVar % Values( cDofs*(l-1)+k)
             END DO
           END IF
         END DO
       END DO
       
       CoeffEps = 1.0d-8 * MAXVAL( ABS( CoeffTable ) )
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
               DO k=1,cDofs
                 VeloVar % Values( cDofs*(i-1)+k ) = VeloVar % Values( cDofs*(i-1)+k ) / CoeffTable( i ) 
               END DO
             END IF
           ELSE
             NormalLoadVar % Values( i ) = 0.0_dp
             SlipLoadVar % Values( i ) = 0.0_dp
             IF( CalculateVelocity ) THEN
               DO k=1,cDofs
                 VeloVar % Values( cDofs*(i-1)+k ) = 0.0_dp
               END DO
             END IF             
           END IF

         END IF
       END DO

       DO i = 1, Projector % NumberOfRows
         j = Projector % InvPerm(i)
         IF( j == 0 ) CYCLE
         
         IF( .NOT. pContact ) THEN
           IF( j > Mesh % NumberOfNodes ) CYCLE
         END IF
         
         k = NormalActiveVar % Perm(j)
         
         IF( NormalActiveVar % Values( k ) < 0.0_dp ) THEN
           IF( CalculateVelocity ) THEN
             DO l=1,cDofs
               VeloVar % Values( cDofs*(k-1)+l ) = 0.0_dp
             END DO
           END IF
         END IF
       END DO
       
       IF( InfoActive(30) ) THEN
         CALL Info('PojectFromSlaveToMaster','Projecting fields')
         CALL VariableValuesRange(NormalLoadVar,'NormalLoadVar',AlwaysSerial=.TRUE.)
         CALL VariableValuesRange(SlipLoadVar,'SlipLoadVar',AlwaysSerial=.TRUE.)
         IF( CalculateVelocity ) THEN
           CALL VariableValuesRange(VeloVar,'VeloVar',AlwaysSerial=.TRUE.)
         END IF
       END IF

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
       CHARACTER(LEN=MAX_NAME_LEN) :: ContactVeloName


       IF(.NOT. ListCheckPresent( BC, 'Dynamic Friction Coefficient') ) RETURN
      
       CALL Info(Caller,'Setting contact friction for boundary',Level=10)

       ContactVeloName = 'Contact Velocity'
       GivenDirection = ListCheckPresent( BC, ContactVeloName )        
       IF( .NOT. GivenDirection .AND. TimeStep == 1 ) THEN
         ContactVeloName = 'Initial Contact Velocity'
         GivenDirection = ListCheckPresent( BC, ContactVeloName ) 
       END IF
       
       IF(.NOT. GivenDirection ) THEN
         IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
           CALL Info(Caller,'Contact velocity not defined: Give "Contact Velocity" or',Level=3)
           CALL Info(Caller,'"Initial Contact Velocity" or use e.g. "Initial Stick Model = True"')
           CALL Fatal(Caller,'Cannot continue with friction without direction!')
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
             IF( cDofs == 3 ) LocalT2 = NTT(:,3)
           ELSE
             Rotated = .FALSE.
             LocalNormal = ContactNormal
             LocalT1 = ContactT1
             IF( cDofs == 3 ) LocalT2 = ContactT2 
           END IF
           
           VeloCoeff = 0.0_dp           
           VeloSign = 1

           IF( GivenDirection ) THEN
             IF( Slave ) THEN
               VeloDir => ListGetConstRealArray( BC, ContactVeloName )
             ELSE
               VeloDir => ListGetConstRealArray( MasterBC, ContactVeloName, Found)
               IF(.NOT. Found ) THEN
                 ! If velocity direction not found in master then use the opposite velocity of the slave
                 VeloDir => ListGetConstRealArray( BC, ContactVeloName )
                 VeloSign = -1
               END IF
             END IF
             VeloCoeff(DofT1) = SUM( VeloDir(1:3,1) * LocalT1 )
             IF( cDofs == 3 ) THEN
               VeloCoeff(DofT2) = SUM( VeloDir(1:3,1) * LocalT2 )
             END IF
           ELSE
             VeloCoeff(DofT1) = VeloVar % Values(cDofs*(k-1)+DofT1) 
             IF(cDofs==3) VeloCoeff(DofT2) = VeloVar % Values(cDofs*(k-1)+DofT2) 
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
           k = cDOFs * (j-1) + DofN 

           k2 = cDOFs * (j-1) + DofT1 
           A % Rhs(k2) = A % Rhs(k2) - VeloCoeff(DofT1) * A % Rhs(k)

           IF( cDofs == 3 ) THEN
             k3 = cDOFs * (j-1) + DofT2
             A % Rhs(k3) = A % Rhs(k3) - VeloCoeff(DofT2) * A % Rhs(k)             
           END IF

           DO l = A % Rows(k),A % Rows(k+1)-1
             DO l2 = A % Rows(k2), A % Rows(k2+1)-1
               IF( A % Cols(l2) == A % Cols(l) ) EXIT
             END DO

             A % Values(l2) = A % Values(l2) - VeloCoeff(DofT1) * A % Values(l)
             
             IF( cDofs == 3 ) THEN
               DO l3 = A % Rows(k3), A % Rows(k3+1)-1
                 IF( A % Cols(l3) == A % Cols(l) ) EXIT
               END DO
               A % Values(l3) = A % Values(l3) - VeloCoeff(DofT2) * A % Values(l)
             END IF
           END DO
         END DO
       END DO
       
       n = COUNT( NodeDone ) 
       CALL Info(Caller,'Number of friction nodes: '//I2S(n),Level=10)
       
       DEALLOCATE( NodeDone )       

       IF( InfoActive(30) ) THEN
         CALL VectorValuesRange(A % Values,SIZE(A % Values),'A-friction')
       END IF
       
     END SUBROUTINE SetSlideFriction
     
   END SUBROUTINE DetermineContact



!> Sets one Dirichlet condition to the desired value
!------------------------------------------------------------------------------
   SUBROUTINE UpdateDirichletDof( A, dof, dval )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    INTEGER :: dof
    REAL(KIND=dp) :: dval

    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF
    
    IF(.NOT. ALLOCATED(A % Dvalues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF
    
    A % Dvalues( dof ) = dval
    A % ConstrainedDOF( dof ) = .TRUE.
    
  END SUBROUTINE UpdateDirichletDof
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE UpdateDirichletDofC( A, dof, cval )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    INTEGER :: dof
    COMPLEX(KIND=dp) :: cval

    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF
    
    IF(.NOT. ALLOCATED(A % Dvalues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF
    
    A % Dvalues( 2*dof-1 ) = REAL( cval )
    A % ConstrainedDOF( 2*dof-1 ) = .TRUE.

    A % Dvalues( 2*dof ) = AIMAG( cval )
    A % ConstrainedDOF( 2*dof ) = .TRUE.
    
  END SUBROUTINE UpdateDirichletDofC
!------------------------------------------------------------------------------



  
!> Releases one Dirichlet condition 
!------------------------------------------------------------------------------
   SUBROUTINE ReleaseDirichletDof( A, dof )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    INTEGER :: dof
    REAL(KIND=dp) :: dval
      
    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF
    
    IF(.NOT.ALLOCATED(A % Dvalues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF
    
    A % ConstrainedDOF( dof ) = .FALSE.
    
  END SUBROUTINE ReleaseDirichletDof
!------------------------------------------------------------------------------


  
!> Release the range or min/max values of Dirichlet values.
!------------------------------------------------------------------------------
  FUNCTION DirichletDofsRange( Solver, Oper ) RESULT ( val ) 
!------------------------------------------------------------------------------
    TYPE(Solver_t), OPTIONAL :: Solver
    CHARACTER(LEN=*), OPTIONAL :: Oper 
    REAL(KIND=dp) :: val
    
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: minv,maxv
    LOGICAL :: FindMin, FindMax
    INTEGER :: i,OperNo
    
    IF( PRESENT( Solver ) ) THEN
      A => Solver % Matrix
    ELSE
      A => CurrentModel % Solver % Matrix
    END IF
    
    val = 0.0_dp
    
    ! Defaulting to range
    OperNo = 0

    IF( PRESENT( Oper ) ) THEN
      IF( Oper == 'range' ) THEN
        OperNo = 0
      ELSE IF( Oper == 'min' ) THEN
        OperNo = 1 
      ELSE IF( Oper == 'max' ) THEN
        OperNo = 2
      ELSE
        CALL Fatal('DirichletDofsRange','Unknown operator: '//TRIM(Oper))
      END IF
    END IF
          
    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      RETURN
    END IF
  
    IF( OperNo == 0 .OR. OperNo == 1 ) THEN
      minv = HUGE( minv ) 
      DO i=1,SIZE( A % ConstrainedDOF )
        IF( A % ConstrainedDOF(i) ) minv = MIN( A % DValues(i), minv ) 
      END DO
      minv = ParallelReduction( minv, 1 ) 
    END IF

    IF( OperNo == 0 .OR. OperNo == 2 ) THEN
      maxv = -HUGE( maxv ) 
      DO i=1,SIZE( A % ConstrainedDOF )
        IF( A % ConstrainedDOF(i) ) maxv = MAX( A % DValues(i), maxv ) 
      END DO
      maxv = ParallelReduction( maxv, 2 ) 
    END IF
    
    IF( OperNo == 0 ) THEN    
      val = maxv - minv
    ELSE IF( OperNo == 1 ) THEN
      val = minv
    ELSE
      val = maxv
    END IF
      
  END FUNCTION DirichletDofsRange
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Show range of vector values mainly for debugging purposes.
!------------------------------------------------------------------------------
  SUBROUTINE VectorValuesRange(x,n,str,AlwaysSerial)
!------------------------------------------------------------------------------    
    REAL(KIND=dp) :: x(:)
    INTEGER :: n
    CHARACTER(LEN=*) :: str
    LOGICAL, OPTIONAL :: AlwaysSerial
    
    INTEGER :: np
    REAL(KIND=dp) :: s(3)
    LOGICAL :: Parallel, Found

    Parallel = ( ParEnv % PEs > 1)
    IF( Parallel ) Parallel = .NOT. CurrentModel % Mesh % SingleMesh 
    IF( Parallel ) THEN
      IF( PRESENT( AlwaysSerial ) ) THEN
        IF( AlwaysSerial ) Parallel = .FALSE.
      END IF
    END IF
    IF( Parallel ) THEN
      IF( ListGetLogical( CurrentModel % Simulation,'Serial Range',Found ) ) THEN
        Parallel = .FALSE.
      END IF
    END IF

    np = n
      
    IF( np == 0 ) THEN
      s(1) = HUGE(s(1))
      s(2) = -HUGE(s(2))
      s(3) = 0.0_dp        
    ELSE          
      s(1) = MINVAL( x(1:np) ) 
      s(2) = MAXVAL( x(1:np) ) 
      s(3) = SUM( x(1:np) ) 
    END IF
      
    IF( Parallel ) THEN
      np = ParallelReduction( np )
      IF( np > 0 ) THEN
        s(1) = ParallelReduction( s(1),1 ) 
        s(2) = ParallelReduction( s(2),2 ) 
        s(3) = ParallelReduction( s(3) )
      END IF
    END IF

    IF( np == 0 ) THEN
      WRITE(Message,*) 'Size of vector is zero: '//TRIM(str)
      CALL Info('VectorValuesRange',Message)              
    ELSE
      !s(3) = s(3) / np
      !WRITE(Message,*) '[min,max,ave] for '//TRIM(str)//':', s
      WRITE(Message,*) '[size,min,max,sum] for '//TRIM(str)//': '//I2S(np)//' ', s
      CALL Info('VectorValuesRange',Message)
    END IF
        
  END SUBROUTINE VectorValuesRange
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Show range of field variable mainly for debugging purposes.
!------------------------------------------------------------------------------
  SUBROUTINE VariableValuesRange(pVar,str,AlwaysSerial)
    TYPE(Variable_t), POINTER :: pVar
    CHARACTER(LEN=*) :: str
    LOGICAL, OPTIONAL :: AlwaysSerial
    
    REAL(KIND=dp), POINTER :: pVals(:)
    INTEGER :: n
    
    IF( ASSOCIATED(pVar) ) THEN
      pVals => pVar % Values
      n = SIZE( pVals )
    ELSE
      pVals => NULL()
      n = 0
    END IF
    CALL VectorValuesRange(pVals,n,str,AlwaysSerial)
    
  END SUBROUTINE VariableValuesRange


  
  SUBROUTINE FindClosestNode(Mesh,Coord,MinDist,MinNode,Parallel,Eps,Perm)
    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: Coord(:)
    REAL(KIND=dp) :: MinDist
    INTEGER :: MinNode
    LOGICAL :: Parallel
    REAL(KIND=dp), OPTIONAL :: Eps
    INTEGER, OPTIONAL :: Perm(:)
    
    INTEGER :: i,j,k,n,dim
    REAL(KIND=dp) :: Dist
    INTEGER, POINTER :: neigh(:)
    
    CALL Info('FindClosestNode','Trying to find the closest node',Level=20)
    
    dim = SIZE( Coord )
    MinDist = HUGE(MinDist)
    MinNode = 0
    
    DO i=1,Mesh % NumberOfNodes
      IF( PRESENT( Perm ) ) THEN
        IF( Perm(i) == 0) CYCLE
      END IF
      Dist = (Mesh % Nodes % x(i) - Coord(1))**2 
      IF(dim >= 2) Dist = Dist + (Mesh % Nodes % y(i) - Coord(2))**2
      IF(dim == 3) Dist = Dist + (Mesh % Nodes % z(i) - Coord(3))**2
      
      IF( PRESENT( Eps ) ) THEN
        IF( Dist > Eps ) CYCLE
      END IF
      
      IF(Dist < MinDist ) THEN
        MinDist = Dist
        MinNode = i
      END IF
    END DO
    
    ! We found minimum square
    MinDist = SQRT(MinDist)    
    
    !PRINT *,'MinNode serial:',ParEnv % MyPe, MinDist, MinNode, &
    !    Mesh % Nodes % x(MinNode), Mesh % Nodes % y(MinNode)
    
    ! In parallel case eliminate all except the nearest node. 
    ! This relies on the fact that for each node partition the 
    ! distance to nearest node is computed accurately. 
    IF( Parallel ) THEN

      ! In parallel apply load only on the owner partition:
      ! On shared nodes penalize the non-owner candidate. 
      ! ---------------------------------------------------
      IF(MinNode>0) THEN
        neigh => Mesh % ParallelInfo % NeighbourList(MinNode) % Neighbours
        ! Find the 1st index active for this model
        n = SIZE(neigh)
        IF(n>1) THEN         
          ! Of the sharing partitions who owns this
          DO i=1,n
            IF(ParEnv % Active(neigh(i)+1)) EXIT
          END DO
          ! And if it is not in this partition then skip it
          IF(i<n) THEN
            IF(neigh(i) /= ParEnv % MyPE) THEN          
              MinDist = HUGE(MinDist)
            END IF
          END IF
        END IF
      END IF

      ! These should be legit candidates.
      ! Let's choose the one with shortest distance
      Dist = MinDist 
      MinDist = ParallelReduction( Dist, 1 )
      IF( ABS( MinDist - Dist ) > EPSILON(Dist) ) THEN
        MinNode = 0
      END IF          
    END IF

    !PRINT *,'MinNode:',MinNode, ParEnv % MyPe
    
    CALL Info('FindClosestNode','Closest node found to be: '//I2S(MinNode),Level=20)
    
    WRITE(Message,'(A,ES12.5)') 'Closest node distance: ',MinDist
    CALL Info('FindClosestNode',Message,Level=20)
    
  END SUBROUTINE FindClosestNode


  SUBROUTINE TargetCoordinatesToTargetNodes( Mesh, ValueList, Success, Perm )
    TYPE(Mesh_t) :: Mesh
    TYPE(ValueList_t), POINTER :: ValueList
    LOGICAL :: Success
    INTEGER, OPTIONAL :: Perm(:)
    
    LOGICAL :: Found, Parallel  
    REAL(KIND=dp), POINTER :: CoordNodes(:,:)
    REAL(KIND=dp) :: Eps, MinDist
    INTEGER :: i,j,NoNodes, NoDims, NofNodesFound
    INTEGER, POINTER :: IndNodes(:)
    
    
    Success = .FALSE.

    Parallel = ( ParEnv % PEs > 1 )
    IF(Mesh % SingleMesh) Parallel = .FALSE.
    
    CoordNodes => ListGetConstRealArray(ValueList,'Target Coordinates',Found)
    IF(.NOT. Found) RETURN

    Eps = ListGetConstReal( ValueList, 'Target Coordinates Eps',Found )
    IF ( .NOT. Found ) THEN
      Eps = HUGE(Eps)
    ELSE
      ! We are looking at square of distance
      Eps = Eps**2
    END IF
    
    NoNodes = SIZE(CoordNodes,1)
    NoDims = SIZE(CoordNodes,2)    
    IF(NoNodes == 0 ) RETURN
    
    ALLOCATE( IndNodes(NoNodes) )
    IndNodes = -1
            
    NofNodesFound = 0
    DO j=1,NoNodes           
      CALL FindClosestNode(Mesh,CoordNodes(j,:),MinDist,i,&
          Parallel,Eps,Perm)
      IF(i>0) THEN
        NofNodesFound = NofNodesFound + 1
        IndNodes(NofNodesFound) = i
      END IF
    END DO
    
    ! If no nodes found, add still an empty list and make sure the 
    ! zero is not treated later on. Otherwise this search would be 
    ! retreated each time. 
    NofNodesFound = MAX(1,NofNodesFound)
    
    CALL ListAddIntegerArray( ValueList,'Target Nodes', &
        NOFNodesFound, IndNodes) 
    
    ! Finally deallocate the temporal vectors
    DEALLOCATE( IndNodes )
    Success = .TRUE.
    
  END SUBROUTINE TargetCoordinatesToTargetNodes

  
!------------------------------------------------------------------------------
!> Set Dirichlet boundary condition for given dof. The conditions are
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
    INTEGER :: BC,i,j,j2,k,l,m,n,nd,p,t,k1,k2,OffSet
    LOGICAL :: GotIt, periodic, OrderByBCNumbering, ReorderBCs
    REAL(KIND=dp), POINTER :: WorkA(:,:,:) => NULL()
    REAL(KIND=dp) ::  s

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver

    LOGICAL :: Conditional
    LOGICAL, ALLOCATABLE :: DonePeriodic(:)
    CHARACTER(:), ALLOCATABLE :: CondName, DirName, PassName, PassCondName, EqName

    INTEGER :: NoNodes,NoDims,bf_id,nlen, NOFNodesFound, dim, &
        bndry_start, bndry_end, Upper
    REAL(KIND=dp), POINTER :: CoordNodes(:,:), Condition(:), Work(:)
    REAL(KIND=dp) :: MinDist,Dist, Eps
    LOGICAL, ALLOCATABLE :: ActivePart(:), ActiveCond(:), ActivePartAll(:)
    TYPE(ValueList_t), POINTER :: ValueList, Params
    LOGICAL :: NodesFound, Passive, OffDiagonal, ApplyLimiter
    LOGICAL, POINTER :: LimitActive(:)
    TYPE(Variable_t), POINTER :: Var

    TYPE(Element_t), POINTER :: Parent

    INTEGER :: ind, ElemFirst, ElemLast, bf, BCstr, BCend, BCinc, dgind
    REAL(KIND=dp) :: SingleVal
    LOGICAL :: AnySingleBC, AnySingleBF
    LOGICAL, ALLOCATABLE :: LumpedNodeSet(:)
    LOGICAL :: NeedListMatrix
    INTEGER :: DirCount
    CHARACTER(*), PARAMETER :: Caller = 'SetDirichletBoundaries'
    LOGICAL, ALLOCATABLE :: CandNodes(:)
    INTEGER, ALLOCATABLE :: LumpedIndx(:)
    INTEGER, POINTER :: PlaneInds(:)
    LOGICAL :: Parallel

    TYPE(NormalTangential_t), POINTER :: NT
    LOGICAL, ALLOCATABLE, SAVE :: NTzeroing_done(:,:)
    INTEGER, ALLOCATABLE, SAVE :: NTelement(:,:)

    REAL(KIND=dp) :: Mult(Model % MaxElementNodes), MaxMult, ParMaxMult, MoveCoeff
    LOGICAL :: GotMult
    INTEGER :: maxind
    
!------------------------------------------------------------------------------
! These logical vectors are used to minimize extra effort in setting up different BCs
!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)
    n = MAX( Model % NumberOfBodyForces,Model % NumberOfBCs)
    IF( n == 0 ) THEN
      CALL Info(Caller,'No BCs or Body Forces present, exiting early...',Level=12)
    END IF

    ALLOCATE( ActivePart(n), ActivePartAll(n), ActiveCond(n))
    CondName = Name(1:nlen) // ' Condition'
    PassName = Name(1:nlen) // ' Passive'
    PassCondName = Name(1:nlen) // ' Condition' // ' Passive'

    OffSet = 0
    OffDiagonal = .FALSE.
    IF( PRESENT( PermOffSet) ) OffSet = PermOffSet
    IF( PRESENT( OffDiagonalMatrix ) ) OffDiagonal = OffDiagonalMatrix

    Mesh => Model % Mesh
    ALLOCATE( Indexes(Mesh % MaxElementDOFs) )

    NT => Model % Solver % NormalTangential
    n = NT % NormalTangentialNOFNodes 
    IF( n > 0 ) THEN
      ! We need to have these available for different components of the same vector.
      ! Hence a dirty compromise between localility and saving values. 
      m = 0
      IF( ALLOCATED( NTElement ) ) THEN
        m = SIZE( NTElement, 1 )
      END IF
      IF( m /= n .AND. m > 0 ) THEN
        DEALLOCATE( NTzeroing_done, NTelement )
      END IF
      IF( m /= n ) THEN      
        ALLOCATE( NTzeroing_done(n,3), NTelement(n,3) )
        NTZeroing_done = .FALSE.
        NTelement = 0
      END IF
    END IF
    
    Parallel = ( ParEnv % PEs > 1 ) .AND. ( .NOT. Mesh % SingleMesh ) 
    
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
       CALL Fatal(Caller,'Periodicity not considered with offset')
     END IF

     ALLOCATE( DonePeriodic( Mesh % NumberOFNodes ) )
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

     IF( InfoActive(12) ) THEN
       CALL Info(Caller,'Number of periodic points set: '&
           //I2S(COUNT(DonePeriodic)),Level=12)
     END IF

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
          CALL Warn(Caller,"Requested 'Dirichlet BC Order' but &
               &not 'Set Dirichlet BCs by BC Numbering', ignoring...")
       ELSE IF(SIZE(BCOrder) /= Model % NumberOfBCs) THEN
          CALL Fatal(Caller,"'Dirichlet BC Order' is the wrong length!")
       END IF
    END IF

    bndry_start = Mesh % NumberOfBulkElements+1
    bndry_end   = bndry_start + Mesh % NumberOfBoundaryElements-1
    DirCount = 0

    ! check and set some flags for nodes belonging to n-t boundaries
    ! getting set by other bcs:
    ! --------------------------------------------------------------
    IF ( NT % NormalTangentialNOFNodes>0 ) THEN
      IF ( OrderByBCNumbering ) THEN
        DO i=1,Model % NumberOfBCs
          BC = i
          IF(ReorderBCs) BC = BCOrder(BC)
          IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
          Conditional = ActiveCond(BC)

          DO t = bndry_start, bndry_end
            Element => Mesh % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                   Model % BCs(BC) % Tag ) CYCLE

            ValueList => Model % BCs(BC) % Values
            Model % CurrentElement => Element

            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p) == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = mGetElementDOFs( Indexes )
            END IF
            CALL CheckNTelement(n,t)
          END DO
        END DO
      ELSE
        DO t = bndry_start, bndry_end
          DO BC=1,Model % NumberOfBCs
            IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
            Conditional = ActiveCond(BC)
          
            Element => Mesh % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                 Model % BCs(BC) % Tag ) CYCLE
          
            ValueList => Model % BCs(BC) % Values
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p) == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = mGetElementDOFs( Indexes )
            END IF
            CALL CheckNTelement(n,t)
          END DO
        END DO
      END IF

      IF ( DOF<= 0 ) THEN
        DO t=bndry_start,bndry_end
          Element => Mesh % Elements(t)
          n = Element % TYPE % NumberOfNodes
          DO j=1,n
            k = NT % BoundaryReorder(Element % NodeIndexes(j))
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
        DO i=1,Model % NumberOfBCs
          BC = i
          IF(ReorderBCs) BC = BCOrder(BC)
          IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
          Conditional = ActiveCond(BC)

          DO t = bndry_start, bndry_end
            Element => Mesh % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                Model % BCs(BC) % Tag ) CYCLE
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p) == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = mGetElementDOFs( Indexes )
            END IF
            ValueList => Model % BCs(BC) % Values
            CALL SetElementValues(n)
          END DO
        END DO
      ELSE
        DO t = bndry_start, bndry_end
          DO BC=1,Model % NumberOfBCs
            IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
            Conditional = ActiveCond(BC)
            
            Element => Mesh % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BC) % Tag ) CYCLE
            
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p)  == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = mGetElementDOFs( Indexes )
            END IF
            ValueList => Model % BCs(BC) % Values
            CALL SetElementValues(n)
          END DO
        END DO
      END IF
    END IF

    BLOCK
      INTEGER,ALLOCATABLE :: ChildBCs(:)
      INTEGER,POINTER :: BCInds(:)
      REAL(KIND=dp) :: zero = 0.0_dp
      IF( ListGetLogical( Model % Solver % Values,'Extruded Child BC Zero',GotIt ) ) THEN
        CALL Info(Caller,'Setting extruded BCs (start & end) to zero!',Level=10)

        ALLOCATE(ChildBCs(2*Model % NumberOfBodies))
        ChildBCs=0

        ! Collect a list of child BCs that were generated when extruding
        m = 0
        DO i=1,Model % NumberOfBodies
          BCInds => ListGetIntegerArray(Model % Bodies(i) % Values,'Extruded Child BCs',GotIt)
          IF(GotIt) THEN
            ChildBCs(m+1:m+2) = BCInds(1:2)
            m=m+2
          END IF          
        END DO        
        IF(m==0) CALL Fatal(Caller,'No "Extruded Child BCs" to set')

        IF( InfoActive(20) ) THEN
          PRINT *,'ChildBCs:',m,ChildBCs(1:m)
        END IF
          
        ! Set the extruded BCs to zero. Note that only this value is available currently.
        DO t = bndry_start, bndry_end
          Element => Mesh % Elements(t)
          IF(ANY(ChildBCs(1:m) == Element % BoundaryInfo % Constraint ) ) THEN          
            Model % CurrentElement => Element
            n = mGetElementDOFs( Indexes, Element, Model % Solver )            
            DO i=1,n
              CALL SetSinglePoint(Indexes(i),DOF,zero,.TRUE.)
            END DO
          END IF
        END DO
      END IF
    END BLOCK
      
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
      EqName = ListGetString( Solver % Values, 'Equation', GotIt )

      IF( Passive ) THEN
        ALLOCATE(PassPerm(Mesh % NumberOfNodes),NodeIndexes(1))
        PassPerm = 0
        DO i=0,Mesh % PassBCCnt-1
          j=Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements-i
          PassPerm(Mesh % Elements(j) % NodeIndexes)=1
        END DO
      END IF
        
      DO t=1,Solver % Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        IF( Element % BodyId <= 0 .OR. Element % BodyId > Model % NumberOfBodies ) THEN
          CALL Warn(Caller,'Element body id beyond body table!')
          CYCLE
        END IF
                    
        bf_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force', GotIt)
        
        IF(.NOT. GotIt) CYCLE
        IF(.NOT. ActivePart(bf_id) .AND. .NOT. ActivePartAll(bf_id)) CYCLE
        Conditional = ActiveCond(bf_id)

        Model % CurrentElement => Element

        GotIt = CheckElementEquation( Model, Element, EqName ) 
        
        n = Element % TYPE % NumberOfNodes
        Indexes(1:n) = Element % NodeIndexes

        ValueList => Model % BodyForces(bf_id) % Values
        IF(.NOT. ASSOCIATED( ValueList ) ) CYCLE
        
        IF (ListGetLogical(ValueList,PassCondName,GotIt)) THEN
          IF (.NOT.CheckPassiveElement(Element)) CYCLE

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

            NodeIndexes(1) = Indexes(j)
            IF(PassPerm(NodeIndexes(1))==0) CALL SetPointValues(1)
          END DO
        ELSE
          CALL SetElementValues(n)
        END IF
        
        ! Set the higher p-dofs related to Dirichlet BC's to zero. 
        IF ( isActivePElement(Element, Solver) ) THEN
          nd = mGetElementDOFs( Indexes, Uelement = Element, USolver = Model % Solver )
          DO i=n+1,nd            
            CALL SetSinglePoint(Indexes(i),DOF,0.0_dp,.TRUE.)            
          END DO
        END IF
        
      END DO
      
      IF(Passive) DEALLOCATE(PassPerm,NodeIndexes)
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
      ! At the first calling the list of coordinates is transformed to list of nodes.
      IF(.NOT. NodesFound) THEN
        IF( ListCheckPresent( ValueList,'Target Coordinates' ) ) THEN
          CALL TargetCoordinatesToTargetNodes( Mesh, ValueList, NodesFound )
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

    ! Check the boundaries and body forces for possible single nodes BCs that are used to fixed 
    ! the domain for undetermined equations. The loop is slower than optimal in the case that there is 
    ! a large amount of different boundaries that have a node to set. 
    !--------------------------------------------------------------------------------------------
    DirName = TRIM(Name)//' Single Node'
    AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
    AnySingleBF = ListCheckPresentAnyBodyForce( Model, DirName )

    IF( AnySingleBC .OR. AnySingleBF ) THEN
      Solver => Model % Solver
      Mesh   => Solver % Mesh

      DO bc = 1,Model % NumberOfBCs  + Model % NumberOfBodyForces    

        ! Make a distinction between BCs and BFs. 
        ! These are treated in the same loop because most of the logic is still the same. 
        IF( bc <= Model % NumberOfBCs ) THEN
          IF(.NOT. AnySingleBC ) CYCLE
          ValueList => Model % BCs(BC) % Values
          ElemFirst =  Mesh % NumberOfBulkElements + 1 
          ElemLast =  Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        ELSE
          IF( .NOT. AnySingleBF ) CYCLE
          ValueList => Model % BodyForces(bc - Model % NumberOfBCs) % Values
          ElemFirst =  1
          ElemLast =  Mesh % NumberOfBulkElements
        END IF

        SingleVal = ListGetCReal( ValueList,DirName, GotIt) 
        IF( .NOT. GotIt ) CYCLE
        ind = ListGetInteger( ValueList,TRIM(Name)//' Single Node Index',GotIt )     
        
        ! On the first time find a one single uniquely defined node for setting 
        ! the value. In parallel it will be an unshared node with the highest possible 
        ! node number 
        IF(.NOT. GotIt ) THEN                  
          ind = 0
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)
            n = Element % TYPE % NumberOfNodes
            NodeIndexes => Element % NodeIndexes

            IF( bc <= Model % NumberOfBCs ) THEN
              IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
            ELSE
              j = Element % BodyId
              IF( j < 0 .OR. j > Model % NumberOfBodies ) CYCLE
              bf = ListGetInteger( Model % Bodies(j) % Values,'Body Force',GotIt)
              IF(.NOT. GotIt) CYCLE
              IF( bc - Model % NumberOfBCs /= bf ) CYCLE
            END IF
            
            DO i=1,n
              j = NodeIndexes(i)
              IF( Perm(j) == 0) CYCLE
              IF( Parallel ) THEN
                IF( SIZE( Mesh % ParallelInfo % NeighbourList(j) % Neighbours) > 1 ) CYCLE               
                IF( Mesh % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPe ) CYCLE               
              END IF
              ind = j 
              EXIT
            END DO
            IF( ind > 0 ) EXIT
          END DO

          k = ind
          IF( Parallel ) THEN
            k = ParallelReduction( ind, 2 ) 
            
            ! Find the maximum partition that owns a suitable node. 
            ! It could be minimum also, just some convection is needed. 
            k = -1
            IF( ind > 0 ) k = ParEnv % MyPe          
            k = ParallelReduction( k, 2 ) 
            IF( k == -1 ) THEN
              CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
            ELSE
              IF( k /= ParEnv % MyPe ) ind = 0                         
              IF( InfoActive(8) ) THEN
                ind = ParallelReduction( ind, 2 )                
                CALL Info(Caller,'Fixing single node '&
                    //I2S(ind)//' at partition '//I2S(k),Level=8)
                IF( k /= ParEnv % MyPe ) ind = 0
              END IF
            END IF
          ELSE
            IF( ind == 0 ) THEN
              CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
            ELSE              
              CALL Info(Caller,'Fixing single node '//I2S(ind),Level=8)
            END IF
          END IF
            
          CALL ListAddInteger( ValueList,TRIM(Name)//' Single Node Index', ind )          
        END IF

        ! Ok, if this is the partition where the single node to eliminate the floating should 
        ! be eliminated then set it here. Index equal to zero tells that we are in a wrong partition.        
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

      ALLOCATE(PassPerm(Mesh % NumberOfNodes),NodeIndexes(1));PassPerm=0
      DO i=0,Mesh % PassBCCnt-1
        j=Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements-i
        PassPerm(Mesh % Elements(j) % NodeIndexes)=1
      END DO

      DO i=1,Solver % NumberOfActiveElements
        Element => Mesh % Elements(Solver % ActiveElements(i))
        IF (CheckPassiveElement(Element)) THEN
          nd = mGetElementDOFs(Indexes,UElement=Element)
          n  = Element % Type % NumberOfNodes
          DO j=1,nd
            k = Indexes(j)
            IF (k<=0) CYCLE

            k=Perm(k)
            IF (k<=0) CYCLE

            IF(Indexes(j)<=SIZE(PassPerm)) THEN
              IF(PassPerm(Indexes(j))==1) CYCLE
            END IF

!           s=0._dp
!           DO l=1,NDOFs
!             m=NDOFs*(k-1)+l
!             s=s+ABS(A % Values(A % Diag(m)))
!           END DO
!           IF (s>EPSILON(s)) CYCLE


            DO l=1,NDOFs
              s=0._dp
              m=NDOFs*(k-1)+l
              s=s+ABS(A % Values(A % Diag(m)))
              IF (s>EPSILON(s)) CYCLE

              m = NDOFs*(k-1)+l
              IF(A % ConstrainedDOF(m)) CYCLE
              CALL SetSinglePoint(k,l,Solver % Variable % Values(m),.FALSE.)
            END DO
          END DO
        END IF
      END DO
      DEALLOCATE(PassPerm,NodeIndexes)
    END IF

    ! Check the boundaries and body forces for possible single nodes BCs that must have a constant
    ! value on that boundary / body force.
    !--------------------------------------------------------------------------------------------
    NeedListMatrix = .FALSE.
    
    DO l = 0, 1
      IF( l == 0 ) THEN
        DirName = TRIM(Name)//' Constant'
      ELSE
        DirName = TRIM(Name)//' Profile'
      END IF
      
      AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
      AnySingleBF = ListCheckPresentAnyBodyForce( Model, DirName )

      IF( AnySingleBC .OR. AnySingleBF ) THEN
        GotMult = (l == 1 )
        ALLOCATE( LumpedNodeSet( SIZE( Perm ) ) )

        IF( AnySingleBC ) CALL Info(Caller,'Found BC constraint: '//TRIM(DirName),Level=6)
        IF( AnySingleBF ) CALL Info(Caller,'Found BodyForce constraint: '//TRIM(DirName),Level=6)

        ! Improve the logic in future
        ! Now we assume that if the "supernode" has been found then also the matrix has the correct topology. 
        IF( AnySingleBC ) THEN
          NeedListMatrix = .NOT. ListCheckPresentAnyBC( Model, TRIM(Name)//' Constant Node Index')
        ELSE 
          NeedListMatrix = .NOT. ListCheckPresentAnyBodyForce( Model, TRIM(Name)//' Constant Node Index')
        END IF

        ! Move the list matrix because of its flexibility. Register the initial topology.
        IF( NeedListMatrix ) THEN
          CALL CRS_ChangeTopology(A, Init=.TRUE.)
          CALL List_toListMatrix(A)
        END IF
      
        DO bc = 1,Model % NumberOfBCs + Model % NumberOfBodyForces
          ! Make a distinction between BCs and BFs. 
          ! These are treated in the same loop because most of the logic is still the same. 
          IF( bc <= Model % NumberOfBCs ) THEN
            IF(.NOT. AnySingleBC ) CYCLE
            ValueList => Model % BCs(BC) % Values
            ElemFirst =  Mesh % NumberOfBulkElements + 1 
            ElemLast =  Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          ELSE
            IF(.NOT. AnySingleBF ) CYCLE
            ValueList => Model % BodyForces(bc - Model % NumberOfBCs) % Values
            ElemFirst =  1
            ElemLast =  Model % NumberOfBulkElements
          END IF

          IF(GotMult) THEN
            IF(.NOT. ListCheckPresent( ValueList,DirName) ) CYCLE
          ELSE            
            IF(.NOT. ListGetLogical( ValueList,DirName, GotIt) ) CYCLE
          END IF
                    
          MoveCoeff = ListGetCReal( ValueList,TRIM(DirName)//' Resistance',GotIt)
          IF(.NOT. GotIt) MoveCoeff = 1.0_dp

          ! This tells us that this has been visited before            
          ind = ListGetInteger( ValueList,TRIM(DirName)//' Node Index',GotIt )               
         
          ! On the first time find a one single uniquely defined node for setting 
          ! the value. In parallel it will be an unshared node with the highest possible 
          ! node number 
          IF( GotIt ) THEN        
            dgind = ListGetInteger( ValueList,TRIM(DirName)//' DG Node Index',GotIt )               
            IF( GotMult ) THEN
              MaxMult = ListGetConstReal( ValueList,TRIM(DirName)//' Max Value',UnfoundFatal=.TRUE.) 
            END IF
          ELSE
            MaxMult = 0.0_dp          
            ind = 0
            dgind = 0
            maxind = 0
            
            DO t = ElemFirst, ElemLast
              Element => Mesh % Elements(t)

              IF( bc <= Model % NumberOfBCs ) THEN
                IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
              ELSE
                bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values,'Body Force',GotIt)
                IF( bc - Model % NumberOfBCs /= bf ) CYCLE
              END IF

              n = Element % TYPE % NumberOfNodes
              NodeIndexes => Element % NodeIndexes
                
              IF(GotMult) THEN
                Mult(1:n) = ListGetReal( ValueList,TRIM(DirName),n,NodeIndexes,UnfoundFatal=.TRUE.)
              END IF

              DO i=1,n
                j = NodeIndexes(i)

                IF(Model % Solver % DG) THEN                                   
                  CALL PickDgIndexes(Element,Indexes)                  
                  IF( Perm(Indexes(i)) == 0) CYCLE
                ELSE
                  IF( Perm(j) == 0) CYCLE
                END IF
                  
                IF( Parallel ) THEN
                  IF( SIZE( Mesh % ParallelInfo % NeighbourList(j) % Neighbours) > 1 ) CYCLE               
                  IF( Mesh % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPe ) CYCLE               
                END IF
                IF( GotMult) THEN
                  ! Find the point with maximum value of the multiplier
                  IF( ABS(Mult(i)) > ABS(MaxMult) ) THEN
                    maxind = j
                    MaxMult = Mult(i)
                  END IF
                ELSE
                  ind = j 
                  IF( Model % Solver % DG ) dgind = Indexes(i)
                  EXIT
                END IF
              END DO
              IF( ind > 0 ) EXIT
            END DO

            IF(GotMult) ind = maxind

            ! Find the maximum partition that owns the node. 
            ! It could be minimum also, just some convention is needed. 
            IF( Parallel ) THEN
              IF( GotMult ) THEN
                ParMaxMult = ABS(MaxMult)
                ParMaxMult = ParallelReduction( ParMaxMult, 2 ) 
                IF(ABS(ABS(MaxMult)-ParMaxMult) > 1.0e-3*ParMaxMult) ind = 0 
              END IF
              k = -1
              IF( ind > 0 ) k = ParEnv % MyPe          
              k = ParallelReduction( k, 2 ) 
              IF( k == -1 ) THEN
                CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
              END IF
              j = 0
              IF( k /= ParEnv % MyPe .AND. ind > 0) THEN
                ind = 0
                dgind = 0
                j = 1 
              END IF
              ! Just a counter for partitions that have hits but do not own the index.
              j = ParallelReduction(j)
              IF(j>0) THEN
                CALL Warn(Caller,'Problems: boundary extends over '//I2S(j+1)//' partitions')
              END IF
            ELSE
              IF( ind == 0 ) THEN
                CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
              END IF
            END IF

            CALL ListAddInteger( ValueList,TRIM(DirName)//' Node Index', ind )
            IF( Model % Solver % DG ) THEN
              CALL ListAddInteger( ValueList,TRIM(DirName)//' DG Node Index', dgind )
            END IF

            IF(GotMult) THEN
              CALL ListAddConstReal( ValueList,TRIM(DirName)//' Max Value', MaxMult )
            END IF
            WRITE(Message,'(A,ES12.3)') TRIM(DirName)//' Max Value', MaxMult
            CALL Info(Caller,Message,Level=15)

            NeedListMatrix = .TRUE.
          END IF

          ! This is probably in parallel a passive partition. 
          IF( ind == 0 ) CYCLE

          ! Ok, now sum up the rows to the corresponding nodal index
          ! We go through elements and need to mark the nodes in order not to visit them twice. 
          LumpedNodeSet = .FALSE.

          ! Actually from this on we do not need the node index if we have a DG field.
          ! Replace it with the DG index.
          IF( Model % Solver % DG ) THEN
            IF( dgind == 0 ) THEN
              CALL Fatal(Caller,'Could not determine DG index associated to node index!')
            END IF
            ind = dgind
          END IF

          ! If we need to do scaling then play with the supernode too!
          IF( Ndofs == 1 ) THEN
            BLOCK
              INTEGER :: k0
              k0 = Offset + NDOFs * (Perm(ind)-1) + DOF
              IF(ABS(MoveCoeff-1.0) > EPSILON(MoveCoeff)) THEN
                CALL MoveRow( A, k0, k0, MoveCoeff )
                b(k0) = MoveCoeff * b(k0)
              END IF
            END BLOCK
          END IF
          ! Supernode has been set, if needed. 
          LumpedNodeSet(ind) = .TRUE.
                    
          IF(.NOT. ALLOCATED(LumpedIndx)) THEN
            ALLOCATE(LumpedIndx(Model % NumberOfBCs + Model % NumberOfBodyForces))
            LumpedIndx = 0
          END IF
          LumpedIndx(bc) = ind
            
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)

            IF( bc <= Model % NumberOfBCs ) THEN
              IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
            ELSE
              bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values,'Body Force',GotIt)
              IF( bc - Model % NumberOfBCs /= bf ) CYCLE
            END IF

            n = Element % TYPE % NumberOfNodes
            IF ( Model % Solver % DG ) THEN              
              CALL PickDgIndexes(Element,Indexes)                  
            ELSE
              Indexes(1:n) = Element % NodeIndexes
            END IF
                
            IF(GotMult) Mult(1:n) = ListGetReal( ValueList,TRIM(Name)//' Profile',n,&
                Element % NodeIndexes,UnfoundFatal=.TRUE.)
                       
            CALL SetLumpedRows(ind,n)
          END DO

          IF( NDOFs == 1 ) THEN
            SingleVal = ListGetCReal( ValueList,TRIM(DirName)//' Flux',GotIt)
            IF( GotIt ) THEN
              t = Offset + Perm(ind)
              b(t) = b(t) + SingleVal
            END IF

            SingleVal = ListGetCReal( ValueList,TRIM(DirName)//' Coefficient',GotIt)
            IF( GotIt ) THEN                        
              t = Offset + Perm(ind)
              CALL AddToMatrixElement(A,t,t,SingleVal) 
            END IF                        
          END IF

          n = COUNT( LumpedNodeSet ) 
          CALL Info(Caller,'Number of lumped nodes set: '//I2S(n),Level=10)
        END DO
      END IF
    END DO
      
    IF( NeedListMatrix ) THEN
      DEALLOCATE( LumpedNodeSet )
      
      ! Revert back to CRS matrix and change to the new topology. 
      CALL List_ToCRSMatrix(A)
      CALL CRS_ChangeTopology(A, Init=.FALSE.)

      CALL Info(Caller,'Modified matrix non-zeros: '&
          //I2S(SIZE( A % Cols )),Level=8)
    END IF
      
    ! We are back to CRS matrix.
    ! If we have a fixed point iteration we may add the flux multiplied by resistance to rhs as well. 
    BLOCK
      INTEGER :: k0
      REAL(KIND=dp) :: prevFlux
      TYPE(ValueList_t), POINTER :: vList
      IF(ALLOCATED(LumpedIndx) ) THEN
        DO bc = 1,SIZE(LumpedIndx) 
          ind = LumpedIndx(bc)
          IF(ind>0) THEN
            k0 = Offset + NDOFs * (Perm(ind)-1) + DOF
            prevFlux = CRS_MatrixRowVectorMultiply(A,Model % Solver % Variable % Values,k0)            

            ! In case this was added, remove it from the flux.
            SingleVal = ListGetCReal( ValueList,TRIM(Name)//' Constant Coefficient',GotIt)
            prevFlux = prevFlux - SingleVal * Model % Solver % Variable % Values(k0)
            
            IF( bc <= Model % NumberOfBCs ) THEN
              Vlist => Model % BCs(bc) % Values
              WRITE(Message,'(A,ES12.3)') 'Previous bc '//I2S(bc)//' lumped flux: ',prevFlux 
            ELSE
              i = bc-Model % NumberOfBCs
              Vlist => Model % BodyForces(i) % Values
              WRITE(Message,'(A,ES12.3)') 'Previous bf '//I2S(i)//' lumped flux: ',prevFlux 
            END IF            
            CALL Info(Caller,Message)
            DirName = TRIM(Name)//' Constant Prev Flux'
            CALL ListAddConstReal(Vlist, DirName, prevFlux )
          END IF
        END DO
      END IF
    END BLOCK

    ! Check the boundaries for a curve bc.
    !--------------------------------------------------------------------------------------------
    BLOCK
      INTEGER :: l,dofs,c1,c2,c3,ivec
      REAL(KIND=dp) :: Coeff,d(3),Dfdx(3),f,x(3),a11,a22
      LOGICAL, ALLOCATABLE :: NodeDone(:)
      TYPE(Variable_t), POINTER :: Dvar
      INTEGER, POINTER :: GivenComps(:)
      INTEGER :: Comps(3)
      LOGICAL :: AnyHingeBC, HingeBC
      REAL(KIND=dp) :: Normal(3),Tan1(3),Tan2(3),xt(2),cfit(7)
      
      
      DirName = TRIM(Name)//' Curve'
      AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
      AnyHingeBC = ListGetLogicalAnyBC( Model, TRIM(Name)//' Hinge' )
      
      IF( AnySingleBC .OR. AnyHingeBC ) THEN
        IF( AnySingleBC ) THEN
          CALL Info(Caller,'Found BC constraint for curve: '//TRIM(DirName),Level=6)
        END IF          
        IF( AnySingleBC ) THEN
          CALL Info(Caller,'Found BC constraint for hinge: '//TRIM(Name)//' Hinge',Level=6)
        END IF
          
        Solver => Model % Solver
                
        Dvar => Solver % Variable
        IF(.NOT. ASSOCIATED( DVar ) ) THEN
          CALL Fatal(Caller,'Solver variable not associated for Curve constraint!')
        END IF
        
        dim = CoordinateSystemDimension()        
        dofs = MIN( dim, DVar % Dofs )
        IF( dofs < 2 ) THEN
          CALL Fatal(Caller,'Curve constraint only makes sense for vector fields!')
        END IF

        ALLOCATE(NodeDone(Mesh % NumberOfNodes))
        NodeDone = .FALSE.
        
        DO bc = 1,Model % NumberOfBCs
          
          ValueList => Model % BCs(BC) % Values

          HingeBC = .FALSE.
          IF( AnyHingeBC ) THEN
            HingeBC = ListGetLogical( ValueList,TRIM(Name)//' Hinge',GotIt)
          END IF
          
          IF(HingeBC ) THEN            
            IF( dim == 2 ) THEN
              CALL CylinderFit(Mesh, ValueList, bc, dim, cfit )
            ELSE IF( dim == 3 ) THEN
              CALL CylinderFit(Mesh, ValueList, bc, dim, cfit )
              Normal = cfit(4:6)
              CALL TangentDirections(Normal,Tan1,Tan2)
            END IF
            !PRINT *,'Hinge Params',cfit
          ELSE           
            IF( .NOT. ListCheckPresent( ValueList,DirName) ) CYCLE
          END IF

          Comps = 0
          GivenComps => ListGetIntegerArray( ValueList, TRIM(Dirname)//' Components',GotIt)
          IF(GotIt) THEN
            dofs = MIN(dofs, SIZE(GivenComps) )
            Comps(1:dofs) = GivenComps(1:dofs)
          ELSE
            DO i=1,dofs
              Comps(i) = i
            END DO
          END IF
                    
          ElemFirst = Mesh % NumberOfBulkElements + 1 
          ElemLast = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
                         
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)

            IF(.NOT. ASSOCIATED(Element % BoundaryInfo) ) CYCLE
            IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
            n = Element % TYPE % NumberOfNodes
            Indexes(1:n) = Element % NodeIndexes
            
            DO i=1,n
              j = Indexes(i)
              IF(NodeDone(j)) CYCLE
              NodeDone(j) = .TRUE.
              k = DVar % Perm(j)

              IF(k==0) CYCLE

              l = DVar % dofs*(k-1)
              d(1:dofs) = DVar % Values(l+1:l+dofs)
              IF(dofs==2) d(3) = 0.0_dp
              
              x(1) = Mesh % Nodes % x(j) 
              x(2) = Mesh % Nodes % y(j)
              x(3) = Mesh % Nodes % z(j)

              x = x + d
              
              IF( HingeBC ) THEN
                ! We can analytically derive the case of 2d hinge
                IF( dim == 2 ) THEN
                  x(1:2) = x(1:2) - cfit(1:2)
                  f = x(1)**2 + x(2)**2 - cfit(3)**2
                  DfDx(1) = 2*x(1)
                  DfDx(2) = 2*x(2)
                ELSE
                  x = x - cfit(1:3)
                  xt(1) = SUM(x * Tan1)
                  xt(2) = SUM(x * Tan2)
                                    
                  f = xt(1)**2 + xt(2)**2 - cfit(7)**2
                  
                  DfDx(1) = 2*xt(1)*Tan1(1) + 2*xt(2)*Tan2(1)
                  DfDx(2) = 2*xt(1)*Tan1(2) + 2*xt(2)*Tan2(2)
                  DfDx(3) = 2*xt(1)*Tan1(3) + 2*xt(2)*Tan2(3)                  
                END IF
              ELSE                
                f = ListGetFunVec( ValueList, DirName, x(1:dofs), dofs, DfDx=dfdx(1:dofs) )
              END IF

              ! Check whether this is an normal-tangential node. If it is then m>0. 
              ! We already know that this has active perm for DVar. 
              m = 0
              IF ( NT % NormalTangentialNOFNodes > 0 ) THEN
                m = NT % BoundaryReorder(j) 
              END IF
                                         
              IF( m == 0 ) THEN                
                ! Let us take the most sensitive component to be the one for
                ! which the curve constraint is applied ensuring maximum diagonal entry.              ´              
                ! Then choose 2nd (and 3rd) components in order. 
                c1 = comps(1)
                DO ivec=2,dofs
                  IF( ABS(DfDx(comps(ivec))) > ABS(DfDx(c1))) c1 = comps(ivec)
                END DO
                DO ivec=1,dofs
                  IF(comps(ivec) /= c1) THEN
                    c2 = comps(ivec)
                    EXIT
                  END IF
                END DO

                IF( dofs == 3 ) THEN
                  c3 = 6 - c1 - c2
                END IF
               
                ! It may happen that the rows are linearly dependent on each other.
                ! Then a solution of type (x+y) for both is not good. Rather use then
                ! (x-y) for the other to have a unique solution. 
                a11 = GetMatrixElement(A,l+c1,l+c1)
                a22 = GetMatrixElement(A,l+c2,l+c2)
                
                ! This is a simple sign rule that avoids two equations being redundant.
                ! Don't multiply too numbers that could be almost zero!
                Coeff = -SIGN(1.0_dp,a11) * SIGN(1.0_dp,a22) * &
                    SIGN(1.0_dp, DfDx(c1)) * SIGN(1.0_dp, DfDx(c2))
              ELSE
                ! For normal-tangential system the normal component (1st one) should
                ! by construction be most sensitive to deviations from the curve. 
                c1 = 1
                c2 = 2
                c3 = 3
                Coeff = 1.0_dp
              END IF
                                                      
              ! Move all the entries from "c1" to "c2" and nullify the row.
              CALL MoveRow( A, l+c1, l+c2, Coeff )
              b(l+c2) = b(l+c2) + Coeff * b(l+c1)
              b(l+c1) = 0.0_dp

              ! In parallel only the owner makes the curve constraint
              IF( ParEnv % PEs > 1 ) THEN
                IF( A % ParallelInfo % NeighbourList(l+1) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
              END IF

              ! Residual mode is never active at this stage
              b(l+c1) = -f              
              b(l+c1) = b(l+c1) + DfDx(1)*d(1) + DfDx(2)*d(2)
              IF(dofs==3) b(l+c1) = b(l+c1) + DfDx(3)*d(3)
              
              IF(m>0) THEN                           
                CALL RotateNTSystem( DfDx, j )
              END IF

              CALL AddToMatrixElement( A, l+c1, l+1, DfDx(1) )
              CALL AddToMatrixElement( A, l+c1, l+2, DfDx(2) )
              IF(dofs==3) CALL AddToMatrixElement( A, l+c1, l+3, DfDx(3) )
            END DO
          END DO
        END DO
        
        n = COUNT( NodeDone ) 
        CALL Info(Caller,'Number of curved nodes set: '//I2S(n),Level=10)        
      END IF
      
    END BLOCK

    ! Check the boundaries and body forces for possible single nodes BCs that must have a constant
    ! value on that boundary / body force.
    !--------------------------------------------------------------------------------------------
    DirName = TRIM(Name)//' Plane'
    AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
    
    IF( AnySingleBC ) THEN
      dim = CoordinateSystemDimension()
      
      ALLOCATE( LumpedNodeSet( SIZE( Perm ) ) )

      CALL Info(Caller,'Found BC constraint: '//TRIM(DirName),Level=6)

      ! Improve the logic in future
      ! Now we assume that if the "supernode" has been found then also the matrix has the correct topology. 
      NeedListMatrix = .NOT. ListCheckPresentAnyBC( Model, TRIM(Name)//' Plane Node Indices')
      
      ! Move the list matrix because of its flexibility
      IF( NeedListMatrix ) THEN
        CALL CRS_ChangeTopology(A,Init=.TRUE.)
        CALL List_toListMatrix(A)
      END IF

      ElemFirst = Mesh % NumberOfBulkElements + 1 
      ElemLast = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      DO bc = 1,Model % NumberOfBCs 

        ValueList => Model % BCs(BC) % Values
        IF( .NOT. ListGetLogical( ValueList,DirName, GotIt) ) CYCLE

        PlaneInds => ListGetIntegerArray( ValueList,TRIM(Name)//' Plane Node Indices',GotIt )     

        IF(.NOT. GotIt ) THEN
          IF(.NOT. ALLOCATED(CandNodes) ) THEN
            ALLOCATE( CandNodes( Mesh % NumberOfNodes ) )        
          END IF
          CandNodes = .FALSE.

          ! Add nodes to the set that are associated with this BC only.
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)            
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
              NodeIndexes => Element % NodeIndexes
              CandNodes(NodeIndexes) = .TRUE.
            END IF
          END DO

          ! Remove nodes from the set that may be set by other BCs also. 
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)            
            IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) THEN
              NodeIndexes => Element % NodeIndexes
              CandNodes(NodeIndexes) = .FALSE.
            END IF
          END DO

          ALLOCATE(PlaneInds(3))
          CALL FindExtremumNodes(Mesh,CandNodes,dim,PlaneInds) 
          
          CALL ListAddIntegerArray( ValueList,TRIM(Name)//' Plane Node Indices',dim, PlaneInds )
          NeedListMatrix = .TRUE.
        END IF

        IF( Parallel ) CALL Warn(Caller,'Node index perhaps not set properly in parallel')
        ! IF( ind == 0 ) CYCLE

        ! Ok, now sum up the rows to the corresponding nodal index
        LumpedNodeSet = .FALSE.

        ! Don't lump the "supernodes" and therefore mark it set already
        LumpedNodeSet(PlaneInds) = .TRUE.

        DO t = ElemFirst, ElemLast
          Element => Mesh % Elements(t)

          IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
            n = Element % TYPE % NumberOfNodes
            Indexes(1:n) = Element % NodeIndexes
            CALL SetRigidRows(PlaneInds,bc,n)
          END IF
        END DO

        n = COUNT( LumpedNodeSet ) 
        CALL Info(Caller,'Number of lumped nodes set: '//I2S(n),Level=10)
      END DO

      IF( NeedListMatrix ) THEN
        DEALLOCATE( LumpedNodeSet )
        
        ! Revert back to CRS matrix
        CALL List_ToCRSMatrix(A)
        CALL CRS_ChangeTopology(A,Init=.FALSE.)
      END IF
    END IF

    
    IF( InfoActive(12) )  THEN
      IF( Parallel ) THEN
        DirCount = ParallelReduction( DirCount ) 
      END IF
      CALL Info(Caller,'Number of dofs set for '//TRIM(Name)//': '&
          //I2S(DirCount),Level=12)
    END IF
      
!------------------------------------------------------------------------------

  CONTAINS

    ! For boundary element in a DG mesh get the DG indexes from itself, or
    ! from the parent elements.
    !-------------------------------------------------------------------------
    SUBROUTINE PickDgIndexes(Element,DgIndexes)
      TYPE(Element_t), POINTER :: Element
      INTEGER :: DGIndexes(:)
      
      TYPE(Element_t), POINTER :: Parent
      INTEGER :: i,j,lr,n,m

      n = Element % Type % NumberOfNodes
       
      IF( ASSOCIATED( Element % DgIndexes ) ) THEN
        DGIndexes(1:n) = Element % DGIndexes(1:n)
        IF(ANY(DgIndexes(1:n) == 0) ) THEN
          ! If this fails, maybe used the 2nd strategy always...
          CALL Fatal('PickDgIndexes','Some of BC element DG indexes is zero!')
        END IF
      ELSE IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
        m = 0
        DO lr=1,2
          IF(lr==1) THEN
            Parent => Element % BoundaryInfo % Left
          ELSE
            Parent => Element % BoundaryInfo % Right
          END IF
          IF(.NOT. ASSOCIATED( Parent ) ) CYCLE
          IF(.NOT. ASSOCIATED( Parent % DGIndexes ) ) CYCLE
          IF(ANY(Perm(Parent % DGIndexes) == 0)) CYCLE
          
          DO i=1,Element % TYPE % NumberOfNodes
            DO j=1,Parent % TYPE % NumberOfNodes
              IF( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
                DGIndexes(i) = Parent % DGIndexes(j)
                m = m+1
                EXIT
              END IF
            END DO
          END DO
          EXIT
        END DO
        IF(m /= n) THEN
          CALL Fatal('PickDgIndexes','Could not find all DG Indexes for BC element')
        END IF
      END IF
      
    END SUBROUTINE PickDgIndexes

    
     ! Check n-t node setting element
     !-------------------------------
    SUBROUTINE CheckNTElement(n,elno)
      INTEGER :: n,elno
      INTEGER :: i,j,k,l,m,dim,kmax
      LOGICAL :: found
      REAL(KIND=dp) :: Condition(n), RotVec(3)
      
      dim = CoordinateSystemDimension()

      IF ( DOF <= 0 ) RETURN
      
      IF ( NT % NormalTangentialNOFNodes == 0 ) RETURN
      IF ( ALL(NT % BoundaryReorder(Indexes(1:n))<1) ) RETURN
      IF ( .NOT. ListCheckPresent(ValueList, Name) ) RETURN
      IF ( ListGetLogical(ValueList,NT % NormalTangentialName,Found) ) RETURN

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
          m = NT % BoundaryReorder(Indexes(j))
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
    SUBROUTINE SetElementValues(n)
!------------------------------------------------------------------------------
      INTEGER :: n
      INTEGER :: i,j,k,l,m,dim,kmax,lmax
      LOGICAL :: CheckNT,found
      REAL(KIND=dp) :: Condition(n), Work(n), RotVec(3)
      
      dim = CoordinateSystemDimension()

      IF ( DOF > 0 ) THEN
        IF (Model % Solver % DG) THEN
          Work(1:n)  = ListGetReal( ValueList, Name, n, Element % NodeIndexes, gotIt )
        ELSE
          Work(1:n)  = ListGetReal( ValueList, Name, n, Indexes, gotIt )
        END IF
        IF ( .NOT. GotIt ) THEN
          Work(1:n)  = ListGetReal( ValueList, Name(1:nlen) // ' DOFs', n, Indexes, gotIt )
        END IF
      ELSE
        CALL ListGetRealArray( ValueList, Name, WorkA, n, Indexes, gotIt )
      END IF
      
      IF ( gotIt ) THEN
        IF ( Conditional ) THEN
          IF (Model % Solver % DG) THEN
            Condition(1:n) = ListGetReal( ValueList, CondName, n, Element % NodeIndexes, gotIt )
          ELSE
            Condition(1:n) = ListGetReal( ValueList, CondName, n, Indexes, gotIt )
          END IF
          Conditional = Conditional .AND. GotIt
        END IF

        !
        ! Check for nodes belonging to n-t boundary getting set by other bcs.
        ! Here we don't need to track p-bubbles as they are not shared by
        ! many BCs. 
        ! -------------------------------------------------------------------
        CheckNT = .FALSE.
        IF ( NT % NormalTangentialNOFNodes>0 .AND. DOF>0 ) THEN
          CheckNT = .TRUE.
          IF ( ALL(NT % BoundaryReorder(Indexes(1:n))<1) ) CheckNT = .FALSE.
          IF ( ListGetLogical(ValueList,NT % NormalTangentialName,Found)) CheckNT=.FALSE.
        END IF
        
        DO j=1,n
          IF ( Conditional ) THEN
            IF( Condition(j) < 0.0_dp ) CYCLE
          END IF

          k = Perm(Indexes(j))
          IF ( k > 0 ) THEN
            
            IF ( DOF>0 ) THEN
              m = 0
              IF ( NT % NormalTangentialNOFNodes>0 ) m = NT % BoundaryReorder(Indexes(j))
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
                  b(lmax) = 0._dp

                  IF( .NOT. OffDiagonal ) THEN
                    b(lmax) = b(lmax) + Work(j) 
                  END IF

                  ! Consider all components of the cartesian vector mapped to the 
                  ! N-T coordinate system. Should this perhaps have scaling included?
                  DirCount = DirCount + 1
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
!------------------------------------------------------------------------------
    SUBROUTINE SetLumpedRows(ind0,n)
!------------------------------------------------------------------------------
      INTEGER :: ind0,n
      INTEGER :: ind,i,j,k,k0,l
      REAL(KIND=dp) :: Coeff
      ! -------------------------------------------------------------------        
      
      DO j=1,n
        ind = Indexes(j)

        IF( LumpedNodeSet(ind) ) CYCLE
        LumpedNodeSet(ind) = .TRUE.

        IF ( DOF > 0 ) THEN
          k0 = Offset + NDOFs * (Perm(ind0)-1) + DOF
          k = OffSet + NDOFs * (Perm(ind)-1) + DOF

          IF( GotMult ) THEN
            Coeff = Mult(j) / MaxMult
          ELSE
            Coeff = 1.0_dp
          END IF
          
          CALL MoveRow( A, k, k0, MoveCoeff )
          b(k0) = b(k0) + MoveCoeff * b(k)

          CALL AddToMatrixElement( A, k, k, 1.0_dp )
          CALL AddToMatrixElement( A, k, k0, -Coeff )
          b(k) = 0.0_dp
        ELSE
          DO l = 1, NDOFs
            k0 = Offset + NDOFs + (Perm(ind0)-1) * DOF
            k = OffSet + NDOFs * (Perm(ind)-1) + l

            IF( GotMult ) THEN
              Coeff = Mult(j) / MaxMult
            ELSE
              Coeff = 1.0_dp
            END IF
              
            CALL MoveRow( A, k, k0, 1.0_dp )
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
!> Set values related to a rigid plane boundary such that any node on the boundary
!> is expressed as linear combinination of the selected three (or two if on line)
!> nodes.
!------------------------------------------------------------------------------
    SUBROUTINE SetRigidRows(inds0,bcind,n)
!------------------------------------------------------------------------------
      INTEGER :: inds0(:)
      INTEGER :: bcind
      INTEGER :: n

      INTEGER :: bcind0 =  0
      INTEGER :: ind,i,j,k,k0
      REAL(KIND=dp) :: Coeff, Weights(3)
      REAL(KIND=dp) :: BaseCoord(3,3),r1(3),r2(3),Coord(3),dCoord(3),Amat(2,2),A0mat(2,2),bvec(2)
      
      SAVE bcind0, BaseCoord, A0mat, r1, r2
!-------------------------------------------------------------------        

      IF(bcind /= bcind0 ) THEN
        BaseCoord = 0.0_dp
        DO i=1,dim
          j = inds0(i)
          BaseCoord(i,1) = Mesh % Nodes % x(j)
          BaseCoord(i,2) = Mesh % Nodes % y(j)
          BaseCoord(i,3) = Mesh % Nodes % z(j)
        END DO
        bcind0 = bcind

        r1 = BaseCoord(2,:) - BaseCoord(1,:)
        Amat(1,1) = SUM(r1*r1)
        
        IF( dim == 3 ) THEN
          r2 = BaseCoord(3,:) - BaseCoord(1,:)
          Amat(1,2) = SUM(r1*r2)
          Amat(2,1) = Amat(1,2)
          Amat(2,2) = SUM(r2*r2)
        END IF

        A0mat = Amat
        bcind0 = bcind
      END IF
                   
      DO j=1,n
        ind = Indexes(j)

        IF( LumpedNodeSet(ind) ) CYCLE
        LumpedNodeSet(ind) = .TRUE.
        
        Coord(1) = Mesh % Nodes % x(ind)
        Coord(2) = Mesh % Nodes % y(ind)
        Coord(3) = Mesh % Nodes % z(ind)

        dCoord = Coord - BaseCoord(1,:)
        
        bvec(1) = SUM( dCoord * r1 )
        IF( dim == 3 ) THEN
          bvec(2) = SUM( dCoord * r2 )
        END IF

        IF( dim == 2 ) THEN
          bvec(1) = bvec(1) / A0mat(1,1)
          Weights(2) = bvec(1)
          Weights(1) = 1.0_dp - Weights(2)
        ELSE
          Amat = A0mat          
          CALL LUSolve(2,Amat,bvec)          
          Weights(2:3) = bvec(1:2)
          Weights(1) = 1.0_dp - SUM(bvec(1:2))
        END IF

        DO l = 1, dim
          k = OffSet + NDOFs * (Perm(ind)-1) + l    

          ! Distribute row in accordance with the weights
          DO m = 1, dim
            k0 = Offset + NDOFs * (Perm(inds0(m))-1) + l          
            Coeff = Weights(m)

            b(k0) = b(k0) + Coeff * b(k)
            IF( m < dim ) THEN
              ! This does not nullify the row
              CALL MoveRow( A, k, k0, Coeff, 1.0_dp )
            ELSE
              ! Now also nullify the row
              CALL MoveRow( A, k, k0, Coeff )
              b(k) = 0.0_dp
            END IF
          END DO

          ! Express the node as linear combination of the base nodes
          DO m = 1,dim
            k0 = Offset + NDOFs * (Perm(inds0(m))-1) + l          
            Coeff = Weights(m)            
            CALL AddToMatrixElement( A, k, k0, -Coeff )
          END DO
          CALL AddToMatrixElement( A, k, k, 1.0_dp )
        END DO

      END DO

!------------------------------------------------------------------------------
    END SUBROUTINE SetRigidRows
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
          ! This can (for example) happen with "Target Node" - keyword used in parallel...
          IF(NodeIndexes(j)>SIZE(Perm) .OR. NodeIndexes(j)<=0) CYCLE

          k = Perm(NodeIndexes(j))
          IF( k == 0 ) CYCLE

          IF ( Conditional .AND. Condition(j) < 0.0d0 ) CYCLE

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


      IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
        ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
        A % ConstrainedDOF = .FALSE.
      END IF
      
      IF(.NOT. ALLOCATED(A % Dvalues)) THEN
        ALLOCATE(A % Dvalues(A % NumberOfRows))
        A % Dvalues = 0._dp
      END IF
      
      k = ind
      IF (ApplyPerm) k = Perm(ind)
      IF( k == 0 ) RETURN
      
      k = OffSet + NDOFs * (k-1) + DOF

      ! Do not add non-zero entries to pure halo nodes which are not associated with the partition.
      ! These are nodes could be created by the -halobc flag in ElmerGrid.
      IF( Parallel ) THEN
        IF( ASSOCIATED( A % ParallelInfo ) ) THEN
          IF( .NOT. ANY( A % ParallelInfo % NeighbourList(k) % Neighbours == ParEnv % MyPe ) ) THEN
            RETURN
          END IF
        END IF
      END IF

      DirCount = DirCount + 1
      
      IF( .NOT. OffDiagonal ) THEN
        A % Dvalues(k) =  val
      END IF
      A % ConstrainedDOF(k) = .TRUE.

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
    LOGICAL :: Parallel
!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)
    BC => Model % BCs(This) % Values

    IF ( ListGetLogical( BC,& 
        'Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      IF( ListGetLogical( BC,'Antisymmetric BC',GotIt ) ) THEN
        Scale = 1.0_dp
      ELSE
        Scale = -1.0_dp
      END IF
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
          A % Dvalues(k) = 0._dp
          A % ConstrainedDOF(k) = .TRUE.
          
          DO l = Projector % Rows(i), Projector % Rows(i+1)-1
            IF ( Projector % Cols(l) <= 0 ) CYCLE
            m = Perm( Projector % Cols(l) )
            IF ( m > 0 ) THEN
              m = NDOFs * (m-1) + DOF
              A % Dvalues(k) = A % Dvalues(k) - Scale * Projector % Values(l) * &
                  Var % Values(m)
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
            'Allocating mortar BCs for solver',Level=8)
        ALLOCATE( Model % Solver % MortarBCs( Model % NumberOfBCs + Model % NumberOfBodyForces ) )
        DO i=1, Model % NumberOfBCs + Model % NumberOfBodyForces
          Model % Solver % MortarBCs(i) % Projector => NULL()
        END DO
      END IF
      
      IF( ASSOCIATED( Projector, &
          Model % Solver % MortarBCs(This) % Projector) ) THEN
        CALL Info('SetPeriodicBoundariesPass1','Using existing projector: '&
            //I2S(This),Level=8)
        RETURN
      END IF
      
      Model % Solver % MortarBCs(This) % Projector => Projector
      CALL Info('SetPeridociBoundariesPass1','Using projector as mortar constraint: '&
          //I2S(This),Level=8)

      MortarBC => Model % Solver % MortarBCs(This)      
      IF( Jump ) THEN
        IF( ASSOCIATED( MortarBC % Diag ) ) THEN
          IF( SIZE( MortarBC % Diag ) < NDofs * Projector % NumberOfRows ) THEN
            DEALLOCATE( MortarBC % Diag ) 
          END IF
        END IF
        IF( .NOT. ASSOCIATED( MortarBC % Diag ) ) THEN
          CALL Info('SetPeriodicBoundariesPass1','Allocating projector mortar diag',Level=10)
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
          CALL Info('SetPeriodicBoundariesPass1','Allocating projector mortar rhs',Level=10)
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
        CALL Info('SetPeriodicBoundariesPass1','Allocating projector mortar perm',Level=10)
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
                          -scale*Projector % Values(l) * A % Values(nn) ) 
                   IF (ASSOCIATED(F % Values)) THEN
                     CALL AddToMatrixElement( F, m, F % Cols(nn), &
                          -scale*Projector % Values(l) * F % Values(nn) )
                   END IF
                END DO
                b(m)=b(m) - scale*Projector % Values(l)*b(k) 
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
      IF( ListGetLogical( BC,'Antisymmetric BC',GotIt ) ) THEN
        Scale = 1.0_dp
      ELSE
        Scale = -1.0_dp
      END IF
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
           val = -weight * coeff 
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
               CALL AddToMatrixElement( A,k,m,val * Projector % Values(l) )
             END IF
          END DO

          b(k) = b(k) - ValueOffset 
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
    CALL Info('SetFrictionBoundaries','Number of friction nodes: '//I2S(n),Level=10)

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
           //I2S(bc_ind),Level=7)
    
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



  
!------------------------------------------------------------------------------
!> Prepare to set Dirichlet conditions for attachment DOFs in the case of
!> component mode synthesis
!------------------------------------------------------------------------------
  SUBROUTINE SetConstraintModesBoundaries( Model, Solver, A, b, &
      Name, NDOFs, Perm )
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model              !< current model structure
    TYPE(Solver_t), TARGET :: Solver    !< current solver structure 
    TYPE(Matrix_t), POINTER :: A        !< global matrix
    REAL(KIND=dp) :: b(:)               !< global RHS vector
    CHARACTER(LEN=*) :: Name            !< name of the dof to be set
    INTEGER :: NDOFs                    !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)          !< The node reordering info, this has been generated at the
                                        !< beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    INTEGER :: i,t,t1,t2,u,j,k,k2,l,l2,n,ent_id,nlen,NormalInd,Ncomplex
    LOGICAL :: Found, IgnoreP, HaveP, Passive, Parallel, DoIt
    TYPE(Solver_t), POINTER :: pSolver
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: Indexes(50),nb,poffset
    INTEGER, POINTER :: NodeIndexes(:), MasterBodies(:)
    INTEGER, ALLOCATABLE :: BCPerm(:)
    INTEGER :: NoModes, NoEntities
    LOGICAL :: ExternalLoop, BFMode, BCMode, LumpedMode, CompMode, &
        EmWaveMode, RhsMode, CoilMode, ComplexMode, SingleMode, Stat, GotDir
    REAL(KIND=dp) :: ModeDir(NDofs)
    REAL(KIND=dp), POINTER :: HelperArray(:,:)
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    CHARACTER(*), PARAMETER :: Caller = 'SetConstraintModesBoundaries'

!------------------------------------------------------------------------------

    CALL Info(Caller,'Checking constraint modes boundaries!',Level=20)
    
    nlen = LEN_TRIM(Name)
    Mesh => Solver % Mesh
    Var => Solver % Variable     

    ! Initially this is -1 and and hence the 2nd call is fast if no modes are present
    IF( Solver % NumberOfConstraintModes == 0 ) RETURN

    ! Check in which list we have the constraint modes set
    BCMode = ListCheckPrefixAnyBC( Model,'Constraint Mode')
    BFMode = ListCheckPrefixAnyBodyForce( Model,'Constraint Mode')
    CompMode = ListCheckPrefixAnyComponent( Model,'Constraint Mode')

    IF( .NOT. (BCMode .OR. BFMode .OR. CompMode) ) THEN
      Solver % NumberOfConstraintModes = 0
      RETURN
    END IF
    
    ExternalLoop = ListGetLogical( Solver % Values,'Nonlinear System Constraint Modes', Found ) .OR. &
        ListGetLogical( Solver % Values,'Steady State Constraint Modes', Found ) .OR. &
        ListGetLogical( Solver % Values,'Run Control Constraint Modes', Found ) .OR. &
        ListGetLogical( Model % Control,'Constraint Modes Analysis', Found )
       
    EmWaveMode = ListGetLogical( Solver % Values,'Constraint Modes EM Wave',Found )
    CoilMode = ListGetLogical( Solver % Values,'Constraint Modes Coils',Found ) 
    RhsMode = ListGetLogical(Solver % Values,'Constraint Modes Rhs',Found ) 
    LumpedMode = ListGetLogical(Solver % Values,'Constraint Modes Lumped',Found )
    
    ! These work on the rhs vector, not Dirichlet values.
    RhsMode = RhsMode .OR. EmWaveMode .OR. CoilMode
    LumpedMode = LumpedMode .OR. EmWaveMode .OR. CoilMode
    
    Element => Mesh % Elements(1)
    pSolver => Solver
    HaveP = isActivePElement(Element,pSolver)
    
    IgnoreP = .FALSE.
    IF( HaveP ) THEN
      IgnoreP = ListGetLogical( Solver % Values,'Ignore Constraint Modes p',Found )  
    END IF      

    CALL Info(Caller,'Setting constraint modes boundaries for variable: '&
        //TRIM(Name),Level=7)

    Parallel = ( ParEnv % PEs > 1 ) .AND. ( .NOT. Mesh % SingleMesh ) 
    
    ! Allocate the indeces for the constraint modes.
    ! We may be revisiting the routine, and the mesh may have changed...
    IF( ASSOCIATED( Var % ConstraintModesIndeces ) ) THEN
      IF( SIZE( Var % ConstraintModesIndeces ) /= A % NumberOfRows ) THEN
        DEALLOCATE( Var % ConstraintModesIndeces )
        IF( RhsMode ) THEN
          IF(ASSOCIATED( Var % ConstraintModesWeights ) ) THEN
            DEALLOCATE( Var % ConstraintModesWeights )
          END IF
        END IF
      END IF
      DoIt = ListGetLogical( Solver % Values,'Constraint Modes Fluxes Norm', Found )
      IF(.NOT. Found) DoIt = ListGetLogical( Solver % Values,'Constraint Modes Matrix Norm', Found )
      IF(DoIt) CALL ListAddConstReal( CurrentModel % Simulation,&
          'res: '//TRIM(Solver % Variable % Name)//' lumped matrix norm',0.0_dp)
    END IF
    IF(.NOT. ASSOCIATED( Var % ConstraintModesIndeces ) ) THEN
      ALLOCATE( Var % ConstraintModesIndeces( A % NumberOfRows ) )
      IF( RhsMode ) THEN
        ALLOCATE( Var % ConstraintModesWeights( A % NumberOfRows ) )
      END IF
    END IF
    
    Var % ConstraintModesIndeces = 0
    IF( RhsMode ) Var % ConstraintModesWeights = 0.0_dp
    
    IF( BCMode ) THEN
      NoEntities = Model % NumberOfBCs
    ELSE IF( BFMode ) THEN
      NoEntities = Model % NumberOfBodyForces
    ELSE IF(CompMode) THEN
      NoEntities = Model % NumberOfComponents
    ELSE
      CALL Fatal(Caller,'Uknown list for constraint modes!')
    END IF
    CALL Info(Caller,'Number of list entities to check: '//I2S(NoEntities),Level=20)
    
    ALLOCATE( BCPerm( NoEntities ) )
    BCPerm = 0    
    j = 0

    DO ent_id = 1,NoEntities
      IF( BCMode ) THEN
        BC => Model % BCs(ent_id) % Values
      ELSE IF(BFMode ) THEN
        BC => Model % BodyForces(ent_id) % Values
      ELSE IF(CompMode) THEN
        BC => Model % Components(ent_id) % Values
      END IF
        
      k = ListGetInteger( BC,'Constraint Mode', Found )
      IF(.NOT. Found ) THEN
        k = ListGetInteger( BC,&
            'Constraint Mode '// Name(1:nlen), Found )
      END IF
      IF( Found ) THEN
        IF( k == 0 ) k = -1  ! Ground gets negative value
        BCPerm(ent_id) = k        
      ELSE
        DoIt = ListGetLogical( BC,'Constraint Modes', Found )
        IF(.NOT. Found ) THEN
          DoIt = ListGetLogical( BC,'Constraint Modes ' // Name(1:nlen), Found ) 
        END IF
        IF(DoIt) THEN
          j = j + 1
          BCPerm(ent_id) = j
        END IF
      END IF
    END DO
    
    j = MAXVAL( BCPerm )
    CALL Info(Caller,'Number of active constraint modes boundaries: '&
        //I2S(j),Level=7)
    IF( j == 0 ) THEN
      CALL Fatal(Caller,&
          'Constraint Modes Analysis requested but no constrained BCs given!')
    END IF

    
    ComplexMode = LIstGetLogical( Solver % Values,'Linear System Complex',Found)
    IF( ComplexMode ) THEN
      Ncomplex = 2
      CALL Info(Caller,'Assuming complex valued system for constraint modes',Level=12)
    ELSE
      Ncomplex = 1
      CALL Info(Caller,'Assuming real valued system for constraint modes',Level=12)
    END IF

    SingleMode = ListGetLogical( Solver % Values,'Constraint Modes Single',Found ) 
    IF(.NOT. Found ) THEN
      SingleMode = ListCheckPresentAnyBC( Model,'Constraint Mode Direction' )            
    END IF
    IF(SingleMode) THEN
      CALL Info(Caller,'Setting constraint modes for all components at once!',Level=12)
    END IF
    
    IF( SingleMode ) THEN
      NoModes = j
    ELSE      
      NoModes = NDOFS * j  / Ncomplex
    END IF
      
    IF( BcMode ) THEN
      t1 = Mesh % NumberOfBulkElements+1
      t2 = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    ELSE
      t1 = 1
      t2 = Mesh % NumberOfBulkElements
    END IF
    
    DO t = t1, t2
      Element => Mesh % Elements(t)      
      
      ModeDir = 1.0_dp
      IF( BCMode ) THEN
        DO ent_id = 1,Model % NumberOfBCs
          IF ( Element % BoundaryInfo % Constraint == Model % BCs(ent_id) % Tag ) EXIT
        END DO
        IF( ent_id > Model % NumberOfBCs ) CYCLE        
        HelperArray => ListGetConstRealArray( Model % BCs(ent_id) % Values,'Constraint Mode Direction',Found)
        IF(Found) ModeDir(1:ndofs) = HelperArray(1:ndofs,1)
      ELSE IF( BFMode ) THEN
        ent_id = ListGetInteger( Model % Bodies( Element % BodyId ) % Values,'Body Force',Found)
        IF( ent_id == 0) CYCLE
        HelperArray => ListGetConstRealArray( Model % BodyForces(ent_id) % Values,'Constraint Mode Direction',Found)
        IF(Found) ModeDir(1:ndofs) = HelperArray(1:ndofs,1)
      ELSE        
        DO ent_id=1,Model % NumberOfComponents
          IF(BCPerm(ent_id) == 0) CYCLE
          MasterBodies => ListGetIntegerArray( Model % Components(ent_id) % Values,'Master Bodies',Found)
          IF(.NOT. Found) CYCLE
          IF( ANY( MasterBodies == Element % BodyId ) ) EXIT          
        END DO
        IF(ent_id > Model % NumberOfComponents ) CYCLE
      END IF

      ! Look-up table for quick determiation
      IF( BCPerm(ent_id) == 0 ) CYCLE

      n = Element % Type % NumberOfNodes

      ! This is used in standard setting of Dirichlet BCs
      nb = mGetElementDOFs( Indexes, Element, Solver )
            
      IF( RhsMode ) THEN
        BLOCK 
          REAL(KIND=dp) :: Basis(n), Weight(n), detJ, s          
          IP = GaussPoints( Element )
          CALL CopyElementNodesFromMesh( Nodes, Solver % Mesh, n, Element % NodeIndexes)
                    
          DO i=1,IP % n 
            stat = ElementInfo( Element, Nodes, IP % U(i), IP % V(i), &
              IP % W(i), detJ, Basis )          
            s = IP % s(i) * DetJ
            DO k=1,NDOFS                            
              Var % ConstraintModesWeights(NDOfs*(Perm(Indexes(1:n))-1)+k) = &
                  Var % ConstraintModesWeights(NDOfs*(Perm(Indexes(1:n))-1)+k) + ModeDir(k) * s * Basis(1:n)
            END DO
          END DO
        END BLOCK
      END IF
            
      ! If for some reason we do not want to set the P dofs to zero
      IF(IgnoreP) nb = n

      !PRINT *,'BCPerm:',MINVAL(BCPerm), MAXVAL(BCPerm), COUNT(BCPerm>0)
      
      ! For vector valued problems treat each component as separate dof
      DO k=1,NDOFs       
        IF( SingleMode ) THEN
          j = BCPerm(ent_id)
        ELSE
          j = NDOFS*(BCPerm(ent_id)-1)+k
        END IF
        DO l=1,nb
          ! The index to constrain
          IF( Perm(Indexes(l)) == 0 ) CYCLE
          l2 = NDOFS*(Perm(Indexes(l))-1)+k         
          Var % ConstraintModesIndeces(l2) = j
        END DO
      END DO
      !PRINT *,'Single range:',SingleMode,MAXVAL(Var % ConstraintModesIndeces), COUNT(Var % ConstraintModesIndeces > 0)
    END DO
    DEALLOCATE(BCPerm)

    ! Some single node or edge could stretch to the surface even though it is not
    ! part of any boundary element in parallel. Hence we need to communicate the tag. 
    IF( Parallel ) THEN
      CALL Info(Caller,'Communicating tags for constraint modes',Level=20)
      CALL CommunicateParallelSystemTag(A % ParallelInfo,Itag = Var % ConstraintModesIndeces,ParOper=2)
    END IF
      
    ! Set the p dofs to negative since we don't ever want to set them to one!
    ! Note that there are some dofs related to ground that are already negative.
    ! Hence ground and p-pubbles are treated alike. 
    IF(HaveP .AND. .NOT. IgnoreP) THEN
      poffset = 2*(NoModes + 1)
      DO l=Mesh % NumberOfNodes+1, SIZE(Perm)
        j = Perm(l)
        IF(j==0) CYCLE
        DO k=1,NDOFs       
          l2 = NDOFS*(j-1)+k                 
          
          ! Subtract a big enough number of the constraint modes so that they are always negative. 
          IF( Var % ConstraintModesIndeces(l2) /= 0 ) THEN
            Var % ConstraintModesIndeces(l2) = Var % ConstraintModesIndeces(l2) - poffset 
          END IF
        END DO
      END DO
    END IF
    
    ! This is just for information.
    ! We show how Dirichlet conditions are split among nodal and p dofs, and ground. 
    IF( InfoActive(12) ) THEN
      DO i=0,NoModes
        j = i
        IF(i==0) j=-1

        k = COUNT( Var % ConstraintModesIndeces == j )
        IF(k>0) THEN
          CALL Info('SetConstaintModesBoundaries',&
              'Mode '//I2S(i)//' has '//I2S(k)//' dofs')
          IF(.NOT. IgnoreP .AND. HaveP ) THEN
            j = j - poffset
            k = COUNT( Var % ConstraintModesIndeces == j )
            IF(k>0) THEN
              CALL Info('SetConstaintModesBoundaries',&
                  'Mode '//I2S(i)//' has additional '//I2S(k)//' p-dofs')
            END IF
          END IF
        END IF
      END DO
    END IF

    
    ! The constraint modes can be either lumped or not.
    ! If they are not lumped then mark each individually
    IF( .NOT. LumpedMode ) THEN
      j = 0
      DO i=1,A % NumberOfRows
        IF( Var % ConstraintModesIndeces(i) > 0 ) THEN
          j = j + 1
          Var % ConstraintModesIndeces(i) = j
        END IF
      END DO
      CALL Info(Caller,'Number of active constraint modes: '&
          //I2S(j),Level=7)
      NoModes = j 
    END IF
        
    ! Manipulate the boundaries such that we need to modify only the r.h.s. in the actual linear solver
    ! Do not manipulate if we are setting fluxes!
    IF( .NOT. (ComplexMode .OR. RhsMode ) ) THEN
      WHERE( Var % ConstraintModesIndeces /= 0 ) 
        A % ConstrainedDOF = .TRUE.
        A % DValues = 0.0_dp
      END WHERE
    END IF
    
    Solver % NumberOfConstraintModes = NoModes

    ! We may want to save the results for postprocessing even when we do this in
    ! one sweep. If we iterate over nonlinear, steady state of run control then automatically
    ! we can have access to all components in saving. 
    IF(.NOT. ExternalLoop ) THEN
      Var % NumberOfConstraintModes = NoModes
      ALLOCATE( Var % ConstraintModes( Var % NumberOfConstraintModes, A % NumberOfRows ) )
      Var % ConstraintModes = 0.0_dp
    END IF
          
    CALL Info(Caller,'All done',Level=10)

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
      A % ConstrainedDOF(PermIndex) = .TRUE.
      A % DValues(PermIndex) = NodeValue      
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichletPoint
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
!> Set distributed loads to vector b.
!------------------------------------------------------------------------------
  SUBROUTINE SetNodalSources( Model, Mesh, SourceName, dofs, Perm, GotSrc, SrcVec )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model  !< The current model structure
    TYPE(Mesh_t)  :: Mesh    !< The current mesh structure
    CHARACTER(LEN=*) :: SourceName   !< Name of the keyword setting the source term
    INTEGER :: DOFs                  !< The total number of DOFs for this equation
    INTEGER :: Perm(:)               !< The node reordering info
    LOGICAL :: GotSrc                !< Did we get something?
    REAL(KIND=dp) :: SrcVec(:)       !< The assemblied source vector
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: i,j,k,t,n,bc,bf,FirstElem,LastElem,nlen
    LOGICAL :: Found,AnyBC,AnyBF,Axisymmetric
    REAL(KIND=dp) :: Coeff
    REAL(KIND=dp), ALLOCATABLE :: FORCE(:,:)
    LOGICAL, ALLOCATABLE :: ActiveBC(:), ActiveBF(:)
    TYPE(ValueList_t), POINTER :: ValueList
    INTEGER, POINTER :: Indexes(:)
    CHARACTER(*), PARAMETER :: Caller = 'SetNodalSources'    
    LOGICAL :: Parallel
    
    nlen = LEN_TRIM(SourceName)
       
    CALL Info(Caller,'Checking for generalized source terms: '&
        //SourceName(1:nlen),Level=15)

    ALLOCATE( ActiveBC(Model % NumberOfBCs), ActiveBF(Model % NumberOfBodyForces) )
    
    ! First make a quick test going through the short boundary condition and
    ! body force lists.
    ActiveBC = .FALSE.
    DO BC=1,Model % NumberOfBCs
      IF(.NOT. ListCheckPresent( Model % BCs(BC) % Values,'Target Boundaries')) CYCLE
      ActiveBC(BC) = ListCheckPrefix( Model % BCs(BC) % Values, SourceName(1:nlen) )
    END DO

    ActiveBF = .FALSE.
    DO bf=1,Model % NumberOFBodyForces
      ActiveBF(bf) = ListCheckPrefix( Model % BodyForces(bf) % Values, SourceName(1:nlen) ) 
    END DO

    AnyBC = ANY(ActiveBC)
    AnyBF = ANY(ActiveBF)    
    
    GotSrc = (AnyBC .OR. AnyBF)
    IF(.NOT. GotSrc ) RETURN

    CALL Info(Caller,'Assembling generalized source terms: '&
        //SourceName(1:nlen),Level=10)

    
    AxiSymmetric = ( CurrentCoordinateSystem() /= Cartesian )

    Parallel = ( ParEnv % PEs > 1 ) .AND. ( .NOT. Mesh % SingleMesh ) 
    
    ! Only loop over BCs and BFs if needed. Here determine the loop.
    FirstElem = HUGE( FirstElem )
    LastElem = 0
    IF(AnyBC) THEN
      FirstElem = Mesh % NumberOfBulkElements + 1
      LastElem = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    END IF
    IF(AnyBF) THEN
      FirstElem = 1
      LastElem = MAX( LastElem, Mesh % NumberOfBulkElements )
    END IF

    n = Mesh % MaxElementNodes
    ALLOCATE( FORCE(dofs,n) )
    FORCE = 0.0_dp
    
    ! Here do the actual assembly loop. 
    DO t=FirstElem, LastElem 
      Element => Mesh % Elements(t)
      Indexes => Element % NodeIndexes
      n = Element % Type % NumberOfNodes

      IF( t > Mesh % NumberOfBulkElements ) THEN
        Found = .FALSE.
        DO BC=1,Model % NumberOfBCs
          IF( .NOT. ActiveBC(BC) ) CYCLE
          IF ( Element % BoundaryInfo % Constraint == Model % BCs(BC) % Tag ) THEN
            Found = .TRUE.
            EXIT
          END IF
        END DO
        IF(.NOT. Found) CYCLE
        ValueList => Model % BCs(BC) % Values
      ELSE                
        bf = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force',Found)        
        IF(.NOT. Found) CYCLE
        IF(.NOT. ActiveBF(bf) ) CYCLE
        ValueList => Model % BodyForces(bf) % Values
      END IF

      ! In parallel we may have halos etc. By default scaling is one.       
      Coeff = ParallelScalingFactor()
      IF(ABS(Coeff) < TINY(Coeff)) CYCLE

      CALL LocalSourceAssembly(Element, dofs, FORCE )

      DO i=1,dofs
        SrcVec(dofs*(Perm(Indexes)-1)+i) = SrcVec(dofs*(Perm(Indexes)-1)+i) + FORCE(i,1:n)
      END DO
    END DO
      
  
  CONTAINS

!------------------------------------------------------------------------------
    FUNCTION ParallelScalingFactor() RESULT ( Coeff ) 
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Coeff
      TYPE(Element_t), POINTER :: P1, P2
      
      ! Default weight
      Coeff = 1.0_dp
      
      IF ( Parallel ) THEN
        IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
          P1 => Element % BoundaryInfo % Left
          P2 => Element % BoundaryInfo % Right
          IF ( ASSOCIATED(P1) .AND. ASSOCIATED(P2) ) THEN
            IF ( P1 % PartIndex /= ParEnv % myPE .AND. &
                P2 % PartIndex /= ParEnv % myPE ) THEN
              Coeff = 0.0_dp            
            ELSE IF ( P1 % PartIndex /= ParEnv % myPE .OR. &
                P2 % PartIndex /= ParEnv % myPE ) THEN
              Coeff = 0.5_dp
            END IF
          ELSE IF ( ASSOCIATED(P1) ) THEN
            IF ( P1 % PartIndex /= ParEnv % myPE ) Coeff = 0.0_dp
          ELSE IF ( ASSOCIATED(P2) ) THEN
            IF ( P2 % PartIndex /= ParEnv % myPE ) Coeff = 0.0_dp
          END IF
        ELSE IF ( Element % PartIndex/=ParEnv % myPE ) THEN
          Coeff = 0.0_dp
        END IF
      END IF

    END FUNCTION ParallelScalingFactor
!------------------------------------------------------------------------------

    
!------------------------------------------------------------------------------
    SUBROUTINE LocalSourceAssembly(Element, dofs, FORCE)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dofs
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: FORCE(:,:)
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),ElemSource(:,:)
    REAL(KIND=dp) :: weight, SourceAtIp, DetJ
    INTEGER, POINTER :: Indexes(:)
    LOGICAL :: Stat,Found
    INTEGER :: i,j,t,m,n,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes

    SAVE Nodes,Basis,ElemSource
!------------------------------------------------------------------------------

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      m = Mesh % MaxElementNodes
      ALLOCATE(ElemSource(dofs,m), Basis(m), Nodes % x(m), &
          Nodes % y(m), Nodes % z(m), STAT=allocstat)      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed in LocalMatrix')
      END IF
    END IF

    IP = GaussPoints( Element, PReferenceElement = .FALSE.)
    Indexes => Element % NodeIndexes
    n = Element % Type % NumberOfNodes

    Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
    Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
    Nodes % z(1:n) = Mesh % Nodes % z(Indexes)
      
    FORCE = 0._dp

    IF( dofs == 1 ) THEN
      ElemSource(1,1:n) = ListGetReal( ValueList,SourceName(1:nlen), n, Indexes )
    ELSE
      j = 0
      DO i=1,dofs
        ElemSource(i,1:n) = ListGetReal( ValueList,&
            SourceName(1:nlen)//' '//I2S(i), n, Indexes, Found )
        IF( Found ) j = j + 1
      END DO
      IF( j == 0 ) CALL Fatal(Caller,'Could not find for any component: '//SourceName(1:nlen) )
    END IF
    
    DO t=1,IP % n
      ! Basis function at the integration point:
      !-----------------------------------------
      stat = ElementInfo( Element, Nodes, &
          IP % U(t), IP % V(t), IP % W(t), detJ, Basis )
      Weight = IP % s(t) * DetJ

      IF ( AxiSymmetric ) THEN
        Weight = Weight * SUM( Nodes % x(1:n) * Basis(1:n) )
      END IF

      DO i=1,dofs
        SourceAtIP = SUM( ElemSource(i,1:n) * Basis(1:n) )
        FORCE(i,1:n) = FORCE(i,1:n) + &
            Weight * Basis(1:n) * SourceAtIp
      END DO     
    END DO
    
  END SUBROUTINE LocalSourceAssembly
  
!------------------------------------------------------------------------------
END SUBROUTINE SetNodalSources
!------------------------------------------------------------------------------

  

!------------------------------------------------------------------------------
!> Sets nodal loads directly to the matrix structure. 
!> The intended use for this is, for example, in multiphysics coupling where
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
    CHARACTER(:), ALLOCATABLE :: LoadName

    INTEGER, POINTER :: IndNodes(:)
    INTEGER :: NoNodes,NoDims,bf_id,nlen,NOFNodesFound
    REAL(KIND=dp), POINTER :: CoordNodes(:,:), DiagScaling(:),MinDist(:)
    REAL(KIND=dp) :: GlobalMinDist,Dist,Eps
    LOGICAL, ALLOCATABLE :: ActivePart(:), ActivePartAll(:), DoneLoad(:)
    LOGICAL :: NodesFound
    TYPE(ValueList_t), POINTER :: ValueList
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(*), PARAMETER :: Caller = 'SetNodalLoads'

    
    LoadName = TRIM(Name) // ' Load'
    nlen = LEN_TRIM(LoadName)
    
    CALL Info(Caller,'Checking for nodal loads for variable: '//TRIM(Name),Level=20)

    n = MAX(Model % NumberOfBCs, Model % NumberOFBodyForces) 
    ALLOCATE( ActivePart(n), ActivePartAll(n) )

    ALLOCATE( Indexes(Model % Solver % Mesh % MaxElementDOFs) )
!------------------------------------------------------------------------------
! Go through the boundaries
!------------------------------------------------------------------------------

    Mesh => Model % Mesh
    
    ActivePart = .FALSE.
    ActivePartAll = .FALSE.
    DO BC=1,Model % NumberOfBCs
      ValueList => Model % BCs(BC) % Values    
      IF (ListCheckPresent(ValueList,'Target Nodes') .OR. &
          ListCheckPresent(ValueList,'Target Coordinates')) CYCLE

      IF(ListCheckPresent(ValueList,'Target Boundaries') .OR. &
          ListCheckPresent(ValueList,'Name') ) THEN
        ActivePart(BC) = ListCheckPresent( Model % BCs(BC) % Values, LoadName )
        ActivePartAll(BC) = ListCheckPresent( &
            Model % BCs(BC) % Values, LoadName(1:nlen) // ' DOFs' )
      END IF
    END DO

    IF ( ANY(ActivePart) .OR. ANY(ActivePartAll) ) THEN
      CALL Info(Caller,'Setting nodal loads on boundaries: '//TRIM(LoadName),Level=9)
      ALLOCATE(DoneLoad( SIZE(b)/NDOFs) )
      DoneLoad = .FALSE.

      DO BC=1,Model % NumberOfBCs
        IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE

        DO t = Mesh % NumberOfBulkElements + 1, &
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

          Element => Mesh % Elements(t)
          IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BC) % Tag ) CYCLE
          
          Model % CurrentElement => Element
          IF ( ActivePart(BC) ) THEN
            n = Element % TYPE % NumberOfNodes
            Indexes(1:n) = Element % NodeIndexes
          ELSE
            n = mGetElementDOFs( Indexes )
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
      CALL Info(Caller,'Setting nodal loads on body force: '//TRIM(LoadName),Level=9)
      IF(.NOT. ALLOCATED(DoneLoad)) ALLOCATE(DoneLoad( SIZE(b)/NDOFs) )      
      DoneLoad = .FALSE.

      DO t = 1, Mesh % NumberOfBulkElements 
        Element => Mesh % Elements(t)
        bf_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force', GotIt)
        
        IF(.NOT. GotIt) CYCLE
        IF(.NOT. ActivePart(bf_id) .AND. .NOT. ActivePartAll(bf_id) ) CYCLE

        Model % CurrentElement => Element
        IF ( ActivePart(bf_id) ) THEN
          n = Element % TYPE % NumberOfNodes
          Indexes(1:n) = Element % NodeIndexes
        ELSE
          n = mGetElementDOFs( Indexes )
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

      ! At the first calling the list of coordinates is transformed to list of nodes.
      IF(.NOT. NodesFound) THEN
        IF( ListCheckPresent( ValueList,'Target Coordinates' ) ) THEN
          CALL TargetCoordinatesToTargetNodes( Model % Mesh, ValueList, NodesFound )
        END IF
      END IF
      
      IF(NodesFound) THEN           
        CALL Info(Caller,'Setting nodal loads on target nodes: '//TRIM(Name),Level=9)
        NodeIndexes => ListGetIntegerArray( ValueList,'Target Nodes')
        n = SIZE(NodeIndexes)

        IF(ANY(NodeIndexes>0)) THEN
          CALL SetPointLoads(n)
        END IF
      END IF

    END DO

    DEALLOCATE( Indexes )

    CALL Info(Caller,'Finished checking for nodal loads',Level=12)


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
               b(k) = b(k) + Work(j) 
             ELSE
               DO l=1,MIN( NDOFs, SIZE(Worka,1) )
                 k1 = NDOFs * (k-1) + l
                 b(k1) = b(k1) + WorkA(l,1,j) 
               END DO
             END IF
           END IF
         END DO
       END IF
       
     END SUBROUTINE SetElementLoads
     
     
     SUBROUTINE SetPointLoads(n)
       INTEGER :: n
       REAL(KIND=dp) :: Work(n)
       LOGICAL :: ImaginaryLoads

       IF(n<=0) RETURN
       ImaginaryLoads = ASSOCIATED(A % RHS_im)

       IF ( DOF > 0 ) THEN
         Work(1:n) = ListGetReal( ValueList, LoadName, n, NodeIndexes, gotIt )
       ELSE
         CALL ListGetRealArray( ValueList, LoadName, WorkA, n, NodeIndexes, gotIt )
       END IF
       
       IF ( GotIt ) THEN
         DO j=1,n
           IF ( NodeIndexes(j) > SIZE(Perm) ) THEN
             CALL Warn('SetPointLoads','Node number too large!')
             CYCLE
           END IF
           IF( NodeIndexes(j) == 0 ) THEN             
             CALL Warn('SetPointLoads','Node number is zero')
             CYCLE
           END IF
         
           k = Perm(NodeIndexes(j))
           IF ( k > 0 ) THEN
             IF ( DOF>0 ) THEN
               k = NDOFs * (k-1) + DOF
               b(k) = b(k) + Work(j) 
             ELSE
               DO l=1,MIN( NDOFs, SIZE(WorkA,1) )
                 k1 = NDOFs * (k-1) + l
                 b(k1) = b(k1) + WorkA(l,1,j) 
               END DO
             END IF
           END IF
         END DO
       END IF

       IF (ImaginaryLoads) THEN
         IF (DOF > 0) THEN
           Work(1:n) = ListGetReal(ValueList, LoadName(1:nlen) // ' im', n, NodeIndexes, gotIt)
         ELSE
           CALL ListGetRealArray(ValueList, LoadName(1:nlen) // ' im', WorkA, n, NodeIndexes, gotIt)
         END IF
         
         IF (GotIt) THEN
           DO j=1,n
             IF ( NodeIndexes(j) > SIZE(Perm) ) THEN
               CALL Warn('SetPointLoads','Node number too large!')
               CYCLE
             END IF             
             IF( NodeIndexes(j) == 0 ) THEN             
               CALL Warn('SetPointLoads','Node number is zero')
               CYCLE
             END IF
         
             k = Perm(NodeIndexes(j))
             IF ( k > 0 ) THEN
               IF (DOF > 0) THEN
                 k = NDOFs * (k-1) + DOF
                 A % RHS_im(k) = A % RHS_im(k) + Work(j) 
               ELSE
                 DO l=1,MIN( NDOFs, SIZE(WorkA,1) )
                   k1 = NDOFs * (k-1) + l
                   A % RHS_im(k1) = A % RHS_im(k1) + WorkA(l,1,j) 
                 END DO
               END IF
             END IF
           END DO
         END IF
       END IF

     END SUBROUTINE SetPointLoads
     
!------------------------------------------------------------------------------
   END SUBROUTINE SetNodalLoads
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> This subroutine seeks for nodes which are adjacent to the given target node
!> and then creates a couple which corresponds to a given torque. If the 
!> optional definition of the director vector d is given, the torque arm should 
!> ideally be parallel to d and the couple created does not have a d-component. 
!> This version may be more convenient when the torque comes from a dimensionally
!> reduced model over a thin body. Without specifying the director, this 
!> subroutine expects a 3-D geometry.
!
! TO DO: - The target nodes can now be defined only by their indices
!        - Add a way to find the director from the specification of a shell model.
!------------------------------------------------------------------------------
   SUBROUTINE SetCoupleLoads(Model, Perm, A, F, Dofs)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     TYPE(Model_t) :: Model                     !< The current model structure
     INTEGER, POINTER, INTENT(IN) :: Perm(:)    !< The permutation of the associated variable
     TYPE(Matrix_t), INTENT(INOUT) :: A         !< The coefficient matrix of the problem
     REAL(KIND=dp), POINTER, INTENT(INOUT) :: F(:) !< The RHS vector of the problem
     INTEGER, INTENT(IN) :: Dofs                !< The DOF count of the associated variable
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(ValueList_t), POINTER :: ValueList

     LOGICAL :: WithDirector
     LOGICAL :: Found, NoUpperNode, NoLowerNode
     
     INTEGER, ALLOCATABLE :: NearNodes(:) 
     INTEGER, POINTER :: NodeIndexes(:)
     INTEGER, POINTER :: Cols(:), Rows(:), Diag(:)
     INTEGER :: Row, TargetNode, TargetInd, BC, TargetCount
     INTEGER :: i, j, k, l, n, p
     INTEGER :: jx, lx, jy, ly, jz, lz 
     INTEGER :: intarray(1)

     REAL(KIND=dp), ALLOCATABLE :: NearCoordinates(:,:), AllDirectors(:,:), Work(:,:)
     REAL(KIND=dp) :: E(3,3)
     REAL(KIND=dp) :: Torque(3)  ! The torque vector with respect to the global frame
     REAL(KIND=dp) :: d(3)       ! Director at a solid-shell/plate interface    
     REAL(KIND=dp) :: ex(3), ey(3), ez(3)
     REAL(KIND=dp) :: e1(3), e2(3), e3(3)
     REAL(KIND=dp) :: T(3), Force(3), v(3)
     REAL(KIND=dp) :: M1, M2, F1, F2, F3
     REAL(KIND=dp) :: res_x, maxres_x, minres_x
     REAL(KIND=dp) :: res_y, maxres_y, minres_y
     REAL(KIND=dp) :: res_z, maxres_z, minres_z
     REAL(KIND=dp) :: rlower, rupper, FVal, MVal
!------------------------------------------------------------------------------
     IF (.NOT. ListCheckPrefixAnyBC(Model, 'Torque')) RETURN

     Mesh => Model % Solver % Mesh

     IF (.NOT. ASSOCIATED(A % InvPerm)) THEN
       ALLOCATE(A % InvPerm(A % NumberOfRows))
       DO i = 1,SIZE(Perm)
         IF (Perm(i) > 0) THEN
           A % InvPerm(Perm(i)) = i
         END IF
       END DO
     END IF

     ex = [1.0d0, 0.0d0, 0.0d0]
     ey = [0.0d0, 1.0d0, 0.0d0]
     ez = [0.0d0, 0.0d0, 1.0d0]
     E(:,1) = ex
     E(:,2) = ey
     E(:,3) = ez

     Diag   => A % Diag
     Rows   => A % Rows
     Cols   => A % Cols

     DO BC=1,Model % NumberOfBCs
       ValueList => Model % BCs(BC) % Values
       IF (.NOT.ListCheckPresent(ValueList, 'Torque 1') .AND. &
           .NOT.ListCheckPresent(ValueList, 'Torque 2') .AND. &
           .NOT.ListCheckPresent(ValueList, 'Torque 3')) CYCLE
       NodeIndexes => ListGetIntegerArray(ValueList, 'Target Nodes', UnfoundFatal=.TRUE.)

       TargetCount = SIZE(NodeIndexes)
       ALLOCATE(Work(3,TargetCount))
       Work(1,1:TargetCount) = ListGetReal(ValueList, 'Torque 1', TargetCount, NodeIndexes, Found)
       Work(2,1:TargetCount) = ListGetReal(ValueList, 'Torque 2', TargetCount, NodeIndexes, Found)
       Work(3,1:TargetCount) = ListGetReal(ValueList, 'Torque 3', TargetCount, NodeIndexes, Found)

       !
       ! Check whether the torque arm is given by the director vector. This option
       ! is not finalized yet. Here the director definition is sought from the BC
       ! definition, while the director might already be available from the specification 
       ! of a shell model.
       !
       IF (.NOT.ListCheckPresent(ValueList, 'Director 1') .AND. &
           .NOT.ListCheckPresent(ValueList, 'Director 2') .AND. &
           .NOT.ListCheckPresent(ValueList, 'Director 3')) THEN
         WithDirector = .FALSE.
       ELSE
         WithDirector = .TRUE.
         ALLOCATE(AllDirectors(3,TargetCount))
         AllDirectors(1,1:TargetCount) = ListGetReal(ValueList, 'Director 1', TargetCount, NodeIndexes, Found)
         AllDirectors(2,1:TargetCount) = ListGetReal(ValueList, 'Director 2', TargetCount, NodeIndexes, Found)
         AllDirectors(3,1:TargetCount) = ListGetReal(ValueList, 'Director 3', TargetCount, NodeIndexes, Found)
       END IF

       DO p=1,TargetCount
         TargetNode = NodeIndexes(p)
         TargetInd = Perm(NodeIndexes(p))
         IF (TargetInd == 0) CYCLE

         !------------------------------------------------------------------------------
         ! Find nodes which can potentially be used to make a representation of couple:
         !------------------------------------------------------------------------------
         Row = TargetInd * Dofs
         n = (Rows(Row+1)-1 - Rows(Row)-Dofs+1)/DOFs + 1
         ALLOCATE(NearNodes(n), NearCoordinates(3,n))

         k = 0
         DO i = Rows(Row)+Dofs-1, Rows(Row+1)-1, Dofs
           j = Cols(i)/Dofs
           k = k + 1
           NearNodes(k) = A % InvPerm(j)
         END DO
         ! PRINT *, 'POTENTIAL NODE CONNECTIONS:'
         ! print *, 'Nodes near target=', NearNodes(1:k)

         !
         ! The position vectors for the potential nodes where forces may be applied:
         !
         NearCoordinates(1,1:n) = Mesh % Nodes % x(NearNodes(1:n)) - Mesh % Nodes % x(TargetNode)
         NearCoordinates(2,1:n) = Mesh % Nodes % y(NearNodes(1:n)) - Mesh % Nodes % y(TargetNode)
         NearCoordinates(3,1:n) = Mesh % Nodes % z(NearNodes(1:n)) - Mesh % Nodes % z(TargetNode)


         IF (WithDirector) THEN
           !
           ! In this case the torque arm should ideally be parallel to the director vector d.
           ! Construct an orthonormal basis, with d giving the third basis vector.
           !
           d = AllDirectors(:,p)
           e3 = d/SQRT(DOT_PRODUCT(d,d))
           v(1:3) = ABS([DOT_PRODUCT(ex,e3), DOT_PRODUCT(ey,e3), DOT_PRODUCT(ez,e3)]) 
           intarray = MINLOC(v)
           k = intarray(1)
           v(1:3) = E(1:3,k)
           e1 = v - DOT_PRODUCT(v,e3)*e3
           e1 = e1/SQRT(DOT_PRODUCT(e1,e1))
           e2 = CrossProduct(e3,e1)
           !
           ! The torque is supposed to have no component in the direction of d, so remove it
           ! and also find the representation of the altered torque with respect to the local basis:
           !
           Torque = Work(:,p)
           v = DOT_PRODUCT(Torque,e3)*e3
           T = Torque - v
           M1 = DOT_PRODUCT(T,e1)
           M2 = DOT_PRODUCT(T,e2)

           !------------------------------------------------------------------------------
           ! Seek torque arms which are closest to be parallel to d:
           !------------------------------------------------------------------------------
           maxres_z = 0.0d0
           minres_z = 0.0d0
           jz = 0
           lz = 0
           DO i=1,n
             IF (NearNodes(i) == TargetNode) CYCLE
             res_z = DOT_PRODUCT(e3(:), NearCoordinates(:,i)) / &
                 SQRT(DOT_PRODUCT(NearCoordinates(:,i), NearCoordinates(:,i)))
             IF (res_z > 0.0d0) THEN
               !
               ! A near node is on +d side
               !
               IF (res_z > maxres_z) THEN
                 jz = NearNodes(i)
                 maxres_z = res_z
               END IF
             ELSE
               !
               ! A near node is on -d side
               !
               IF (res_z < minres_z) THEN
                 lz = NearNodes(i)
                 minres_z = res_z
               END IF
             END IF
           END DO

           !
           ! Calculate arm lengths with respect to the coordinate axis parallel to d:
           !
           NoUpperNode = .FALSE.
           NoLowerNode = .FALSE.
           IF (jz == 0 .OR. ABS(maxres_z) < AEPS) THEN
             NoUpperNode = .TRUE.
           ELSE
             rupper = DOT_PRODUCT(e3(:), [ Mesh % Nodes % x(jz) - Mesh % Nodes % x(TargetNode), &
                 Mesh % Nodes % y(jz) - Mesh % Nodes % y(TargetNode), &
                 Mesh % Nodes % z(jz) - Mesh % Nodes % z(TargetNode) ])
             ! print *, 'THE NODE ON +d SIDE = ', JZ
             ! print *, 'TORQUE ARM = ', rupper
           END IF

           IF (lz == 0 .OR. ABS(minres_z) < AEPS) THEN
             NoLowerNode = .TRUE.
           ELSE
             rlower = DOT_PRODUCT(-e3(:), [ Mesh % Nodes % x(lz) - Mesh % Nodes % x(TargetNode), &
                 Mesh % Nodes % y(lz) - Mesh % Nodes % y(TargetNode), &
                 Mesh % Nodes % z(lz) - Mesh % Nodes % z(TargetNode) ])
             ! print *, 'THE NODE ON -d SIDE = ', LZ
             ! print *, 'TORQUE ARM = ', rlower
           END IF

           IF (NoUpperNode .OR. NoLowerNode) THEN
             CALL Warn('SetCoupleLoads', 'A couple BC would need two nodes on opposite sides')
           ELSE
             !
             ! The torque generated from point loads as M1 * e1 + M2 * e2 = (r e3) x (f1 * e1 - f2 * e2) = 
             ! (r*f2)* e1 + (r*f1)* e2
             !
             F2 = M1/(rupper + rlower)
             F1 = M2/(rupper + rlower)
             Force = F1 * e1 - F2 * e2
             !
             ! Finally compute the components of force with respect to the global frame and
             ! add to the RHS: 
             !
             F1 = DOT_PRODUCT(Force,ex)
             F2 = DOT_PRODUCT(Force,ey)
             F3 = DOT_PRODUCT(Force,ez)

             k = Perm(jz)
             F((k-1)*Dofs+1) = F((k-1)*Dofs+1) + F1
             F((k-1)*Dofs+2) = F((k-1)*Dofs+2) + F2
             IF (Dofs > 2) F((k-1)*Dofs+3) = F((k-1)*Dofs+3) + F3
             k = Perm(lz)
             F((k-1)*Dofs+1) = F((k-1)*Dofs+1) - F1
             F((k-1)*Dofs+2) = F((k-1)*Dofs+2) - F2
             IF (Dofs > 2) F((k-1)*Dofs+3) = F((k-1)*Dofs+3) - F3
           END IF

         ELSE
           !------------------------------------------------------------------------------
           ! Seek torque arms which are closest to be parallel to the global coordinate
           ! axes: 
           !------------------------------------------------------------------------------
           maxres_x = 0.0d0
           minres_x = 0.0d0
           maxres_y = 0.0d0
           minres_y = 0.0d0
           maxres_z = 0.0d0
           minres_z = 0.0d0
           jx = 0
           lx = 0
           jy = 0
           ly = 0
           jz = 0
           lz = 0
           DO i=1,n
             IF (NearNodes(i) == TargetNode) CYCLE

             IF (ABS(Torque(3)) > AEPS) THEN
               res_x = DOT_PRODUCT(ex(:), NearCoordinates(:,i)) / &
                   SQRT(DOT_PRODUCT(NearCoordinates(:,i), NearCoordinates(:,i)))
               IF (res_x > 0.0d0) THEN
                 !
                 ! A near node is on +E_X side
                 !
                 IF (res_x > maxres_x) THEN
                   jx = NearNodes(i)
                   maxres_x = res_x
                 END IF
               ELSE
                 !
                 ! A near node is on -E_X side
                 !
                 IF (res_x < minres_x) THEN
                   lx = NearNodes(i)
                   minres_x = res_x
                 END IF
               END IF
             END IF

             IF (ABS(Torque(1)) > AEPS) THEN
               res_y = DOT_PRODUCT(ey(:), NearCoordinates(:,i)) / &
                   SQRT(DOT_PRODUCT(NearCoordinates(:,i), NearCoordinates(:,i)))
               IF (res_y > 0.0d0) THEN
                 !
                 ! A near node is on +E_Y side
                 !
                 IF (res_y > maxres_y) THEN
                   jy = NearNodes(i)
                   maxres_y = res_y
                 END IF
               ELSE
                 !
                 ! A near node is on -E_Y side
                 !
                 IF (res_y < minres_y) THEN
                   ly = NearNodes(i)
                   minres_y = res_y
                 END IF
               END IF
             END IF

             IF (ABS(Torque(2)) > AEPS) THEN
               res_z = DOT_PRODUCT(ez(:), NearCoordinates(:,i)) / &
                   SQRT(DOT_PRODUCT(NearCoordinates(:,i), NearCoordinates(:,i)))
               IF (res_z > 0.0d0) THEN
                 !
                 ! A near node is on +E_Z side
                 !
                 IF (res_z > maxres_z) THEN
                   jz = NearNodes(i)
                   maxres_z = res_z
                 END IF
               ELSE
                 !
                 ! A near node is on -E_Z side
                 !
                 IF (res_z < minres_z) THEN
                   lz = NearNodes(i)
                   minres_z = res_z
                 END IF
               END IF
             END IF
           END DO

           IF (ABS(Torque(1)) > AEPS) THEN
             !------------------------------------------------------------------------------
             ! Calculate arm lengths with respect to the Y-axis:
             !------------------------------------------------------------------------------
             NoUpperNode = .FALSE.
             NoLowerNode = .FALSE.
             IF (jy == 0) THEN
               NoUpperNode = .TRUE.
             ELSE
               rupper = DOT_PRODUCT(ey(:), [ Mesh % Nodes % x(jy) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(jy) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(jy) - Mesh % Nodes % z(TargetNode) ])
             END IF

             IF (ly == 0) THEN
               NoLowerNode = .TRUE.
             ELSE
               rlower = DOT_PRODUCT(-ey(:), [ Mesh % Nodes % x(ly) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(ly) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(ly) - Mesh % Nodes % z(TargetNode) ])
             END IF

             !------------------------------------------------------------------------------
             ! Finally, create a couple which tends to cause rotation about the X-axis 
             ! provided nodes on both sides have been identified
             !------------------------------------------------------------------------------
             IF (NoUpperNode .OR. NoLowerNode) THEN
               CALL Warn('SetCoupleLoads', 'A couple BC would need two nodes on opposite Y-sides')
             ELSE
               !
               ! The torque M_X E_X = (r E_Y) x (f E_Z), with the force f>0 applied on +E_Y side:
               !
               MVal = Torque(1)
               FVal = Mval/(rupper + rlower)
               k = Perm(jy)
               F((k-1)*Dofs+3) = F((k-1)*Dofs+3) + Fval
               k = Perm(ly)
               F((k-1)*Dofs+3) = F((k-1)*Dofs+3) - Fval
             END IF
           END IF

           IF (ABS(Torque(2)) > AEPS) THEN
             !
             ! Calculate arm lengths with respect to the Z-axis:
             !
             NoUpperNode = .FALSE.
             NoLowerNode = .FALSE.
             IF (jz == 0) THEN
               NoUpperNode = .TRUE.
             ELSE
               rupper = DOT_PRODUCT(ez(:), [ Mesh % Nodes % x(jz) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(jz) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(jz) - Mesh % Nodes % z(TargetNode) ])
             END IF

             IF (lz == 0) THEN
               NoLowerNode = .TRUE.
             ELSE
               rlower = DOT_PRODUCT(-ez(:), [ Mesh % Nodes % x(lz) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(lz) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(lz) - Mesh % Nodes % z(TargetNode) ])
             END IF

             IF (NoUpperNode .OR. NoLowerNode) THEN
               CALL Warn('SetCoupleLoads', 'A couple BC would need two nodes on opposite Z-sides')
             ELSE
               !
               ! The torque M_Y E_Y = (r E_Z) x (f E_X), with the force f>0 applied on +E_Z side:
               !
               MVal = Torque(2)
               FVal = Mval/(rupper + rlower)
               k = Perm(jz)
               F((k-1)*Dofs+1) = F((k-1)*Dofs+1) + Fval
               k = Perm(lz)
               F((k-1)*Dofs+1) = F((k-1)*Dofs+1) - Fval
             END IF
           END IF

           IF (ABS(Torque(3)) > AEPS) THEN
             !
             ! Calculate arm lengths with respect to the X-axis:
             !
             NoUpperNode = .FALSE.
             NoLowerNode = .FALSE.
             IF (jx == 0) THEN
               NoUpperNode = .TRUE.
             ELSE
               rupper = DOT_PRODUCT(ex(:), [ Mesh % Nodes % x(jx) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(jx) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(jx) - Mesh % Nodes % z(TargetNode) ])
             END IF

             IF (lx == 0) THEN
               NoLowerNode = .TRUE.
             ELSE
               rlower = DOT_PRODUCT(-ex(:), [ Mesh % Nodes % x(lx) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(lx) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(lx) - Mesh % Nodes % z(TargetNode) ])
             END IF

             IF (NoUpperNode .OR. NoLowerNode) THEN
               CALL Warn('SetCoupleLoads', 'A couple BC would need two nodes on opposite Y-sides')
             ELSE
               !
               ! The torque M_Z E_Z = (r E_X) x (f E_Y), with the force f>0 applied on +E_X side:
               !
               MVal = Torque(3)
               FVal = Mval/(rupper + rlower)
               k = Perm(jx)
               F((k-1)*Dofs+1) = F((k-1)*Dofs+1) - Fval
               k = Perm(lx)
               F((k-1)*Dofs+1) = F((k-1)*Dofs+1) + Fval
             END IF
           END IF
         END IF

         DEALLOCATE(NearNodes, NearCoordinates)
       END DO
       DEALLOCATE(Work)
       IF (WithDirector) DEALLOCATE(AllDirectors)
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE SetCoupleLoads
!------------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------------
  SUBROUTINE CommunicateDirichletBCs(A)
  !-------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A

     REAL(KIND=dp), ALLOCATABLE :: d_e(:,:), g_e(:)
     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), fneigh(:), ineigh(:)
     INTEGER :: i,j,k,l,n,nn,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)

     IF( ParEnv % PEs<=1 ) RETURN

     ALLOCATE( fneigh(ParEnv % PEs), ineigh(ParEnv % PEs) )

     nn = 0
     ineigh = 0
     DO i=0, ParEnv % PEs-1
       k = i+1
       IF(.NOT.ParEnv % Active(k) ) CYCLE
       IF(i==ParEnv % myPE) CYCLE
       IF(.NOT.ParEnv % IsNeighbour(k) ) CYCLE
       nn = nn + 1
       fneigh(nn) = k
       ineigh(k) = nn
     END DO

     n = COUNT(A % ConstrainedDOF .AND. A % ParallelInfo % GInterface)
     ALLOCATE( s_e(n, nn ), r_e(n) )
     ALLOCATE( d_e(n, nn ), g_e(n) )

     CALL CheckBuffer( nn*3*n )

     ii = 0
     DO i=1, A % NumberOfRows
       IF(A % ConstrainedDOF(i) .AND. A % ParallelInfo % GInterface(i) ) THEN
          DO j=1,SIZE(A % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = A % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k == ParEnv % MyPE ) CYCLE
            k = k + 1
            k = ineigh(k)
            IF ( k> 0) THEN
              ii(k) = ii(k) + 1
              d_e(ii(k),k) = A % DValues(i)
              s_e(ii(k),k) = A % ParallelInfo % GlobalDOFs(i)
            END IF
          END DO
       END IF
     END DO

     DO i=1, nn
       j = fneigh(i) 

       CALL MPI_BSEND( ii(i),1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD,ierr )
       IF( ii(i) > 0 ) THEN
         CALL MPI_BSEND( s_e(1:ii(i),i),ii(i),MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,ierr )
         CALL MPI_BSEND( d_e(1:ii(i),i),ii(i),MPI_DOUBLE_PRECISION,j-1,112,ELMER_COMM_WORLD,ierr )
       END IF
     END DO

     DO i=1, nn
       j = fneigh(i)
       CALL MPI_RECV( n,1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD, status,ierr )
       IF ( n>0 ) THEN
         IF( n>SIZE(r_e)) THEN
           DEALLOCATE(r_e,g_e)
           ALLOCATE(r_e(n),g_e(n))
         END IF

         CALL MPI_RECV( r_e,n,MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,status,ierr )
         CALL MPI_RECV( g_e,n,MPI_DOUBLE_PRECISION,j-1,112,ELMER_COMM_WORLD, status,ierr )
         DO j=1,n
           k = SearchNode( A % ParallelInfo, r_e(j), Order=A % ParallelInfo % Gorder )
           IF ( k>0 ) THEN
             IF(.NOT. A % ConstrainedDOF(k)) THEN
               A % Dvalues(k) = g_e(j)
               A % ConstrainedDOF(k) = .TRUE.
             END IF
           END IF
         END DO
       END IF
     END DO
     DEALLOCATE(s_e, r_e, d_e, g_e)
  !-------------------------------------------------------------------------------
  END SUBROUTINE CommunicateDirichletBCs
  !-------------------------------------------------------------------------------

  
!-------------------------------------------------------------------------------
  SUBROUTINE EnforceDirichletConditions( Solver, A, b, OffDiagonal ) 
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: b(:)
    LOGICAL, OPTIONAL :: OffDiagonal
    
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: ScaleSystem, DirichletComm, Found, NoDiag
    REAL(KIND=dp) :: dval, s
    INTEGER :: i,j,k,n,n2,n3
    CHARACTER(*), PARAMETER :: Caller = 'EnforceDirichletConditions'
    LOGICAL :: Parallel, DoDiagScale
    
    
    Params => Solver % Values

    IF(.NOT. ALLOCATED( A % ConstrainedDOF ) ) THEN
      CALL Info(Caller,'ConstrainedDOF not associated, returning...',Level=8)
      RETURN
    END IF

    Parallel = ( ParEnv % PEs > 1 ) .AND. ( .NOT. Solver % Mesh % SingleMesh ) 

    n = COUNT( A % ConstrainedDOF )
    IF( Parallel ) n = ParallelReduction( n )
      
    IF( n == 0 ) THEN
      CALL Info(Caller,'No Dirichlet conditions to enforce, exiting!',Level=10)
      RETURN
    ELSE
      CALL Info(Caller,'Enforcing total of '//I2S(n)//' Dirichlet conditions.',Level=10)
    END IF    
        
    ! Communicate the Dirichlet conditions for parallel cases since there may be orphans      
    IF ( Parallel ) THEN
      DirichletComm = ListGetLogical( CurrentModel % Simulation, 'Dirichlet Comm', Found)
      IF(.NOT. Found) DirichletComm = .TRUE.
      IF( DirichletComm) THEN
        CALL Info('EnforceDirichletConditions','Communicating Dirichlet conditions',Level=20)
        CALL CommunicateDirichletBCs(A)
      END IF
    END IF
    
    IF( PRESENT( OffDiagonal ) ) THEN
      NoDiag = OffDiagonal
    ELSE
      NoDiag = .FALSE.
    END IF

    IF( NoDiag ) THEN
      ScaleSystem = .FALSE.
    ELSE
      ScaleSystem = ListGetLogical(Params,'Linear System Dirichlet Scaling',Found)
      IF(.NOT.Found) THEN
        ScaleSystem = ListGetLogical(Params,'Linear System Scaling',Found)
        IF(.NOT.Found) ScaleSystem=.TRUE.
      END IF
    END IF

    DoDiagScale = .FALSE.
    IF( ScaleSystem ) THEN
      CALL Info(Caller,'Applying Dirichlet conditions using scaled diagonal',Level=8)
      CALL ScaleLinearSystem(Solver,A,b,ApplyScaling=.FALSE.) !,scalingStr='diagonal')
      DoDiagScale = ASSOCIATED(A % DiagScaling)
    END IF
    
    n2 = 0
    IF( Parallel ) THEN
      n = 0
      DO i=1,A % NumberOfRows
        IF( A % ConstrainedDOF(i) ) THEN
          IF( A % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
            n = n + 1
          ELSE
            n2 = n2 + 1
          END IF
        END IF
      END DO
      n = ParallelReduction( n )
      n2 = ParallelReduction( n2 ) 
    END IF
        
    ! Eliminate all entries in matrix that may be eliminated in one sweep
    ! If this is an offdiagonal entry this cannot be done.  
    IF ( A % Symmetric .AND. .NOT. NoDiag ) CALL CRS_ElimSymmDirichlet(A,b)
 
    
    DO k=1,A % NumberOfRows

      IF ( A % ConstrainedDOF(k) ) THEN
        
        dval = A % Dvalues(k) 

        s = 1.0_dp
        IF( ScaleSystem ) THEN
          SELECT CASE(A % ScalingMethod ) 
          CASE( 1 )
            s = A % DiagScaling(k)            
            IF( ABS(s) <= TINY(s) ) s = 1.0_dp
            s = 1._dp / s**2
          CASE( 2 )
            s = 1.0_dp / A % DiagScaling(k)
          CASE( 3 ) 
            s = A % AveScaling
          END SELECT
        END IF
          
        CALL ZeroRow(A, k)

        ! Off-diagonal entries for a block matrix are neglected since the code will
        ! also go through the diagonal entries where the r.h.s. target value will be set.
        IF(.NOT. NoDiag ) THEN
#if 0
          ! It is not clear whether we want to set the Dirichlet conditions only
          ! for the owner partition, or to all partitions. The latter may perhaps
          ! help when using partition-specific preconditioners...
          IF( Parallel ) THEN
            IF( SIZE( A % ParallelInfo % NeighbourList(k) % Neighbours) > 1 ) THEN
              IF( A % ParallelInfo % NeighbourList(k) % Neighbours(1) /= ParEnv % MyPe ) THEN
                b(k) = 0.0_dp
                CYCLE
              END IF
            END IF
          END IF
#endif          
          CALL SetMatrixElement(A,k,k,s)
          b(k) = s * dval
        END IF

      END IF
    END DO

    ! Deallocate scaling since otherwise it could be misused out of context
    IF (DoDiagScale) DEALLOCATE( A % DiagScaling ) 
        
    CALL Info(Caller,'Dirichlet conditions enforced for dofs: '//I2S(n), Level=6)
    IF(n2>0) CALL Info(Caller,'Dirichlet conditions shared count: '//I2S(n2), Level=12)
    
  END SUBROUTINE EnforceDirichletConditions
!-------------------------------------------------------------------------------
   

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
    TYPE(Element_t), POINTER :: Element
    INTEGER :: i,j,k,n,np,t,ierr,iter, proc
    LOGICAL :: GotIt, Found, Conditional
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Indexes(:)
    REAL(KIND=dp), ALLOCATABLE :: Condition(:)
    TYPE buff_t
      INTEGER, ALLOCATABLE :: buff(:)
    END TYPE buff_t
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER, TARGET :: pIndexes(12)
    INTEGER, POINTER :: nlist(:)
    TYPE(Buff_t), ALLOCATABLE, TARGET :: n_index(:)
    INTEGER, ALLOCATABLE :: n_count(:), gbuff(:)
    TYPE(Variable_t), POINTER :: DispVar
    LOGICAL :: pDisp
    CHARACTER(*), PARAMETER :: Caller = 'CheckNormalTangentialBoundary'
!------------------------------------------------------------------------------

    ! need an early initialization to average normals across partitions:
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

    pDisp = .FALSE.
    NULLIFY( DispVar )
    DO i=1,Model % NumberOfSolvers
      IF( ListGetLogical( Model % Solvers(i) % Values,'Use p Normals',GotIt ) ) THEN
        DispVar => Model % Solvers(i) % Variable
        pDisp = .TRUE.
        IF( SIZE( DispVar % Perm ) > n ) THEN
          n = SIZE( DispVar % Perm )        
        END IF
        EXIT
      END IF
    END DO
            
    IF ( ASSOCIATED( BoundaryReorder ) ) THEN
      IF ( SIZE(BoundaryReorder) < n ) DEALLOCATE( BoundaryReorder )
    END IF
    
    IF ( .NOT. ASSOCIATED( BoundaryReorder ) ) THEN
      CALL Info( Caller,'Allocating BoundaryOrder of size: '//I2S(n),Level=12)
      IF( pDisp ) THEN
        CALL Info(Caller,'Creating normals for p-dofs as well!',Level=12)
      END IF
      ALLOCATE( BoundaryReorder(n) )
    END IF
    
    BoundaryReorder = 0
    
!------------------------------------------------------------------------------
    DO t=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
                  Mesh % NumberOfBoundaryElements

      Element => Model % Elements(t)
      IF ( Element % TYPE % ElementCode == 101 )  CYCLE

      Indexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes      
      ALLOCATE( Condition(n)  )
      
      DO i=1,Model % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
          IF ( ListGetLogical( Model % BCs(i) % Values,VariableName, gotIt) ) THEN
            Found = ListGetLogical( Model % BCs(i) % Values, &
                TRIM(VariableName) // ' Rotate',gotIt)
            IF ( Found .OR. .NOT. GotIt ) THEN
              Condition(1:n) = ListGetReal( Model % BCs(i) % Values, &
                  TRIM(VariableName) // ' Condition', n, Indexes, Conditional )
              
              DO j=1,n
                IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE                
                k = Indexes(j)
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
          IF (.NOT.Mesh % ParallelInfo % GInterface(i) ) CYCLE

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
          IF (.NOT.Mesh % ParallelInfo % GInterface(i) ) CYCLE

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
                800, ELMER_COMM_WORLD, ierr )
           IF ( n_count(i)>0 ) &
             CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                 801, ELMER_COMM_WORLD, ierr )
        END IF
      END DO

      DO i=1,ParEnv % PEs
        IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff)

        IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
           CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                800, ELMER_COMM_WORLD, status, ierr )
           IF ( n>0 ) THEN
             ALLOCATE( gbuff(n) )
             proc = status(MPI_SOURCE)
             CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                 801, ELMER_COMM_WORLD, status, ierr )

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


    ! We add the normals for p-elements at 2nd stage to get higher indexes for them.
    ! This way the lower indexes should be the same as for non p-elements.
    ! Also these dofs do not to be communicated.
    !-------------------------------------------------------------------------------
    IF( pDisp ) THEN
      IF( ListCheckPresentAnyBC( Model,TRIM(VariableName) // ' Condition') ) THEN
        CALL Fatal(Caller,'Cannot deal with conditional n-t condition and p-elements')
      END IF      
      
      DO t=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
          Mesh % NumberOfBoundaryElements
        
        Element => Model % Elements(t)
        IF ( Element % TYPE % ElementCode < 200 )  CYCLE
        
        n = Element % TYPE % NumberOfNodes
        np = mGetElementDOFs(pIndexes,Element,USolver=DispVar % Solver)
        Indexes => pIndexes

        DO i=1,Model % NumberOfBCs
          IF ( Element % BoundaryInfo % Constraint == &
              Model % BCs(i) % Tag ) THEN
            IF ( ListGetLogical( Model % BCs(i) % Values,VariableName, gotIt) ) THEN
              Found = ListGetLogical( Model % BCs(i) % Values, &
                  TRIM(VariableName) // ' Rotate',gotIt)
              IF ( Found .OR. .NOT. GotIt ) THEN
                DO j=n+1,np
                  k = Indexes(j)
                  IF ( BoundaryReorder(k)==0 ) THEN
                    NumberOfBoundaryNodes = NumberOfBoundaryNodes + 1
                    BoundaryReorder(k) = NumberOfBoundaryNodes
                  END IF
                END DO
              END IF
            END IF
          END IF
        END DO
      END DO
    END IF

    CALL Info(Caller,'Number of normal-tangential dofs: '&
        //I2S(NumberOfBoundaryNodes),Level=10)
    
!------------------------------------------------------------------------------

    n = 0
    IF( ASSOCIATED( BoundaryNormals ) ) THEN
      n = SIZE( BoundaryNormals, 1)
    END IF

    IF( n > 0 .AND. NumberOfBoundaryNodes /= n ) THEN
      DEALLOCATE( BoundaryNormals, BoundaryTangent1, BoundaryTangent2 )
    END IF

    IF ( NumberOfBoundaryNodes == 0 ) THEN        
      DEALLOCATE( BoundaryReorder ) 
    ELSE
      IF( n /= NumberOfBoundaryNodes ) THEN
        ALLOCATE( BoundaryNormals(NumberOfBoundaryNodes,3), & 
            BoundaryTangent1(NumberOfBoundaryNodes,3), &
            BoundaryTangent2(NumberOfBoundaryNodes,3) )
      END IF

      BoundaryNormals  = 0.0_dp
      BoundaryTangent1 = 0.0_dp
      BoundaryTangent2 = 0.0_dp
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
    INTEGER :: NumberOfBoundaryNodes, dim
    REAL(KIND=dp), POINTER :: BoundaryNormals(:,:),BoundaryTangent1(:,:), &
        BoundaryTangent2(:,:)
    CHARACTER(LEN=*) :: VariableName
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: i,j,k,l,m,n,np,t,iBC,ierr,proc,i1,i2,k1,k2
    LOGICAL :: GotIt, Found, PeriodicNormals, Conditional
    REAL(KIND=dp) :: s,Bu,Bv,Nrm(3),Basis(32),DetJ
    INTEGER, POINTER :: Indexes(:)
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
    LOGICAL :: MassConsistent, LhsSystem, RotationalNormals
    LOGICAL, ALLOCATABLE :: LhsTangent(:),RhsTangent(:)
    INTEGER :: LhsConflicts, NormalConflicts, ConflictCount
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Origin(3),Axis(3)
    INTEGER, TARGET :: pIndexes(12)
    REAL(KIND=dp), POINTER :: Pwrk(:,:)
    LOGICAL :: GotOrigin,GotAxis,OneSidedNormals,pDisp
    LOGICAL :: NtBoss, AnyNtBoss, ThisBoss
    INTEGER :: NtBossCount
    LOGICAL, ALLOCATABLE :: NtBossTag(:)
    TYPE(Variable_t), POINTER :: DispVar
    CHARACTER(*), PARAMETER :: Caller = 'AverageBoundaryNormals'
    !------------------------------------------------------------------------------

    CALL Info(Caller,'Setting boundary normals for n-t conditions',Level=8)
    
    pDisp = .FALSE.
    NULLIFY( DispVar )
    DO i=1,Model % NumberOfSolvers
      IF( ListGetLogical( Model % Solvers(i) % Values,'Use p Normals',GotIt ) ) THEN
        DispVar => Model % Solvers(i) % Variable
        pDisp = .TRUE.
        EXIT
      END IF
    END DO

    Mesh => Model % Mesh
    ElementNodes % x => x
    ElementNodes % y => y
    ElementNodes % z => z


    ! Tag all nodes that have priority over conflicting normal-tangential BCs.
    AnyNtBoss = ListGetLogicalAnyBC( Model,'Normal-Tangential Priority')
    IF(AnyNtBoss) THEN
      ALLOCATE(NtBossTag(Mesh % NumberOfNodes) ) 
      NtBossTag = .FALSE.
      DO t=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
          Mesh % NumberOfBoundaryElements
        Element => Model % Elements(t)
        DO i=1,Model % NumberOfBCs
          IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
            ThisBoss = ListGetLogical( Model % BCs(i) % Values,'Normal-Tangential Priority', gotIt)
            IF( ThisBoss ) NtBossTag(Element % NodeIndexes) = .TRUE.
            EXIT
          END IF
        END DO
      END DO
      n = COUNT( NtBossTag )
      CALL Info(Caller,'Number of nodes on with normal-tangential priority is: '//I2S(n),Level=10)
      NtBossCount = 0
    END IF

    
    Mesh => Model % Mesh
    NrmVar => VariableGet( Mesh % Variables, 'Normals' )
    
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
      
      IF ( NumberOfBoundaryNodes>0 ) THEN
        ALLOCATE( n_comp(SIZE(BoundaryReorder)) )
        n_comp = 0

        BoundaryNormals = 0._dp
        ConflictCount = 0

        DO t=Model % NumberOfBulkElements + 1, Model % NumberOfBulkElements + &
                      Model % NumberOfBoundaryElements
          Element => Model % Elements(t)
          IF ( Element % TYPE  % ElementCode < 200 ) CYCLE

          n = Element % TYPE % NumberOfNodes
          IF(pDisp) THEN
            np = mGetElementDOFs(pIndexes,Element,USolver=DispVar % Solver)
            Indexes => pIndexes
          ELSE
            np = n
            Indexes => Element % NodeIndexes
          END IF

          ElementNodes % x(1:n) = Model % Nodes % x(Indexes(1:n))
          ElementNodes % y(1:n) = Model % Nodes % y(Indexes(1:n))
          ElementNodes % z(1:n) = Model % Nodes % z(Indexes(1:n))

          ALLOCATE(Condition(n))

          DO i=1,Model % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
              BC => Model % BCs(i) % Values

              IF ( ListGetLogical( BC, VariableName, gotIt) ) THEN
                ThisBoss = .FALSE.
                IF( AnyNtBoss ) THEN
                  ThisBoss = ListGetLogical(BC,'Normal-Tangential Priority', gotIt)
                END IF
                                
                Found = ListGetLogical( BC, TRIM(VariableName) // ' Rotate',gotIt)
                IF ( Found .OR. .NOT. Gotit ) THEN
                  MassConsistent = ListGetLogical( BC,'Mass Consistent Normals',gotIt)
                  RotationalNormals = ListGetLogical(BC,'Rotational Normals',gotIt)

                  IF( RotationalNormals ) THEN
                    Pwrk => ListGetConstRealArray(BC,'Normals Origin',GotOrigin )
                    IF( GotOrigin ) THEN
                      IF( SIZE(Pwrk,1) /= 3 .OR. SIZE(Pwrk,2) /= 1 ) THEN
                        CALL Fatal(Caller,'Size of > Normals Origin < should be 3!')
                      END IF
                      Origin = Pwrk(1:3,1)
                    END IF
                    Pwrk => ListGetConstRealArray(BC,'Normals Axis',GotAxis )
                    IF( GotAxis ) THEN
                      IF( SIZE(Pwrk,1) /= 3 .OR. SIZE(Pwrk,2) /= 1 ) THEN
                        CALL Fatal(Caller,'Size of > Normals Axis < should be 3!')
                      END IF
                      Axis = Pwrk(1:3,1)
                      ! Normalize axis is it should just be used for the direction
                      Axis = Axis / SQRT( SUM( Axis*Axis ) )
                    END IF
                  END IF
                  
                  Condition(1:n) = ListGetReal( BC,&
                       TRIM(VariableName) // ' Condition', n, Indexes, Conditional )

                  DO j=1,n
                    IF ( Conditional ) THEN
                      IF( Condition(j) < 0._dp ) CYCLE
                    END IF

                    IF( AnyNtBoss ) THEN
                      IF(NtBossTag(Indexes(j)) .AND. .NOT. ThisBoss) THEN
                        NtBossCount = NtBossCount + 1
                        CYCLE
                      END IF
                    END IF                    
                    
                    k = BoundaryReorder( Indexes(j) )
                    IF (k>0) THEN
                      nrm = 0._dp
                      IF (MassConsistent) THEN
                        IF(j>n) CYCLE
                        CALL IntegMassConsistent(j,n,nrm)
                      ELSE IF( RotationalNormals ) THEN
                        nrm(1) = ElementNodes % x(j)
                        nrm(2) = ElementNodes % y(j)
                        nrm(3) = ElementNodes % z(j)

                        !PRINT *,'nrm:',j,nrm
                        
                        IF( GotOrigin ) nrm = nrm - Origin
                        IF( GotAxis ) THEN
                          nrm = nrm - SUM( nrm * Axis ) * Axis
                        ELSE ! Default axis is (0,0,1)
                          nrm(3) = 0.0_dp
                        END IF

                        nrm = nrm / SQRT( SUM( nrm * nrm ) )
                      ELSE
                        Bu = Element % TYPE % NodeU(j)
                        Bv = Element % TYPE % NodeV(j)
                        nrm = NormalVector(Element,ElementNodes,Bu,Bv,.TRUE.)
                      END IF

                      l = n_comp(Indexes(j))
                      n_comp(Indexes(j)) = l + 1
                      IF( l > 0 ) THEN
                        IF( SUM( BoundaryNormals(k,:) * nrm ) < -EPSILON(s) ) THEN
                          ConflictCount = ConflictCount + 1
                          !CALL Warn(Caller,'Node '//I2S(Indexes(j))//' has conflicting normal directions!')
                        END IF
                      END IF

                      BoundaryNormals(k,:) = BoundaryNormals(k,:) + nrm
                    END IF
                  END DO
                END IF
              END IF
            END IF
          END DO
          DEALLOCATE(Condition)
        END DO

        IF(AnyNtBoss ) THEN
          CALL Info(Caller,'Number of priority nodes for normal-tangential dofs: '&
              //I2S(NtBossCount),Level=10)
        END IF
                
        IF( ConflictCount > 0 ) THEN
          CALL Info(Caller,'There are '//I2S(ConflictCount)//' conflicting normal directions!',Level=8)
        END IF
        
        ! Here we go through the periodic projectors and average the normals
        ! such that the normals are the same where the nodes are the same.
        !--------------------------------------------------------------------
        DO iBC=1,Model % NumberOfBCs
          Projector => Model % BCs(iBC) % PMatrix
          IF ( .NOT. ASSOCIATED( Projector ) ) CYCLE

          ! This is the legacy periodic projector.
          ! The mortars etc. should be treated differently. 
          IF( Projector % ProjectorType /= PROJECTOR_TYPE_NODAL ) CYCLE
          BC => Model % BCs(iBC) % Values
                              
          ! TODO: consistent normals, if rotations given:
          ! ---------------------------------------------
          Rot => ListGetConstRealArray(BC,'Periodic BC Rotate', Found )
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
          
          ! Here we are projecting with transpose of the projector which is not
          ! really exact generally, but is usually better than not considering the values
          ! at all!
          !-------------------------------------------------------------------------------
          OneSidedNormals = ListGetLogical(BC,'One Sided Normals',Found ) 
          IF(.NOT. OneSidedNormals ) THEN
            DO i=1,Projector % NumberOfRows
              k = BoundaryReorder(Projector % InvPerm(i))
              IF ( k <= 0 ) CYCLE
              DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                IF ( Projector % Cols(l) <= 0 ) CYCLE
                m = BoundaryReorder(Projector % Cols(l))
                IF ( m>0 ) THEN
                  BoundaryNormals(k,:) = BoundaryNormals(k,:) + &
                      Projector % Values(l) * BoundaryNormals(m,:)
                END IF
              END DO
            END DO
          END IF
            
          ! Ok, now we need to nullify the values so that we can apply the projector
          ! in the next sequence. This used to be done before without the upper part
          !--------------------------------------------------------------------------
          DO i=1,Projector % NumberOfRows
            k = BoundaryReorder(Projector % InvPerm(i))
            IF ( k <= 0 ) CYCLE
            DO l=Projector % Rows(i),Projector % Rows(i+1)-1
              IF ( Projector % Cols(l) <= 0 ) CYCLE
              m = BoundaryReorder(Projector % Cols(l))
              IF ( m>0 ) BoundaryNormals(m,:) = 0.0_dp
            END DO
          END DO
        END DO

        ! Here we use the values of the master side to have same values
        ! on the slave side as well. So this is hierarchical.
        !----------------------------------------------------------------
        DO iBC=1,Model % NumberOfBCs
           Projector => Model % BCs(iBC) % PMatrix
           IF ( .NOT. ASSOCIATED( Projector ) ) CYCLE
           IF( Projector % ProjectorType /= PROJECTOR_TYPE_NODAL ) CYCLE
          
           ! TODO: consistent normals, if rotations given:
           ! ---------------------------------------------
           BC => Model % BCs(iBC) % Values
           Rot => ListGetConstRealArray(BC,'Periodic BC Rotate', Found )
           IF ( Found .AND. ASSOCIATED(Rot) ) THEN
             IF ( ANY(Rot/=0) ) CYCLE
           END IF

           DO i=1,Projector % NumberOfRows
              k = BoundaryReorder(Projector % InvPerm(i))
              IF ( k <= 0 ) CYCLE
              DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                IF ( Projector % Cols(l) <= 0 ) CYCLE
                m = BoundaryReorder(Projector % Cols(l))
                IF ( m > 0 ) THEN
                  BoundaryNormals(m,:) = BoundaryNormals(m,:) + &
                      Projector % Values(l) * BoundaryNormals(k,:)
                END IF
              END DO
            END DO
        END DO
      END IF

      ! Communicate normals in parallel case so that they are consistent
      ! over the interfaces
      !-----------------------------------------------------------------
      IF (ParEnv % PEs>1 ) THEN
        ALLOCATE( n_count(ParEnv% PEs),n_index(ParEnv % PEs) )
        n_count = 0

        IF ( NumberOfBoundaryNodes>0 ) THEN
          DO i=1,Mesh % NumberOfNodes
            IF (BoundaryReorder(i)<=0 .OR. n_comp(i)<=0 ) CYCLE
            IF (.NOT.Mesh % ParallelInfo % GInterface(i) ) CYCLE
  
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
            IF (.NOT.Mesh % ParallelInfo % GInterface(i) ) CYCLE

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
                900, ELMER_COMM_WORLD, ierr )
            IF ( n_count(i)>0 ) THEN
              CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                  901, ELMER_COMM_WORLD, ierr )
              CALL MPI_BSEND( n_index(i) % normals, 3*n_count(i), MPI_DOUBLE_PRECISION, &
                    i-1,  902, ELMER_COMM_WORLD, ierr )
            END IF
          END IF
        END DO
        DO i=1,ParEnv % PEs
          IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff, n_index(i) % Normals)

          IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
             CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                    900, ELMER_COMM_WORLD, status, ierr )
             IF ( n>0 ) THEN
               proc = status(MPI_SOURCE)
               ALLOCATE( gbuff(n), nbuff(3*n) )
               CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                   901, ELMER_COMM_WORLD, status, ierr )

               CALL MPI_RECV( nbuff, 3*n, MPI_DOUBLE_PRECISION, proc, &
                    902, ELMER_COMM_WORLD, status, ierr )

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
            CALL Warn(Caller,'BC '//I2S(i)//' is both N-T master and slave!')
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
          Indexes => Element % NodeIndexes
          
          DO i=1,Model % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
              IF( NtMasterBC(i) ) LhsTangent( Indexes ) = .TRUE.
              IF( NtSlaveBC(i) ) RhsTangent( Indexes ) = .TRUE.
              EXIT
            END IF
          END DO
        END DO

        LhsConflicts = COUNT( LhsTangent .AND. RhsTangent )
        IF( LhsConflicts > 0 ) THEN
          CALL Warn(Caller,&
              'There are '//I2S(LhsConflicts)//' nodes that could be both rhs and lhs!')
        END IF
      END IF
            
      
      ! Normalize the normals and compute the tangent directions.
      !----------------------------------------------------------
      DO i=1,Model % NumberOfNodes
        k = BoundaryReorder(i) 
        IF ( k > 0 ) THEN
          s = SQRT( SUM( BoundaryNormals(k,:)**2 ) )
          IF ( s > TINY(s) ) THEN
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
          ELSE
            CALL Warn(Caller,'Suspiciously small normal for node: '//I2S(i))
          END IF
        END IF        
      END DO
      

      ! Inherit the normal direction for 2nd order p-elements from the nodes.
      !----------------------------------------------------------------------
      IF( pDisp ) THEN
        DO t=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
            Mesh % NumberOfBoundaryElements
          
          Element => Model % Elements(t)
          IF ( Element % TYPE % ElementCode < 200 )  CYCLE
          
          n = Element % TYPE % NumberOfNodes
          np = mGetElementDOFs(pIndexes,Element,USolver=DispVar % Solver)
          Indexes => pIndexes

          IF(ALL(BoundaryReorder(Indexes(1:n)) == 0 ) ) CYCLE
          
          DO i=n+1,np
            i1 = i-n
            i2 = i+1-n
            IF(i2>n) i2=1

            k1 = BoundaryReorder(Indexes(i1)) 
            IF(k1==0) CYCLE

            k2 = BoundaryReorder(Indexes(i2))
            IF(k2==0) CYCLE
            
            k = BoundaryReorder(Indexes(i)) 

            BoundaryNormals(k,:) = ( BoundaryNormals(k1,:) + BoundaryNormals(k2,:) ) / 2

            ! Even though the two normals have unit length their mean may not have unit length. 
            s = SQRT( SUM( BoundaryNormals(k,:)**2 ) )
            IF ( s > TINY(s) ) THEN
              BoundaryNormals(k,:) = BoundaryNormals(k,:) / s
            ELSE
              CALL Warn(Caller,'Starnegly small normal for dofs: '//I2S(i))
            END IF
              
            IF( dim > 2 ) THEN
              BoundaryTangent1(k,:) = ( BoundaryTangent1(k1,:) + BoundaryTangent1(k2,:) ) / 2
              BoundaryTangent2(k,:) = ( BoundaryTangent2(k1,:) + BoundaryTangent2(k2,:) ) / 2

              s = SQRT( SUM( BoundaryTangent1(k,:)**2 ) )
              IF ( s > TINY(s) ) BoundaryTangent1(k,:) = BoundaryTangent1(k,:) / s

              s = SQRT( SUM( BoundaryTangent2(k,:)**2 ) )
              IF ( s > TINY(s) ) BoundaryTangent2(k,:) = BoundaryTangent2(k,:) / s
            END IF
          END DO
        END DO
      END IF

      ! Save the normals and tangents as fields if requested. 
      !----------------------------------------------------------
      IF( ListGetLogical( Model % Simulation,'Save Averaged Normals',Found ) ) THEN
        CALL Info(Caller,'Saving averaged boundary normals to variable: Averaged Normals')
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

    IF( InfoActive(25) ) THEN
      DO i=1,3        
        CALL VectorValuesRange(BoundaryNormals(:,i),NumberOfBoundaryNodes,'Normal '//I2S(i))
      END DO
      IF( dim > 2 ) THEN
        DO i=1,3        
          CALL VectorValuesRange(BoundaryTangent1(:,i),NumberOfBoundaryNodes,'Tangent1 '//I2S(i))
          CALL VectorValuesRange(BoundaryTangent2(:,i),NumberOfBoundaryNodes,'Tangent2 '//I2S(i))
        END DO
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
     TYPE(Solver_t), TARGET :: Solver  !< Solver to be initialized.
!------------------------------------------------------------------------------
     CHARACTER(:), ALLOCATABLE :: Method
     LOGICAL :: GotIt, ExtrapolateInTime
     INTEGER :: i, Order,ndofs
     REAL(KIND=dp), POINTER CONTIG :: SaveValues(:)
     TYPE(Matrix_t), POINTER :: A
     TYPE(Variable_t), POINTER :: Var
     TYPE(Solver_t), POINTER :: pSolver
     REAL(KIND=dp) :: eFact
     REAL(KIND=dp), ALLOCATABLE :: NewVals(:)
     
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
         IF (Solver % TimeOrder > 1) THEN
           Method = 'bossak'
         ELSE
           CALL Warn( 'InitializeTimestep', &
               'Timestepping method defaulted to IMPLICIT EULER' )

           Solver % Beta = 1.0_dp
           Method = 'implicit euler'
         END IF
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

         CASE('adams-bashforth')
           Solver % Beta = 0.0d0

         CASE('adams-moulton')
           Solver % Beta = 1.0d0

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

         CASE('generalized-alpha','bossak')

         CASE DEFAULT 
           WRITE( Message, * ) 'Unknown timestepping method: ',Method
           CALL Fatal( 'InitializeTimestep', Message )
       END SELECT

     END IF

     ndofs = Solver % Matrix % NumberOfRows
     Var => Solver % Variable

     eFact = 0.0_dp
     ExtrapolateInTime = ListGetLogical(CurrentModel % Simulation, &
              'Timestep extrapolation', GotIt  )
     IF( ExtrapolateInTime ) THEN
       eFact = ListGetCReal(CurrentModel % Simulation, 'Timestep extrapolation factor', GotIt )
       IF(.NOT.GotIt) eFact = 1._dp
     END IF

     IF ( Solver % TimeOrder>1 .OR. method/='bdf') THEN
       IF( Solver % TimeOrder == 2 ) THEN         
         Solver % Beta = -1 ! use this to select method later (for now)
         SELECT CASE(Method)
         CASE('generalized-alpha')
         Solver % Beta = ListGetConstReal(Solver % Values, 'Generalized-alpha rinf', GotIt)
         IF ( .NOT. GotIt ) &
           Solver % Beta = ListGetConstReal( CurrentModel % Simulation, &
                         'Generalized-alpha rinf', GotIt )
         IF( .NOT.GotIt ) Solver % Beta = 0.75_dp

         CASE('bossak')
           Solver % Alpha = ListGetConstReal(Solver % Values, 'Bossak Alpha', GotIt)
           IF ( .NOT. GotIt ) &
             Solver % Alpha = ListGetConstReal(CurrentModel % Simulation,'Bossak Alpha',GotIt)
           IF(.NOT.GotIt) Solver % Alpha = -0.05_dp

         CASE DEFAULT
           Solver % Alpha = ListGetConstReal(Solver % Values, 'Bossak Alpha', GotIt)
           IF ( .NOT. GotIt ) &
             Solver % Alpha = ListGetConstReal(CurrentModel % Simulation,'Bossak Alpha',GotIt)
           IF(.NOT.GotIt) Solver % Alpha = -0.05_dp
         END SELECT
       ELSE
         IF (Solver % DoneTime==1 .AND.Solver % Beta/= 0.0_dp) Solver % Beta = 1.0_dp
       END IF
       
       SELECT CASE( Solver % TimeOrder )
       CASE(1)
         Order = MIN(Solver % DoneTime, Solver % Order)
         DO i=Order, 2, -1
           Var % PrevValues(:,i) = Var % PrevValues(:,i-1)
         END DO

         IF( ExtrapolateInTime ) THEN
           NewVals = (1+eFact)*Var % Values-eFact*Var % PrevValues(:,1)
           Var % PrevValues(:,1) = Var % Values 
           Var % Values = NewVals
         ELSE
           Var % PrevValues(:,1) = Var % Values
         END IF
           
         Solver % Matrix % Force(:,2) = Solver % Matrix % Force(:,1)
         
       CASE(2)
         Var % PrevValues(:,3) = Var % Values
         Var % PrevValues(:,4) = Var % PrevValues(:,1)
         Var % PrevValues(:,5) = Var % PrevValues(:,2)
         Var % PrevValues(:,7) = Var % PrevValues(:,6)
       END SELECT
     ELSE
       Order = MIN(Solver % DoneTime, Solver % Order)
       DO i=Order, 2, -1
         Var % PrevValues(:,i) = Var % PrevValues(:,i-1)
       END DO

       IF( ExtrapolateInTime ) THEN
         NewVals = (1+eFact)*Var % Values-eFact*Var % PrevValues(:,1)
         Var % PrevValues(:,1) = Var % Values 
         Var % Values = NewVals
       ELSE
         Var % PrevValues(:,1) = Var % Values
       END IF
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
         CALL MatrixVectorMultiply( A, Var % Values, A % BulkResidual )
         A % Values => SaveValues
         A % BulkResidual = A % BulkResidual - A % BulkRhs
       END IF
     END IF


     ! Advance also the exported variables if they happen to be time-dependent
     ! They only have normal prevvalues, when writing this always 2. 
     BLOCK
       INTEGER :: VarNo,n
       LOGICAL :: Found
       CHARACTER(LEN=MAX_NAME_LEN) :: str, var_name
       
       VarNo =0      
       DO WHILE( .TRUE. )
         VarNo = VarNo + 1         
         str = ComponentName( 'exported variable', VarNo )    
         
         var_name = ListGetString( Solver % Values, str, Found )    
         IF(.NOT. Found) EXIT
         
         CALL VariableNameParser( var_name ) 
         
         Var => VariableGet( Solver % Mesh % Variables, Var_name )
         IF( .NOT. ASSOCIATED(Var)) CYCLE
         IF( .NOT. ASSOCIATED(Var % PrevValues) ) CYCLE

         n = SIZE( Var % PrevValues,2 )
         DO i=n,2,-1
           Var % PrevValues(:,i) = Var % PrevValues(:,i-1)
         END DO
         Var % PrevValues(:,1) = Var % Values
       END DO
     END BLOCK

     ! This is after the timestep initialization since then we can recommunicate the cyclic
     ! solution for PrevValues in time-parallel simulation.
     IF( ListGetLogical( Solver % Values,'Store Cyclic System',GotIt ) ) THEN      
       pSolver => Solver
       CALL StoreCyclicSolution(pSolver)
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
    INTEGER :: Order
    LOGICAL :: Found
    CHARACTER(:), ALLOCATABLE :: Method, Simulation
!------------------------------------------------------------------------------

    IF ( Solver % Matrix % FORMAT == MATRIX_LIST ) THEN
      CALL List_toCRSMatrix(Solver % Matrix)
    END IF

    Simulation = ListGetString( CurrentModel % Simulation, 'Simulation Type' )
    IF ( Simulation == 'transient' ) THEN
      Method = ListGetString( Solver % Values, 'Timestepping Method' )
      Order = MIN(Solver % DoneTime, Solver % Order)

      IF ( Order <= 0 .OR. Solver % TimeOrder /= 1 .OR. Method=='bdf' ) RETURN

      IF ( Solver % Beta /= 0.0_dp ) THEN
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
    CHARACTER(:), ALLOCATABLE :: tmpname
    INTEGER :: i
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var,Var1, PrimVar
!------------------------------------------------------------------------------
    Mesh => TopMesh

    PrimVar => VariableGet( PrimaryMesh % Variables, Name, ThisOnly=.TRUE.)
    IF ( .NOT.ASSOCIATED( PrimVar) ) RETURN

    DO WHILE( ASSOCIATED(Mesh) )
      ! Make the same variable invalid in all other meshes.
      IF ( .NOT.ASSOCIATED( PrimaryMesh, Mesh) ) THEN
        Var => VariableGet( Mesh % Variables, Name, ThisOnly=.TRUE.)
        IF ( ASSOCIATED( Var ) ) THEN
          Var % Valid = .FALSE.
          Var % PrimaryMesh => PrimaryMesh
        END IF

        IF ( PrimVar % DOFs > 1 ) THEN
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
      Mesh => Mesh % Next
    END DO 

    ! Tell that values have changed in the primary mesh.
    ! Interpolation can then be activated if we request the same variable in the
    ! other meshes. 
    PrimVar % ValuesChanged = .TRUE.
    IF ( PrimVar % DOFs > 1 ) THEN
      DO i=1,PrimVar % DOFs
        tmpname = ComponentName( Name, i )
        Var => VariableGet( PrimaryMesh % Variables, tmpname, .TRUE. )
        IF ( ASSOCIATED(Var) ) Var % ValuesChanged = .TRUE.
      END DO
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
     TYPE(NormalTangential_t), POINTER :: NT
!------------------------------------------------------------------------------

     NT => CurrentModel % Solver % NormalTangential
     
     IF ( NT % NormalTangentialNOFNodes <= 0 ) RETURN
     
     IF( NodeNumber > SIZE( NT % BoundaryReorder ) ) THEN
       CALL Fatal('RotateNTSystem',&
           'Index '//I2S(NodeNumber)//' beyond BoundaryReorder size '//I2S(SIZE(NT % BoundaryReorder)))
     END IF
     
     k = NT % BoundaryReorder(NodeNumber)
     IF ( k <= 0 ) RETURN

     dim = CoordinateSystemDimension()
     IF ( dim < 3 ) THEN
       Bu = Vec(1)
       Bv = Vec(2)
       Vec(1) =  NT % BoundaryNormals(k,1)*Bu + NT % BoundaryNormals(k,2)*Bv
       Vec(2) = -NT % BoundaryNormals(k,2)*Bu + NT % BoundaryNormals(k,1)*Bv
     ELSE
       Bu = Vec(1)
       Bv = Vec(2)
       Bw = Vec(3)

       RM(:,1) = NT % BoundaryNormals(k,:)
       RM(:,2) = NT % BoundaryTangent1(k,:)
       RM(:,3) = NT % BoundaryTangent2(k,:)

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
    TYPE(NormalTangential_t), POINTER :: NT
!------------------------------------------------------------------------------

    NT => CurrentModel % Solver % NormalTangential

    IF ( NT % NormalTangentialNOFNodes <= 0 ) RETURN

    dim = CoordinateSystemDimension()    
    IF( ndofs < dim ) RETURN

    
    DO i=1,SIZE(NT % BoundaryReorder)
       k = NT % BoundaryReorder(i)
       IF ( k <= 0 ) CYCLE
       j = Perm(i)
       IF ( j <= 0 ) CYCLE

       IF ( dim < 3 ) THEN
          Bu = Solution(NDOFs*(j-1)+1)
          Bv = Solution(NDOFs*(j-1)+2)

          Solution(NDOFs*(j-1)+1) = NT % BoundaryNormals(k,1)*Bu + NT % BoundaryNormals(k,2)*Bv
          Solution(NDOFs*(j-1)+2) = -NT % BoundaryNormals(k,2)*Bu + NT % BoundaryNormals(k,1)*Bv

       ELSE
          Bu = Solution(NDOFs*(j-1)+1)
          Bv = Solution(NDOFs*(j-1)+2)
          Bw = Solution(NDOFs*(j-1)+3)
 
          RM(:,1) = NT % BoundaryNormals(k,:)
          RM(:,2) = NT % BoundaryTangent1(k,:)
          RM(:,3) = NT % BoundaryTangent2(k,:)

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
     TYPE(NormalTangential_t), POINTER :: NT
!------------------------------------------------------------------------------

     NT => CurrentModel % Solver % NormalTangential

     IF ( NT % NormalTangentialNOFNodes <= 0 ) RETURN

     dim = CoordinateSystemDimension()
     IF ( ndofs < dim ) RETURN

     
     DO i=1,SIZE(NT % BoundaryReorder)
       k = NT % BoundaryReorder(i)
       IF ( k <= 0 ) CYCLE
       j = Perm(i)
       IF ( j <= 0 ) CYCLE

       IF ( dim < 3 ) THEN
         Bu = Solution(NDOFs*(j-1)+1)
         Bv = Solution(NDOFs*(j-1)+2)

         Solution(NDOFs*(j-1)+1) = NT % BoundaryNormals(k,1) * Bu - &
                         NT % BoundaryNormals(k,2) * Bv

         Solution(NDOFs*(j-1)+2) = NT % BoundaryNormals(k,2) * Bu + &
                         NT % BoundaryNormals(k,1) * Bv
       ELSE
         Bu = Solution(NDOFs*(j-1)+1)
         Bv = Solution(NDOFs*(j-1)+2)
         Bw = Solution(NDOFs*(j-1)+3)

         RM(1,:) = NT % BoundaryNormals(k,:)
         RM(2,:) = NT % BoundaryTangent1(k,:)
         RM(3,:) = NT % BoundaryTangent2(k,:)

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
    TYPE(NormalTangential_t), POINTER :: NT
!------------------------------------------------------------------------------

    A = 0._dp
    k = 0 

    NT => CurrentModel % Solver % NormalTangential
    
    IF (NT % NormalTangentialNOFNodes > 0) THEN
      k = NT % BoundaryReorder(n)
    END IF
      
    IF (k > 0) THEN
      Rotated = .TRUE.
      dim = CoordinateSystemDimension()
      IF (dim==2) THEN
        A(1,1) = NT % BoundaryNormals(k,1)
        A(1,2) = -NT % BoundaryNormals(k,2)
        A(2,1) = NT % BoundaryNormals(k,2)
        A(2,2) = NT % BoundaryNormals(k,1)
        A(3,3) = 1._dp
      ELSE
        A(:,1) = NT % BoundaryNormals(k,:)
        A(:,2) = NT % BoundaryTangent1(k,:)
        A(:,3) = NT % BoundaryTangent2(k,:)
      END IF
    ELSE
      Rotated = .FALSE.
      A(1,1)=1._dp
      A(2,2)=1._dp
      A(3,3)=1._dp            
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
    
    INTEGER :: NormDim, NormDofs, Dofs,i,j,k,n,nn,totn,PermStart
    INTEGER, POINTER :: NormComponents(:)
    INTEGER, ALLOCATABLE :: iPerm(:)
    REAL(KIND=dp) :: Norm, nscale, val
    LOGICAL :: Stat, ComponentsAllocated, ConsistentNorm
    REAL(KIND=dp), POINTER :: x(:)
    REAL(KIND=dp), ALLOCATABLE, TARGET :: y(:)
    LOGICAL :: Parallel
    LOGICAL, ALLOCATABLE :: PassiveDof(:)
    INTEGER, POINTER :: Perm(:)
    
    CALL Info('ComputeNorm','Computing norm of solution',Level=10)

    IF(PRESENT(values)) THEN
      x => values
    ELSE
      x => Solver % Variable % Values
    END IF
    
    Parallel = Solver % Parallel
    
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

    ALLOCATE( PassiveDof(0:Dofs-1) )
    IF( NormDofs < Dofs ) THEN
      PassiveDof = .TRUE.
      DO i=1,NormDofs
        PassiveDof(NormComponents(i)-1) = .FALSE.        
      END DO
    ELSE
      PassiveDof = .FALSE.
    END IF
          
    
    n = nin
    totn = 0
    
    IF( Parallel ) THEN
      ConsistentNorm = ListGetLogical(Solver % Values,'Nonlinear System Consistent Norm',Stat)
      IF (ConsistentNorm) CALL Info('ComputeNorm','Using consistent norm in parallel',Level=10)
    ELSE
      ConsistentNorm = .FALSE.
    END IF

    Perm => Solver % Variable % Perm
    
    PermStart = ListGetInteger(Solver % Values,'Norm Permutation',Stat)
    IF ( Stat .AND. PermStart > 1 ) THEN
      ALLOCATE(iPerm(SIZE(Perm))); iPerm=0
      n = 0
      DO i=PermStart,SIZE(iPerm)
        IF ( Perm(i)>0 ) THEN
          n = n + 1
          iPerm(n) = Perm(i)
        END IF
      END DO
      
      ALLOCATE(y(n))
      y = x(iPerm(1:n))
      x => y
      DEALLOCATE(iPerm)
    END IF

    IF( ListGetLogical( Solver % Values,'Nonlinear System Nodal Norm', Stat ) ) THEN
      nn = Solver % Mesh % NumberOfNodes
      i = MAXVAL(Perm(1:nn))
      j = MINVAL(Perm(nn+1:SIZE(Perm)))
      IF(j==0) THEN
        CONTINUE
      ELSE IF(j>i) THEN      
        n = Dofs * i
        CALL Info('ComputeNorm','Considering only the nodal entries in norm computation!',Level=7)
      ELSE
        CALL Warn('ComputeNorm','Nodal norm is only available if all the p-dofs follow the last nodal dof!')
      END IF
    END IF
    

    IF( ConsistentNorm ) THEN
      ! In consistent norm we have to skip the dofs not owned by the partition in order
      ! to count each dof only once. 
      Norm = 0.0_dp
      
      SELECT CASE(NormDim)

      CASE(0) 
        DO j=1,n
          IF(PassiveDof(MODULO(j-1,Dofs))) CYCLE
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = MAX( Norm, ABS( val ) )
          totn = totn + 1
        END DO

      CASE(1)
        DO j=1,n
          IF(PassiveDof(MODULO(j-1,Dofs))) CYCLE
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = Norm + ABS(val)
          totn = totn + 1
        END DO

      CASE(2)          
        DO j=1,n
          IF(PassiveDof(MODULO(j-1,Dofs))) CYCLE
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = Norm + val**2
          totn = totn + 1
        END DO

      CASE DEFAULT
        DO j=1,n
          IF(PassiveDof(MODULO(j-1,Dofs))) CYCLE
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = Norm + val**NormDim 
          totn = totn + 1
        END DO
      END SELECT
      
      totn = ParallelReduction(totn) 
      IF(totn == 0) GOTO 10

      nscale = 1.0_dp * totn
      
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
    
    ELSE IF( NormDofs < Dofs ) THEN
      Norm = 0.0_dp
      totn = ParallelReduction(n) 
      IF(totn == 0) GOTO 10

      nscale = NormDOFs*totn/(1._dp*DOFs)

      SELECT CASE(NormDim)
      CASE(0)
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = MAX(Norm, MAXVAL( ABS(x(j::Dofs))) )
        END DO
        IF( Parallel ) THEN
          Norm = ParallelReduction(Norm,2)
        END IF
      CASE(1)
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = Norm + SUM( ABS(x(j::Dofs)) )
        END DO
        IF( Parallel ) THEN
          Norm = ParallelReduction(Norm)/nscale
        ELSE
          Norm = Norm/nscale
        END IF
      CASE(2)
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = Norm + SUM( x(j::Dofs)**2 )
        END DO
        IF( Parallel ) THEN
          Norm = SQRT(ParallelReduction(Norm)/nscale)
        ELSE
          Norm = SQRT(Norm/nscale)
        END IF
      CASE DEFAULT
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = Norm + SUM( x(j::Dofs)**NormDim )
        END DO
        IF( Parallel ) THEN
          Norm = (ParallelReduction(Norm)/nscale)**(1.0d0/NormDim)
        ELSE
          Norm = (Norm/nscale)**(1.0d0/NormDim)
        END IF
      END SELECT
            
    ELSE IF( Parallel ) THEN
      totn = ParallelReduction(n)
      IF (totn == 0) THEN
        CALL Warn('ComputeNorm','Requested norm of a variable with no Dofs')
        Norm = 0.0_dp
      ELSE
        nscale = 1.0_dp * totn        
        val = 0.0_dp
        SELECT CASE(NormDim)
        CASE(0)
          IF (n>0) val = MAXVAL(ABS(x(1:n)))
          Norm = ParallelReduction(val,2)
        CASE(1)
          IF (n>0) val = SUM(ABS(x(1:n)))
          Norm = ParallelReduction(val)/nscale
        CASE(2)
          IF (n>0) val = SUM(x(1:n)**2)
          Norm = SQRT(ParallelReduction(val)/nscale)
        CASE DEFAULT
          IF (n>0) val = SUM(x(1:n)**NormDim)
          Norm = (ParallelReduction(val)/nscale)**(1.0d0/NormDim)
        END SELECT
      END IF

    ELSE
      Norm = 0.0_dp
      IF(n==0) GOTO 10 
      
      SELECT CASE(NormDim)
      CASE(0)
        Norm = MAXVAL(ABS(x(1:n)))
      CASE(1)
        Norm = SUM(ABS(x(1:n)))/n
      CASE(2)
        Norm = SQRT(SUM(x(1:n)**2)/n)
      CASE DEFAULT
        Norm = (SUM((x(1:n)**NormDim)/n))**(1.0_dp/NormDim)
      END SELECT
    END IF
    
10  IF( ComponentsAllocated ) THEN
      DEALLOCATE( NormComponents ) 
    END IF
!------------------------------------------------------------------------------
  END FUNCTION ComputeNorm
!------------------------------------------------------------------------------


  SUBROUTINE UpdateDependentObjects( Solver, SteadyState )

    TYPE(Solver_t), TARGET :: Solver
    LOGICAL :: SteadyState

    TYPE(ValueList_t), POINTER :: SolverParams
    LOGICAL :: Found, DoIt
    REAL(KIND=dp) :: dt
    TYPE(Variable_t), POINTER :: dtVar, VeloVar, pVar
    CHARACTER(:), ALLOCATABLE :: str
    INTEGER, POINTER :: UpdateComponents(:)
    INTEGER :: k
    CHARACTER(*), PARAMETER :: Caller = 'UpdateDependentObjects'
  
    SolverParams => Solver % Values
    
    IF( SteadyState ) THEN
      CALL Info(Caller,'Updating objects depending on primary field in steady state',Level=20)
    ELSE
      CALL Info(Caller,'Updating objects depending on primary field in nonlinear system',Level=20)
    END IF
    
    
    ! The update of exported variables on nonlinear or steady state level.
    ! In nonlinear level the nonlinear iteration may depend on the updated values.
    ! Steady-state level is often sufficient if the dependendence is on some other solver.
    !-----------------------------------------------------------------------------------------    
    IF( SteadyState ) THEN
      DoIt = ListGetLogical( SolverParams,&
          'Update Exported Variables', Found )
    ELSE        
      DoIt = ListGetLogical( SolverParams,&
          'Nonlinear Update Exported Variables',Found )
    END IF
    IF( DoIt ) THEN
      CALL Info(Caller,'Updating exported variables',Level=20)
      CALL UpdateExportedVariables( Solver )	
    END IF

    
    ! Update components that depende on the solution of the solver.
    ! Nonlinear level allows some nonlinear couplings within the solver. 
    !-----------------------------------------------------------------------------------------
    IF( SteadyState ) THEN
      UpdateComponents => ListGetIntegerArray( SolverParams, &
          'Update Components', DoIt )
    ELSE
      UpdateComponents => ListGetIntegerArray( SolverParams, &
          'Nonlinear Update Components', DoIt )
    END IF
    IF( DoIt ) THEN
      CALL Info(Caller,'Updating components',Level=20)
      CALL UpdateDependentComponents( UpdateComponents )	
    END IF

    ! Compute derivative of solution with time i.e. velocity 
    ! For 2nd order schemes there is direct pointer to the velocity component
    ! Thus only 1st order schemes need to be computed.
    !-----------------------------------------------------------------------------------------
    DoIt = .FALSE.
    IF( SteadyState ) THEN
      DoIt = ListGetLogical( SolverParams,'Calculate Velocity',Found )
    ELSE
      DoIt = ListGetLogical( SolverParams,'Nonlinear Calculate Velocity',Found )
    END IF
    
    IF( DoIt ) THEN
      CALL Info(Caller,'Updating variable velocity')
      IF( .NOT. ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        CALL Warn(Caller,'Cannot calculate velocity without previous values!')
      ELSE IF( Solver % TimeOrder == 1) THEN
        dtVar => VariableGet( Solver % Mesh % Variables, 'timestep size' )
        dt = dtVar % Values(1)
        
        pVar => Solver % Variable 
        k = INDEX(pVar % Name,'[')-1
        IF( k > 0 ) THEN
          str = pVar % Name(1:k)//' Velocity'
        ELSE
          str = TRIM(pVar % Name)//' Velocity'
        END IF       
        VeloVar => VariableGet( Solver % Mesh % Variables, str )        
        VeloVar % Values = (pVar % Values - pVar % PrevValues(:,1)) / dt
      END IF
    END IF
    
    ! Finally compute potentially velocities related to exported variables.
    ! Do this on nonlinear level only when 'Nonlinear Calculate Velocity' is set true.
    !-----------------------------------------------------------------------------------------
    IF( SteadyState .OR. DoIt ) THEN          
      CALL DerivateExportedVariables( Solver )  
    END IF       

    IF( SteadyState ) THEN
      CALL Info(Caller,'Finished updating steady state objects!',Level=32)
    ELSE
      CALL Info(Caller,'Finished updating nonlinear objects!',Level=32)
    END IF

  END SUBROUTINE UpdateDependentObjects


  
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
    INTEGER :: i, n, nn, RelaxAfter, IterNo, MinIter, MaxIter, dofs
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp), POINTER :: b(:), x(:), r(:)
    REAL(KIND=dp), POINTER :: x0(:)
    REAL(KIND=dp) :: Norm, PrevNorm, rNorm, bNorm, Change, PrevChange, Relaxation, tmp(1),dt, &
        Tolerance, MaxNorm, eps, Ctarget, Poffset, nsum, dpsum
    INTEGER, TARGET  ::  Dnodes(1)
    INTEGER, POINTER :: Indexes(:)
    TYPE(Variable_t), POINTER :: iterVar, VeloVar, dtVar, WeightVar
    CHARACTER(:), ALLOCATABLE :: str,SolverName,ConvergenceType

    LOGICAL :: Stat, ConvergenceAbsolute, Relax, RelaxBefore, DoIt, Skip, &
        SkipConstraints, ResidualMode, RelativeP, NodalNorm
    TYPE(Matrix_t), POINTER :: MMatrix
    REAL(KIND=dp), POINTER CONTIG :: Mx(:), Mb(:), Mr(:)
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: TmpXVec, TmpRVec, TmpRHSVec
    INTEGER :: ipar(1)
    TYPE(ValueList_t), POINTER :: SolverParams
    CHARACTER(*), PARAMETER :: Caller = 'ComputeChange'
    LOGICAL :: Parallel, SingleMesh, x0Allocated
    
    SolverParams => Solver % Values
    RelativeP = .FALSE.
    SingleMesh = Solver % Mesh % SingleMesh

    IF(.NOT. ASSOCIATED(Solver % Variable) ) THEN
      CALL Info(Caller,'Solver variable not found for: '&
           //TRIM(ListGetString(SolverParams,'equation')),Level=10)
      RETURN
    END IF
    
    Parallel = Solver % Parallel
      
    IF(PRESENT(Matrix)) THEN
      A => Matrix
    ELSE
      A => Solver % Matrix
    END IF
            
    IF(SteadyState) THEN	
      Skip = ListGetLogical( SolverParams,'Skip Compute Steady State Change',Stat)
      IF( Skip ) THEN
        CALL Info(Caller,'Skipping the computation of steady state change',Level=15)
        RETURN
      END IF
        
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
      Relax = Relax .AND. ABS(Relaxation-1.0_dp) > EPSILON(Relaxation)

      iterVar => VariableGet( Solver % Mesh % Variables, 'coupled iter' )
      IterNo = NINT( iterVar % Values(1) )
      IF( Relax ) THEN
        RelaxAfter = ListGetInteger(SolverParams,'Steady State Relaxation After',Stat)
        IF( Stat .AND. RelaxAfter >= IterNo ) Relax = .FALSE.
      END IF	

      NodalNorm = ListGetLogical(SolverParams,'Steady State Nodal Norm',Stat)
      
      RelaxBefore = .TRUE.
      IF(Relax) THEN
        RelaxBefore = ListGetLogical( SolverParams, &
            'Steady State Relaxation Before', Stat )      
        IF (.NOT. Stat ) RelaxBefore = .TRUE.
      END IF

      ! Steady state system has never any constraints
      SkipConstraints = .FALSE.
      
    ELSE
      iterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
      IF( ASSOCIATED( iterVar ) ) THEN
        IterNo = NINT( iterVar % Values(1) )
        Solver % Variable % NonlinIter = IterNo

        Skip = ListGetLogical( SolverParams,'Skip Advance Nonlinear iter',Stat)
        IF( .NOT. Skip )  iterVar % Values(1) = IterNo + 1 

        IF( .NOT. Solver % NewtonActive ) THEN
          i = ListGetInteger( SolverParams, 'Nonlinear System Newton After Iterations',Stat )
          IF( Stat .AND. i <= IterNo ) Solver % NewtonActive = .TRUE.
        END IF
      END IF
        
      Skip = ListGetLogical( SolverParams,'Skip Compute Nonlinear Change',Stat)

      IF(Skip) THEN
        CALL Info(Caller,'Skipping the computation of nonlinear change',Level=15)
        RETURN
      END IF
        
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
      Relax = Relax .AND. ( ABS( Relaxation - 1.0_dp) > EPSILON( Relaxation ) )
      IF( Relax ) THEN
        RelaxAfter = ListGetInteger(SolverParams,'Nonlinear System Relaxation After',Stat)
        IF( Stat .AND. RelaxAfter >= Solver % Variable % NonlinIter ) Relax = .FALSE.

        RelativeP = ListGetLogical( SolverParams,'Relative Pressure Relaxation',Stat) 
        IF( RelativeP) CALL Info(Caller,'Using relative pressure relaxation',Level=10)
      END IF
      
      SkipConstraints = ListGetLogical( SolverParams,&
          'Nonlinear System Convergence Without Constraints',Stat) 

      NodalNorm = ListGetLogical(SolverParams,'Nonlinear System Nodal Norm',Stat)
      
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
    
    IF( SkipConstraints ) THEN
      n = MIN( n, A % NumberOfRows )
    END IF
      
    ! If requested (for p-elements) only use the dofs associated to nodes. 
    ! One should not optimize bandwidth if this is desired. 
    IF( NodalNorm ) THEN
      i = MAXVAL( Solver % Variable % Perm(1:Solver % Mesh % NumberOfNodes ) )
      n = MIN(n,i*Solver % Variable % Dofs)
    END IF
    
    Stat = .FALSE.
    x0 => NULL()
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
      ELSE
        x0 => Solver % Variable % Values
        Stat = .TRUE.
      END IF
    END IF
    
    x0allocated = .FALSE.
    IF(Stat .AND. .NOT. SkipConstraints ) THEN
      IF (SIZE(x0) /= SIZE(x)) THEN
         CALL Info(Caller,'WARNING: Possible mismatch in length of vectors ('&
            //I2S(SIZE(x))//' vs. '//I2S(SIZE(x0))//')!',Level=10)

         ! Try to account for differing sizes of the x & x0 (effectively "Relaxation Factor=1"
         ! for the additional "x"-vector dofs."ResidualMode" -case unhandeled.
         IF (n > SIZE(x0) .AND. .NOT. ResidualMode) THEN
           BLOCK
              REAL(KIND=dp), ALLOCATABLE :: y(:)
              INTEGER :: n0
              n0 = size(x0)
              y = x0
              ALLOCATE(x0(n))
              x0allocated = .TRUE.
              x0(1:n0) = y; x0(n0+1:n)=x(n0+1:n)
           END BLOCK
         END IF
       END IF
    END IF

    ! This ensures that the relaxation does not affect the mean of the pressure
    IF( RelativeP ) THEN
      dofs = Solver % Variable % Dofs

      dpsum = SUM(x(dofs:n:dofs)) - SUM(x0(dofs:n:dofs)) 
      nsum = 1.0_dp * n / dofs

      dpsum = ParallelReduction( dpsum ) 
      nsum = ParallelReduction( nsum )

      Poffset = (1-Relaxation) * dpsum / nsum
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
        IF( RelativeP ) x(dofs:n:dofs) = x(dofs:n:dofs) + Poffset
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
    IF( Norm /= Norm ) THEN
      PRINT *,'Norm:',Norm,PrevNorm, n
      CALL NumericalError(Caller,'Norm of solution appears to be NaN')
    END IF

    IF( SteadyState ) THEN
      MaxNorm = ListGetCReal( SolverParams, &
          'Steady State Max Norm', Stat )
    ELSE
      MaxNorm = ListGetCReal( SolverParams, &
          'Nonlinear System Max Norm', Stat )
    END IF    

    IF( Stat ) THEN
      CALL Info(Caller,Message)
      CALL NumericalError(Caller,'Norm of solution exceeded given bounds')
    END IF
      
    SELECT CASE( ConvergenceType )
        
    CASE('residual')
      !--------------------------------------------------------------------------
      ! x is solution of A(x0)x=b(x0) thus residual should really be r=b(x)-A(x)x
      ! Instead we use r=b(x0)-A(x0)x0 which unfortunately is one step behind.
      !--------------------------------------------------------------------------
      IF(PRESENT(RHS)) THEN
        b => RHS
      ELSE
        b => A % rhs
      END IF
      
      ALLOCATE(r(n))
      r=0._dp

      IF (Parallel) THEN
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
      b => A % rhs     
      
      IF (Parallel) THEN

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
      IF( .NOT. ConvergenceAbsolute ) THEN
        IF( Norm + PrevNorm > 0.0) THEN
          Change = Change * 2.0_dp/ (Norm+PrevNorm)
        END IF
      END IF
      DEALLOCATE(r)      

    CASE('norm')

      Change = ABS( Norm-PrevNorm )
      IF( .NOT. ConvergenceAbsolute .AND. Norm + PrevNorm > 0.0) THEN
        Change = Change * 2.0_dp/ (Norm+PrevNorm)
      END IF
      
    CASE DEFAULT
      CALL Warn(Caller,'Unknown convergence measure: '//TRIM(ConvergenceType))    
      
    END SELECT


    ! This could be a multislice case, for example. We don't want each slice to have
    ! different iteration count so we need to check the max norm of the partitions
    ! even for multislice case. 
    IF( SingleMesh ) THEN
      IF( ListGetInteger( CurrentModel % Simulation,'Number of Slices', Stat ) > 1 ) THEN
        CALL Info(Caller,'Communicating maximum norm in single mesh operation',Level=10)
        Change = ParallelReduction( Change, 2 )
      END IF
    END IF

    
    !--------------------------------------------------------------------------
    ! Check for convergence: 0/1
    !--------------------------------------------------------------------------
    IF(SteadyState) THEN
      PrevChange = Solver % Variable % SteadyChange
      Solver % Variable % SteadyChange = Change
      Tolerance = ListGetCReal( SolverParams,'Steady State Convergence Tolerance',Stat)
      IF( Stat ) THEN
        IF( Change <= Tolerance ) THEN
          Solver % Variable % SteadyConverged = 1
        ELSE
          Solver % Variable % SteadyConverged = 0
        END IF          
      END IF
      
      Tolerance = ListGetCReal( SolverParams,'Steady State Divergence Limit',Stat)
      IF( Stat .AND. Change > Tolerance ) THEN
        IF( IterNo > 1 .AND. Change > PrevChange ) THEN
          CALL Info(Caller,'Steady state iteration diverged over tolerance',Level=5)
          Solver % Variable % SteadyConverged = 2
        END IF
      END IF
      
      Tolerance = ListGetCReal( SolverParams,'Steady State Exit Condition',Stat)
      IF( Stat .AND. Tolerance > 0.0 ) THEN
        CALL Info(Caller,'Nonlinear iteration condition enforced by exit condition',Level=6)
        Solver % Variable % SteadyConverged = 3
      END IF

    ELSE
      PrevChange = Solver % Variable % NonlinChange 
      Solver % Variable % NonlinChange = Change
      Solver % Variable % NonlinConverged = 0

      MaxIter = ListGetInteger( SolverParams,'Nonlinear System Max Iterations',Stat)            
      
      Tolerance = ListGetCReal( SolverParams,'Nonlinear System Convergence Tolerance',Stat)
      IF( Stat ) THEN
        IF( Change <= Tolerance ) THEN
          Solver % Variable % NonlinConverged = 1
        ELSE IF( IterNo >= MaxIter ) THEN
          IF( ListGetLogical( SolverParams,'Nonlinear System Abort Not Converged',Stat ) ) THEN
            CALL Fatal(Caller,'Nonlinear iteration did not converge to tolerance')
          ELSE
            CALL Info(Caller,'Nonlinear iteration did not converge to tolerance',Level=6)
            ! Solver % Variable % NonlinConverged = 2            
          END IF
        END IF
      END IF

      Tolerance = ListGetCReal( SolverParams,'Nonlinear System Divergence Limit',Stat)
      IF( Stat .AND. Change > Tolerance ) THEN        
        IF( ( IterNo > 1 .AND. Change > PrevChange ) .OR. ( IterNo >= MaxIter ) ) THEN
          IF( ListGetLogical( SolverParams,'Nonlinear System Abort Diverged',Stat ) ) THEN
            CALL Fatal(Caller,'Nonlinear iteration diverged over limit')
          ELSE
            CALL Info(Caller,'Nonlinear iteration diverged over limit',Level=6)
            Solver % Variable % NonlinConverged = 2
          END IF
        END IF
      END IF

      Tolerance = ListGetCReal( SolverParams,'Nonlinear System Exit Condition',Stat)
      IF( Stat .AND. Tolerance > 0.0 ) THEN
        CALL Info(Caller,'Nonlinear iteration condition enforced by exit condition',Level=6)
        Solver % Variable % NonlinConverged = 3
      END IF
      
      IF( Solver % Variable % NonlinConverged > 1 ) THEN
        MinIter = ListGetInteger( SolverParams,'Nonlinear System Min Iterations',Stat)
        IF( Stat .AND. IterNo < MinIter ) THEN
          CALL Info(Caller,'Enforcing continuation of iteration',Level=6)
          Solver % Variable % NonlinConverged = 0
        END IF
      END IF
      
      IF( .NOT. Solver % NewtonActive ) THEN
        Tolerance = ListGetCReal( SolverParams, 'Nonlinear System Newton After Tolerance',Stat )
        IF( Stat .AND. Change < Tolerance ) Solver % NewtonActive = .TRUE.
      END IF     
    END IF

    
    IF(Relax .AND. .NOT. RelaxBefore) THEN
      x(1:n) = (1-Relaxation)*x0(1:n) + Relaxation*x(1:n)
      IF( RelativeP ) x(dofs:n:dofs) = x(dofs:n:dofs) + Poffset
      Solver % Variable % Norm = ComputeNorm(Solver,n,x)
    END IF
    
    ! Steady state output is done in MainUtils
    SolverName = ListGetString( SolverParams, 'Equation',Stat)
    IF(.NOT. Stat) SolverName = TRIM(Solver % Variable % Name)
 
    IF(SteadyState) THEN        
      WRITE( Message, '(a,g15.8,g15.8,a)') &
         'SS (ITER='//i2s(IterNo)//') (NRM,RELC): (',Norm, Change,&
          ' ) :: '// TRIM(SolverName)
    ELSE
      WRITE( Message, '(a,g15.8,g15.8,a)') &
         'NS (ITER='//i2s(IterNo)//') (NRM,RELC): (',Norm, Change,&
          ' ) :: '// TRIM(SolverName)
    END IF
    CALL Info( Caller, Message, Level=3 )

    ! This provides a way to directly save the convergence data into an external
    ! file making it easier to follow the progress of Elmer simulation in other software.
    !------------------------------------------------------------------------------------    
    IF( ListGetLogical( CurrentModel % Simulation,'Convergence Monitor',Stat ) ) THEN
      CALL WriteConvergenceInfo()  
    END IF
    
    ! Optional a posteriori scaling for the computed fields
    ! May be useful for some floating systems where one want to impose some intergral 
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
      IF( Parallel ) THEN
        CALL Fatal(Caller,'Setting average value not implemented in parallel!')
      END IF
      Ctarget = ListGetCReal( SolverParams,'Average Solution Value',Stat)      
      str = ListGetString( SolverParams,'Average Solution Weight Variable',Stat)
      IF( Stat ) THEN
        WeightVar => VariableGet( Solver % Mesh % Variables, str )
        IF( .NOT. ASSOCIATED( WeightVar ) ) THEN
          CALL Fatal(Caller,'> Average Solution Weight < missing: '//TRIM(str))
        END IF
        IF( SIZE(x) /= SIZE(WeightVar % Values ) ) THEN
          CALL Fatal(Caller,'Field and weight size mismatch: '//TRIM(str))          
        END IF
        Ctarget = Ctarget - SUM( WeightVar % Values * x ) / SUM( WeightVar % Values )
      ELSE
        Ctarget = Ctarget - SUM(x) / SIZE(x)
      END IF
      x = x + Ctarget
    END IF


    ! Calculate derivative a.k.a. sensitivity
    DoIt = .FALSE.
    IF( SteadyState ) THEN
      DoIt = ListGetLogical( SolverParams,'Calculate Derivative',Stat )
    ELSE
      DoIt = ListGetLogical( SolverParams,'Nonlinear Calculate Derivative',Stat )
    END IF

    IF( DoIt ) THEN
      IF( SteadyState ) THEN
        iterVar => VariableGet( Solver % Mesh % Variables, 'coupled iter' )
        IterNo = NINT( iterVar % Values(1) )
      ELSE
        iterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
        IterNo = NINT( iterVar % Values(1) )
      END IF
      
      Eps = 1.0_dp
      IF( IterNo > 1 ) THEN
        dtVar => VariableGet( Solver % Mesh % Variables, 'derivative eps' )
        IF( ASSOCIATED( dtVar ) ) THEN
          eps = dtVar % Values(1)
          Stat = .TRUE.
        ELSE
          eps = ListGetCReal( SolverParams,'derivative eps',Stat)
        END IF
        IF(.NOT. Stat) THEN
          Eps = 1.0_dp
          CALL Info(Caller,'Derivative Eps not given, using one',Level=7)
        END IF
      END IF
      
      str = GetVarname(Solver % Variable) // ' Derivative'
      VeloVar => VariableGet( Solver % Mesh % Variables, str )
      IF( ASSOCIATED( VeloVar ) ) THEN
        CALL Info(Caller,'Computing variable:'//TRIM(str),Level=7)
        VeloVar % Values = (x - x0) / eps
      ELSE
        CALL Warn(Caller,'Derivative variable not present')
      END IF

    END IF
    
    IF(.NOT. SteadyState ) THEN    
      CALL UpdateDependentObjects( Solver, .FALSE. )        
    END IF

    IF(x0allocated) DEALLOCATE(x0)

  CONTAINS

    SUBROUTINE WriteConvergenceInfo()

      INTEGER :: ConvInds(5),ConvUnit
      CHARACTER(:), ALLOCATABLE :: ConvFile
      LOGICAL, SAVE :: ConvVisited = .FALSE.

      IF( ParEnv % MyPe /= 0 ) RETURN

      ConvFile = ListGetString(CurrentModel % Simulation,&
          'Convergence Monitor File',Stat)
      IF(.NOT. Stat) ConvFile = 'convergence.dat'

      IF( ConvVisited ) THEN
        OPEN(NEWUNIT=ConvUnit, FILE=ConvFile,STATUS='old',POSITION='append')
      ELSE
        OPEN(NEWUNIT=ConvUnit, File=ConvFile)
        WRITE(ConvUnit,'(A)') '! solver  ss/ns  timestep  coupled  nonlin  norm  change'
        ConvVisited = .TRUE.
      END IF

      ConvInds = 0
      ConvInds(1) = Solver % SolverId

      IF( SteadyState ) ConvInds(2) = 1 

      iterVar => VariableGet( Solver % Mesh % Variables, 'timestep' )
      ConvInds(3) = NINT( iterVar % Values(1) )

      iterVar => VariableGet( Solver % Mesh % Variables, 'coupled iter' )
      ConvInds(4) = NINT( iterVar % Values(1) )

      iterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
      ConvInds(5) = NINT( iterVar % Values(1) )

      WRITE(ConvUnit,'(5I8,2G16.8)') ConvInds,Norm,Change
      CLOSE(ConvUnit)
      
    END SUBROUTINE WriteConvergenceInfo
        
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeChange
!------------------------------------------------------------------------------
    



!------------------------------------------------------------------------------
!> Adaptive version for getting gaussian integration points.
!> Also saves some time in initializations.
!> Note: the routine uses the pointer to Solver to check whether definitions
!> need to be remade. 
!----------------------------------------------------------------------------------------------

  FUNCTION GaussPointsAdapt( Element, Solver, PReferenceElement, EdgeBasis ) RESULT(IntegStuff)

    IMPLICIT NONE
    TYPE(Element_t) :: Element
    TYPE(Solver_t), OPTIONAL, TARGET :: Solver
    LOGICAL, OPTIONAL :: PReferenceElement           !< For switching to the p-version reference element
    LOGICAL, OPTIONAL :: EdgeBasis                   !< Choosing IP points for edge basis function
    TYPE( GaussIntegrationPoints_t ) :: IntegStuff   !< Structure holding the integration points

    CHARACTER(:), ALLOCATABLE :: VarName, GaussDef
    TYPE(Solver_t), POINTER :: pSolver, prevSolver => NULL()
    TYPE(Variable_t), POINTER :: IntegVar
    INTEGER :: AdaptOrder, AdaptNp, Np, RelOrder
    REAL(KIND=dp) :: MinLim, MaxLim, MinV, MaxV, V
    LOGICAL :: UseAdapt, Found,ElementalRule
    INTEGER :: i,n,ElementalNp(8),prevVisited = -1
    LOGICAL :: Debug, InitDone, pRef, IsBC, prevIsBC 
    INTEGER :: EdgeBasisDegree
    
    SAVE prevSolver, UseAdapt, MinLim, MaxLim, IntegVar, AdaptOrder, AdaptNp, RelOrder, Np, &
        ElementalRule, ElementalNp, prevVisited, pRef, prevIsBC 

    IF( PRESENT( Solver ) ) THEN
      pSolver => Solver
    ELSE
      pSolver => CurrentModel % Solver
    END IF
       
    ! Initialize again when we go from not bulk to bulk, or vice verse
    IsBC = ASSOCIATED( Element % BoundaryInfo )
   
    InitDone = ASSOCIATED( pSolver, prevSolver ) .AND. &
        ( prevVisited == pSolver % TimesVisited ) .AND. (.NOT. (IsBC .NEQV. PrevIsBC) )
    
    IF( .NOT. InitDone ) THEN
      PrevIsBC = IsBC 

      RelOrder = ListGetInteger( pSolver % Values,'Relative Integration Order',Found )
      AdaptNp = 0
      Np = ListGetInteger( pSolver % Values,'Number of Integration Points',Found )

      ! Elemental explicit rule will dominate over all other rules
      GaussDef = ListGetString( pSolver % Values,'Element Integration Points',ElementalRule )
      IF( ElementalRule ) THEN
        CALL ElementalGaussRules( GaussDef )
      END IF
      
      VarName = ListGetString( pSolver % Values,'Adaptive Integration Variable',UseAdapt )
      IF( UseAdapt ) THEN
        CALL Info('GaussPointsAdapt','Using adaptive gaussian integration rules',Level=7)
        IntegVar => VariableGet( pSolver % Mesh % Variables, VarName )
        IF( .NOT. ASSOCIATED( IntegVar ) ) THEN
          CALL Fatal('GaussPointsAdapt','> Adaptive Integration Variable < does not exist')
        END IF
        IF( IntegVar % TYPE /= Variable_on_nodes ) THEN
          CALL Fatal('GaussPointsAdapt','Wrong type of integration variable!')
        END IF
        MinLim = ListGetCReal( pSolver % Values,'Adaptive Integration Lower Limit' )
        MaxLim = ListGetCReal( pSolver % Values,'Adaptive Integration Upper Limit' )
        AdaptNp = ListGetInteger( pSolver % Values,'Adaptive Integration Points',Found )
        IF(.NOT. Found ) THEN
          AdaptOrder = ListGetInteger( pSolver % Values,'Adaptive Integration Order',Found )        
        END IF
        IF(.NOT. Found ) AdaptOrder = 1
        !PRINT *,'Adaptive Integration Strategy:',MinV,MaxV,AdaptOrder,AdaptNp
      END IF

      pRef = .FALSE.
      
      EdgeBasisDegree = 0
      IF( PRESENT(EdgeBasis) ) THEN
        IF( EdgeBasis ) THEN
          EdgeBasisDegree = 1          
          IF( ListGetLogical( pSolver % Values,'Quadratic Approximation',Found ) ) THEN
            EdgeBasisDegree = 2
            pRef = .TRUE.
          END IF
          IF( ListGetLogical(pSolver % Values, 'Use Piola Transform', Found) ) THEN
            pRef = .TRUE.
          END IF

          ! If elemental rule has not been given then use special edge basis rules
          ! to overrule any other rule for the gauss points. 
          IF(.NOT. ElementalRule ) THEN
            ! We can alter between the two explicit rules of edge basis using relative integration order. 
            IF( RelOrder /= 0 ) THEN
              IF( RelOrder == 1 .AND. EdgeBasisDegree == 1 ) THEN
                EdgeBasisDegree = 2
              ELSE IF( RelOrder == -1 .AND. EdgeBasisDegree == 2 ) THEN
                EdgeBasisDegree = 1
              ELSE
                CALL Warn('GaussPointsAdapt','Relative integration order does not have any effect for Edge Basis')
              END IF
            END IF
            ! This is inherited info from the rules in: EdgeElementGaussPoints
            IF( EdgeBasisDegree == 2 ) THEN
              ElementalNp = [1,2,6,9,11,27,18,27]
            ELSE IF( Pref ) THEN
              ElementalNp = [1,2,3,9,11,27,18,27]
            ELSE
              ElementalNp = [1,2,3,4,4,27,6,8]
            END IF
            ElementalRule = .TRUE.
          END IF

          IF( UseAdapt ) THEN
            CALL Fatal('GaussPointsAdapt','Adaptive rules not yet compatible with EdgeBasis')
          END IF
        END IF
      ELSE
        ! Apart from edge elements we may have p-elements defined.
        ! If not specified check from the current solver. 
        IF( PRESENT(pReferenceElement) ) THEN
          pRef = pReferenceElement 
        ELSE
          pRef = isActivePElement(Element,pSolver)
        END IF
      END IF
      
      prevSolver => pSolver
      prevVisited = pSolver % TimesVisited
    END IF
    
    IF( ElementalRule ) THEN
      ! Elemental explicit rule always has the prevalance
      Np = ElementalNp( Element % TYPE % ElementCode / 100 )
    ELSE IF( UseAdapt ) THEN
      RelOrder = 0
      Np = 0

      n = Element % TYPE % NumberOfNodes        
      MinV = MINVAL( IntegVar % Values( IntegVar % Perm( Element % NodeIndexes(1:n) ) ) )
      MaxV = MAXVAL( IntegVar % Values( IntegVar % Perm( Element % NodeIndexes(1:n) ) ) )
      
      IF( .NOT. ( MaxV < MinLim .OR. MinV > MaxLim ) ) THEN
        RelOrder = AdaptOrder
        Np = AdaptNp
      END IF
    END IF

    !IF( Debug ) PRINT *,'Adapt',UseAdapt,Element % ElementIndex, n,MaxV,MinV,MaxLim,MinLim,Np,RelOrder
    
    IF( Np > 0 ) THEN
      IntegStuff = GaussPoints( Element, Np = Np, PReferenceElement = pRef )
    ELSE IF( RelOrder /= 0 ) THEN
      IntegStuff = GaussPoints( Element, RelOrder = RelOrder, PReferenceElement = pRef )
    ELSE      
      IntegStuff = GaussPoints( Element, PReferenceElement = pRef )
    END IF

    !IF( Debug ) PRINT *,'Adapt real nodes',IntegStuff % n

  CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE ElementalGaussRules(GaussDef)
!------------------------------------------------------------------------------
      CHARACTER(LEN=*) :: GaussDef
!------------------------------------------------------------------------------
      INTEGER  :: i,j,k,n,m,iostat
      

      n = LEN_TRIM(GaussDef)
      ElementalNp = 0

      !PRINT *,'gauss def:',GaussDef(1:n)

      DO i=2,8        
        j = 0
        
        SELECT CASE( i )
        CASE( 2 )
          j =  INDEX( GaussDef(1:n), '-line' ) ! position of string "-line"
          m = 5 ! length of string "-line" after which the integer should follow
        CASE( 3 )
          j =  INDEX( GaussDef(1:n), '-tri' ) 
          m = 4
        CASE( 4 )
          j =  INDEX( GaussDef(1:n), '-quad' )
          m = 5
        CASE( 5 )
          j =  INDEX( GaussDef(1:n), '-tetra' )
          m = 6
        CASE( 6 )
          j =  INDEX( GaussDef(1:n), '-pyramid' )
          m = 8
        CASE( 7 )
          j =  INDEX( GaussDef(1:n), '-prism' )
          m = 6
        CASE( 8 )
          j =  INDEX( GaussDef(1:n), '-brick' )
          m = 6
        END SELECT
        
        IF( j > 0 ) THEN
          READ( GaussDef(j+m:n), *, IOSTAT = iostat ) k
          IF( iostat /= 0 ) THEN
            CALL Fatal('ElementGaussRules','Problems reading integer from: '//TRIM(GaussDef(j+m:n)))
          END IF
          ElementalNp(i) = k
        END IF
      END DO

      !PRINT *,'Elemental Gauss Rules:',ElementalNp
      
!------------------------------------------------------------------------------
    END SUBROUTINE ElementalGaussRules
!------------------------------------------------------------------------------
        
  END FUNCTION GaussPointsAdapt
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
    TYPE(ValueList_t), POINTER :: SolverParams
    CHARACTER(:), ALLOCATABLE :: SolverName, ConvergenceType, FileName


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
        CALL Fatal('CheckStepSize','Unknown CostMode: '//I2S(SearchMode))
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
            'xiter '//I2S(iter),1,x2,Solver % Variable % Perm ) 
        PRINT *,'Xiter range:',MINVAL(x2),MAXVAL(x2)
        NULLIFY(x2)
        
!        NULLIFY( x2 ) 
!        ALLOCATE( x2(n/2) ) 
!        x2 = x(2::2)
!        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
!            'yiter '//I2S(iter),1,x2,Solver % Variable % Perm ) 
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
      CALL Fatal('CheckStepSize','Unknown SearchMode: '//I2S(SearchMode))
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
          Change = Change * 2.0_dp / ( Norm + PrevNorm )
        END IF
      END IF

      Solver % Variable % NonlinChange = Change
      Solver % Variable % Norm = Norm
     
      IF( Solver % Variable % NonlinChange <  NonlinTol ) THEN
        Solver % Variable % NonlinConverged = 1
      END IF
      
      SolverName = ListGetString( SolverParams, 'Equation',Stat)
      IF(.NOT. Stat) SolverName = TRIM(Solver % Variable % Name)
            
      IterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter')
      m = NINT(IterVar % Values(1))
      
      ! This replaces the standard error output usually written by the ComputeChange
      WRITE( Message, '(a,g15.8,g15.8,a)') &
          'NS (ITER='//i2s(m)//') (NRM,RELC): (',Norm, Change, &
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
         
          CALL Warn('BisectSearch','Bisection method improved but faced local maximum')
          ReduceStep = .FALSE.
        ELSE 
          IF( MINVAL ( r ) < Residual0 ) THEN
            CALL Warn('BisectSearch','Bisection method improved but faced local maximum')
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
!> Apply Anderson acceleration to the solution of nonlinear system.
!> Also may apply acceleration to the linear system. 
!------------------------------------------------------------------------------
  SUBROUTINE NonlinearAcceleration(A,x,b,Solver,PreSolve,NoSolve)    
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) CONTIG :: b(:),x(:)
    TYPE(Solver_t) :: Solver
    LOGICAL :: PreSolve
    LOGICAL, OPTIONAL :: NoSolve
    !------------------------------------------------------------------------------
    ! We have a special structure for the iterates and residuals so that we can
    ! cycle over the pointers instead of the values. 
    TYPE AndersonVect_t
      LOGICAL :: Additive
      REAL(KIND=dp), POINTER :: Iterate(:), Residual(:), Ax(:)
      INTEGER :: tag
    END TYPE AndersonVect_t
    TYPE(AndersonVect_t), ALLOCATABLE :: AndersonBasis(:), AndersonTmp
    INTEGER :: AndersonInterval, ItersCnt, AndersonVecs, VecsCnt, iter, n,i,j,k
    TYPE(Variable_t), POINTER :: iterV, Svar
    REAL(KIND=dp), ALLOCATABLE :: Alphas(:),AxTable(:,:),TmpVec(:) 
    REAL(KIND=dp) :: Nrm, AndersonRelax
    LOGICAL :: Found, DoRelax, KeepBasis, Visited = .FALSE., Parallel    
    INTEGER :: LinInterval
    INTEGER :: PrevSolverId = -1
    
    SAVE AndersonBasis, TmpVec, Alphas, ItersCnt, AndersonInterval, VecsCnt, AndersonVecs, &
        PrevSolverId, AxTable, AndersonRelax, DoRelax, Visited, KeepBasis, LinInterval
        
    IF( PreSolve ) THEN
      CALL Info('NonlinearAcceleration','Performing pre-solution steps',Level=8)
    ELSE
      CALL Info('NonlinearAcceleration','Performing post-solution steps',Level=8)
    END IF

    Parallel = ( ParEnv % PEs > 1 ) 
        
    iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
    iter = NINT(iterV % Values(1))

    IF(PRESENT(NoSolve)) NoSolve = .FALSE.
    
    n = A % NumberOfRows
          
    IF(.NOT. Visited ) THEN
      PrevSolverId = Solver % SolverId
      CALL Info('NonlinearAcceleration','Allocating structures for solution history',Level=8)

      AndersonInterval = ListGetInteger( Solver % Values,&
          'Nonlinear System Acceleration Interval',Found)      
      LinInterval = ListGetInteger( Solver % Values,&
          'Linear System Acceleration Interval',Found)      

      AndersonVecs = MAX( AndersonInterval, LinInterval )
      IF( AndersonVecs == 0 ) THEN
        CALL Fatal('NonlinearAcceleration','Both acceleration intervals are zero!')
      END IF
            
      AndersonRelax = ListGetCReal( Solver % Values,&
          'Nonlinear System Acceleration Relaxation',DoRelax)
      KeepBasis = ListGetLogical( Solver % Values,&
          'Nonlinear System Acceleration Keep Vectors',Found)            

      ItersCnt = 0    ! relates to "AndersonInterval"
      VecsCnt = 0     ! relates to "AndersonVecs"
      
      IF(.NOT. ALLOCATED( AndersonBasis ) ) THEN
        ALLOCATE( AndersonBasis(AndersonVecs) )
        DO i=1,AndersonVecs
          ALLOCATE( AndersonBasis(i) % Residual(n), &
              AndersonBasis(i) % Iterate(n) )
          AndersonBasis(i) % Residual = 0.0_dp
          AndersonBasis(i) % Iterate = 0.0_dp          
        END DO
        ALLOCATE( TmpVec(n), Alphas(AndersonVecs) )
      END IF
      Visited = .TRUE.
    END IF
    
    IF( PrevSolverId /= Solver % SolverId ) THEN
      CALL Fatal('NonlinearAcceleration','Current implementation only supports one solver!')
    END IF
      
    
    IF( PreSolve ) THEN           
      IF( iter == 1 ) THEN
        ItersCnt = 0
        IF( .NOT. KeepBasis ) VecsCnt = 0
      END IF

      ItersCnt = ItersCnt + 1
      VecsCnt = VecsCnt + 1

      ! Calculate the residual of the matrix equation
      ! Here 'x' comes before being modified hence A(x) is consistent. 
      CALL MatrixVectorMultiply( A, x, TmpVec )
      TmpVec = TmpVec - b

      ! Add the iterate and residual to the basis vectors.
      ! This is fast as we operate with pointers mainly.
      AndersonTmp = AndersonBasis(AndersonVecs)        
      DO i=AndersonVecs,2,-1
        AndersonBasis(i) = AndersonBasis(i-1)
      END DO
      AndersonBasis(1) = AndersonTmp
      AndersonBasis(1) % Residual = TmpVec
      AndersonBasis(1) % Iterate = x 

      ! Pure Anderson sweep is done every AndersonInterval iterations if we have full basis.
      IF(.NOT. DoRelax .AND. AndersonInterval > 0 ) THEN
        IF( VecsCnt >= AndersonVecs .AND. ItersCnt >= AndersonInterval ) THEN
          CALL AndersonMinimize( )
          ItersCnt = 0
          IF(PRESENT(NoSolve)) NoSolve = .TRUE.
          RETURN
        END IF
      END IF
      
      IF( LinInterval > 0 ) THEN
        CALL AndersonGuess()
      END IF
    ELSE
      ! Relaxation strategy is done after each linear solve.
      IF( DoRelax ) THEN
        CALL Info('NonlinearAcceleration','Minimizing residual using history data',Level=6)
        CALL AndersonMinimize( )
      END IF
    END IF

  CONTAINS 


    !------------------------------------------------------------------------------
    FUNCTION Mydot( n, x, y ) RESULT(s)
      !------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp)  :: s
      REAL(KIND=dp) CONTIG :: x(:)
      REAL(KIND=dp) CONTIG, OPTIONAL :: y(:)
      !------------------------------------------------------------------------------
      IF ( .NOT. Parallel ) THEN
        IF( PRESENT( y ) ) THEN
          s = DOT_PRODUCT( x(1:n), y(1:n) )
        ELSE
          s = DOT_PRODUCT( x(1:n), x(1:n) )
        END IF
      ELSE
        IF( PRESENT( y ) ) THEN
          s = ParallelDot( n, x, y )
        ELSE
          s = ParallelDot( n, x, x )
        END IF
      END IF
      !------------------------------------------------------------------------------
    END FUNCTION Mydot
    !------------------------------------------------------------------------------


    !------------------------------------------------------------------------------
    SUBROUTINE Mymv( A, x, b, Update )
      !------------------------------------------------------------------------------
      REAL(KIND=dp) CONTIG :: x(:), b(:)
      TYPE(Matrix_t), POINTER :: A
      LOGICAL, OPTIONAL :: Update
      !------------------------------------------------------------------------------
      IF ( .NOT. Parallel ) THEN
        CALL CRS_MatrixVectorMultiply( A, x, b )
      ELSE
        IF ( PRESENT( Update ) ) THEN
          CALL ParallelMatrixVector( A,x,b,Update,ZeroNotOwned=.TRUE. )
        ELSE
          CALL ParallelMatrixVector( A,x,b,ZeroNotOwned=.TRUE. )
        END IF
      END IF
      !------------------------------------------------------------------------------
    END SUBROUTINE Mymv
    !------------------------------------------------------------------------------

    
    ! Given set of basis vectors and residuals find a new suggestion for solution.
    ! Either use as such or combine it to solution when relaxation is used.
    ! This is applied to boost nonlinear iteration. 
    !------------------------------------------------------------------------------
    SUBROUTINE AndersonMinimize()
      INTEGER ::m, n, AndersonMinn
      REAL(KIND=dp) :: rr, rb
      
      m = MIN( ItersCnt, AndersonInterval )      
      
      AndersonMinN = ListGetInteger( Solver % Values,&
          'Nonlinear System Acceleration First Iteration',Found )
      IF(.NOT. (Found .OR. DoRelax)) AndersonMinN = AndersonInterval
            
      ! Nothing to do 
      IF( m < AndersonMinN ) RETURN
      
      ! If size of our basis is just one, there is not much to do...
      ! We can only perform classical relaxation. 
      IF( m == 1 ) THEN
        x = AndersonRelax * x + (1-AndersonRelax) * AndersonBasis(1) % Iterate
        RETURN
      END IF
      
      ! If we are converged then the solution should already be the 1st component.
      ! Hence use that as the basis. 
      Alphas(1) = 1.0_dp     
      TmpVec = AndersonBasis(1) % Residual
      
      ! Minimize the residual
      n = SIZE( AndersonBasis(1) % Residual ) 
      DO k=2,m
        rr = MyDot( n, AndersonBasis(k) % Residual ) 
        rb = MyDot( n, AndersonBasis(k) % Residual, TmpVec )         
        Alphas(k) = -rb / rr 
        TmpVec = TmpVec + Alphas(k) * AndersonBasis(k) % Residual
      END DO

      ! Normalize the coefficients such that the sum equals unity
      ! This way for example, Dirichlet BCs will be honored.
      Alphas = Alphas / SUM( Alphas(1:m) )

      IF( InfoActive(10) ) THEN
        DO i=1,m
          WRITE(Message,'(A,I0,A,ES12.3)') 'Alpha(',i,') = ',Alphas(i)
          CALL Info('NonlinearAcceleration',Message)
        END DO
      END IF
              
      ! Create the new suggestion for the solution vector
      ! We take part of the suggested new solution vector 'x' and
      ! part of minimized residual that was used in anderson acceleration.
      IF( DoRelax ) THEN
        Alphas = Alphas * (1-AndersonRelax)
        x = AndersonRelax * x
        DO k=1,m
          x = x + Alphas(k) * AndersonBasis(k) % Iterate
        END DO
      ELSE
        x = Alphas(1) * AndersonBasis(1) % Iterate
        DO k=2,m
          x = x + Alphas(k) * AndersonBasis(k) % Iterate
        END DO
      END IF
        
    END SUBROUTINE AndersonMinimize

    
    ! Given set of basis vectors and a linear system
    ! find a combincation of the vectors that minimizes the norm of the linear
    ! system. This may be used to provide a better initial guess for a linear system.
    !--------------------------------------------------------------------------------
    SUBROUTINE AndersonGuess()
      INTEGER :: AndersonMinN

      REAL(KIND=dp), POINTER, SAVE ::Betas(:), Ymat(:,:)
      LOGICAL, SAVE :: AllocationsDone = .FALSE.
      INTEGER :: i,j,m
      
      IF(.NOT. AllocationsDone ) THEN
        m = LinInterval
        DO i=1,LinInterval
          ALLOCATE( AndersonBasis(i) % Ax(n) )
          AndersonBasis(i) % Ax = 0.0_dp
        END DO
        ALLOCATE(Betas(m),Ymat(m,m))
        AllocationsDone = .TRUE.
      END IF
      
      m = MIN( VecsCnt, LinInterval )      

      ! Calculate the residual of the matrix equation
      DO i=1,m
        CALL Mymv( A, AndersonBasis(i) % Iterate, AndersonBasis(i) % Ax )
      END DO

      DO i=1,m
        DO j=i,m
          Ymat(i,j) = SUM( AxTable(:,i) * AxTable(:,j) )
          Ymat(j,i) = Ymat(i,j)
        END DO
        Betas(i) = SUM( AxTable(:,i) * b )
      END DO
      
      CALL LUSolve(m, YMat(1:m,1:m), Betas(1:m) )

      IF( InfoActive(10) ) THEN
        DO i=1,m
          WRITE(Message,'(A,I0,A,ES12.3)') 'Beta(',i,') = ',Betas(i)
          CALL Info('NonLinearAcceleration',Message)
        END DO
      END IF
                                
      x = Betas(m) * AndersonBasis(m) % Iterate
      DO i=1,m-1
        x = x + Betas(i) * AndersonBasis(i) % Iterate
      END DO

    END SUBROUTINE AndersonGuess
    
  END SUBROUTINE NonlinearAcceleration
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
    CHARACTER(:), ALLOCATABLE :: IntVarName, MaskName
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
      IntVarName = TRIM(VarName)
    ELSE IF ( WeightAtBoundary ) THEN
      IntVarName = GetVarName(Solver % Variable) // ' Boundary Weights'
    ELSE
      IntVarName = GetVarName(Solver % Variable) // ' Weights'
    END IF
    WeightsVar => VariableGet( Mesh % Variables, IntVarName )

    IF( PRESENT( VarName ) ) THEN
      MaskName = TRIM(VarName)
    ELSE
      MaskName = "Calculate " // TRIM(IntVarName)
    END IF
      
    IF( WeightAtBoundary ) THEN
      ElemStart = Mesh % NumberOfBulkElements + 1
      ElemFin = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      UseMask = ListCheckPresentAnyBC( CurrentModel, MaskName )
    ELSE
      ElemStart = 1
      ElemFin = Mesh % NumberOfBulkElements 
      UseMask = ListCheckPresentAnyBodyForce( CurrentModel, MaskName )
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
        Solution = 0.0_dp
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

    CALL Info('CalculateNodalWeights',&
        'Computing weights for solver to variable: '//TRIM(IntVarName),Level=6)
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
          IF( .NOT. ListCheckPresent( ElemParams, MaskName ) ) CYCLE
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

    IF( ParEnv % PEs > 1 ) THEN
      IF( ASSOCIATED( Solver % Matrix ) ) THEN
        CALL ParallelSumVector(Solver % Matrix, Weights )
      ELSE
        CALL Warn('CalculateNodalWeights','Cannot communicate weights without parallel matrix!')
      END IF
    END IF
      
    DEALLOCATE(Basis, ElementNodes % x, ElementNodes % y, &
        ElementNodes % z, LocalIndexes )

    CALL Info('CalculateNodalWeights','All done',Level=8)

  END SUBROUTINE CalculateNodalWeights
!------------------------------------------------------------------------------


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
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3),Coeff
    LOGICAL :: Found, Stat, BodyElem


    CoordSys = CurrentCoordinateSystem()

    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Warn('CalculateEntityWeights','Mesh not associated!')
      RETURN
    END IF

    CALL Info('CalculateEntityWeights','Computing weights for the mesh entities',Level=6)
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
        DO bc_id = 1,NoBC
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

      ElementNodes % x = 0.0_dp
      ElementNodes % y = 0.0_dp
      ElementNodes % z = 0.0_dp
      
      ElementNodes % x(1:n) = Mesh % Nodes % x(Indexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(Indexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(Indexes)

      IntegStuff = GaussPoints( Element, Element % TYPE % GaussPoints, &
          EdgeBasis = .FALSE., PReferenceElement = .FALSE. )
      
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
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      END IF
      IF( NoBF > 0 ) THEN
        tmp_weights(1:NoBF ) = bf_weights
        CALL MPI_ALLREDUCE( tmp_weights, bf_weights, NoBF, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      END IF
      IF( NoBodies > 0 ) THEN
        tmp_weights(1:NoBodies ) = body_weights
        CALL MPI_ALLREDUCE( tmp_weights, body_weights, NoBodies, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      END IF
      IF( NoMat > 0 ) THEN
        tmp_weights(1:NoMat ) = mat_weights
        CALL MPI_ALLREDUCE( tmp_weights, mat_weights, NoMat, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
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
!>  Scale linear system with different strategies.
!------------------------------------------------------------------------------
  SUBROUTINE ScaleLinearSystem(Solver,A,b,x,DiagScaling, & 
          ApplyScaling,RhsScaling,ConstraintScaling,ScalingStr)
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), OPTIONAL :: b(:), x(:)
    REAL(KIND=dp), OPTIONAL, TARGET :: DiagScaling(:)
    LOGICAL, OPTIONAL :: ApplyScaling, RhsScaling,ConstraintScaling
    CHARACTER(*), OPTIONAL :: ScalingStr
    INTEGER :: n

    CHARACTER(:), ALLOCATABLE :: str
    LOGICAL :: Parallel, Found, ComplexMatrix
    
    n = A % NumberOfRows
    Parallel = Solver % Parallel    
    
    IF( ListGetLogical( Solver % Values,'Linear System Pseudo Complex',Found ) ) THEN
      ComplexMatrix = .TRUE.
    ELSE
      ComplexMatrix = A % COMPLEX
    END IF

    IF( PRESENT( ScalingStr ) ) THEN
      str = ScalingStr
    ELSE IF( ListGetLogical(Solver % Values,'Linear System Row Equilibration',Found ) ) THEN
      str = 'row equilibration'
    ELSE
      str = ListGetString( Solver % Values,'Linear System Scaling Method',Found )
      IF(.NOT. Found) str = 'diagonal'
    END IF
    
    SELECT CASE( str ) 
    CASE('diagonal')
      CALL ScaleLinearSystemDiagonal()

    CASE('row equilibration','rowsum')
      CALL RowEquilibration(Solver, A, b, Parallel, ApplyScaling )
      
    CASE('constant')
      CALL ScaleLinearSystemConstant()

    CASE('none')
      CALL Info('ScaleLinearSystem','No scaling will be applied!',Level=12)
      RETURN
      
    CASE DEFAULT
      CALL Fatal('ScaleLinearSystem','Unknown scaling method: '//TRIM(str))
    END SELECT
     
          
  CONTAINS

    !-------------------------------------------------------------
    !> Scale system Ax = b as (DAD)y = Db, where
    !> D = 1/SQRT(Diag(A)) and y = D^-1 x. By default, if the scaled
    !> system is rewritten as A'y = b', this subroutine also performs
    !> an additional scaling so that the final form is given by
    !> A'(y/|b'|) = b'/|b'|. Whether the last step is taken depends
    !> on the optional argument RhsScaling.
    !-------------------------------------------------------------    
    SUBROUTINE ScaleLinearSystemDiagonal()

      INTEGER :: i,j
      REAL(KIND=dp) :: bnorm,s
      COMPLEX(KIND=dp) :: DiagC
      LOGICAL :: DoRHS, DoCM
      REAL(KIND=dp), POINTER  :: Diag(:)      
      TYPE(Matrix_t), POINTER :: CM

      
      A % ScalingMethod = 1                
      
      IF( PRESENT( DiagScaling ) ) THEN
        CALL Info('ScaleLinearSystem','Reusing existing > DiagScaling < vector',Level=12)
        Diag => DiagScaling 
      ELSE
        IF(.NOT. ASSOCIATED(A % DiagScaling)) THEN
        CALL Info('ScaleLinearSystem','Creating > DiagScaling < vector of size '//I2S(n),Level=10)
          ALLOCATE( A % DiagScaling(n) ) 
        ELSE
          CALL Info('ScaleLinearSystem','Recomputing > DiagScaling < vector of size '//I2S(n),Level=12)
        END IF
        Diag => A % DiagScaling
        Diag(1:n) = 0._dp

        IF ( ComplexMatrix ) THEN
          CALL Info('ScaleLinearSystem','Assuming complex matrix while scaling',Level=20)

          !$OMP PARALLEL DO &
          !$OMP SHARED(Diag, A, N) &
          !$OMP PRIVATE(i, j) &
          !$OMP DEFAULT(NONE)
          DO i=1,n,2
            j = A % Diag(i)
            IF(j>0) THEN
              Diag(i)   = A % Values(j)
              Diag(i+1) = A % Values(j+1)
            ELSE
              Diag(i) = 0._dp
              Diag(i+1) = 0._dp
            END IF
          END DO
          !$OMP END PARALLEL DO
        ELSE
          CALL Info('ScaleLinearSystem','Assuming real valued matrix while scaling',Level=25)

          !$OMP PARALLEL DO &
          !$OMP SHARED(Diag, A, N) &
          !$OMP PRIVATE(i, j) &
          !$OMP DEFAULT(NONE)
          DO i=1,n
            j = A % Diag(i)
            IF (j>0) Diag(i) = A % Values(j)
          END DO
          !$OMP END PARALLEL DO
        END IF

        IF ( Parallel ) THEN
          CALL Info('ScaleLinearSystem','Performing parallel summation of > DiagScaling < vector',Level=20)
          CALL ParallelSumVector(A, Diag)
        END IF

        IF ( ComplexMatrix ) THEN
          !$OMP PARALLEL DO &
          !$OMP SHARED(Diag, A, N) &
          !$OMP PRIVATE(i, j, DiagC, s) &
          !$OMP DEFAULT(NONE)
          DO i=1,n,2
            DiagC = CMPLX(Diag(i),-Diag(i+1),KIND=dp)

            s = SQRT( ABS( DiagC ) )
            IF( s > TINY(s) ) THEN 
              Diag(i)   = 1.0_dp / s
              Diag(i+1) = 1.0_dp / s
            ELSE
              Diag(i)   = 1.0_dp
              Diag(i+1) = 1.0_dp
            END IF
          END DO
          !$OMP END PARALLEL DO
        ELSE
          s = 0.0_dp
          ! TODO: Add threading
          IF (ANY(ABS(Diag) <= TINY(bnorm))) s=1
          IF(Parallel) s = ParallelReduction(s,2) 

          IF(s > TINY(s) ) THEN 
            DO i=1,n
              IF ( ABS(Diag(i)) <= TINY(bnorm) ) THEN
                Diag(i) = SUM( ABS(A % Values(A % Rows(i):A % Rows(i+1)-1)) )
              ELSE
                j = A % Diag(i)
                IF (j>0) Diag(i) = A % Values(j)
              END IF
            END DO
            IF ( Parallel ) CALL ParallelSumVector(A, Diag)
          END IF

          !$OMP PARALLEL DO &
          !$OMP SHARED(Diag, N, bnorm) &
          !$OMP PRIVATE(i) &
          !$OMP DEFAULT(NONE)
          DO i=1,n
            IF ( ABS(Diag(i)) > TINY(bnorm) ) THEN
              Diag(i) = 1.0_dp / SQRT(ABS(Diag(i)))
            ELSE
              Diag(i) = 1.0_dp
            END IF
          END DO
          !$OMP END PARALLEL DO
        END IF
      END IF


      ! Optionally we may just create the diag and leave the scaling undone
      !--------------------------------------------------------------------
      IF( PRESENT( ApplyScaling ) ) THEN
        IF(.NOT. ApplyScaling ) THEN
          CALL Info('ScaleLinearSystem','Application of scaling skipped!',Level=20)
          RETURN
        END IF
      END IF

      CALL Info('ScaleLinearSystem','Scaling matrix values',Level=20)

      !$OMP PARALLEL &
      !$OMP SHARED(Diag, A, N) &
      !$OMP PRIVATE(i,j) &
      !$OMP DEFAULT(NONE)

      !$OMP DO
      DO i=1,n
        DO j = A % Rows(i), A % Rows(i+1)-1
          A % Values(j) = A % Values(j) * &
              ( Diag(i) * Diag(A % Cols(j)) )
        END DO
      END DO
      !$OMP END DO NOWAIT

      IF ( ASSOCIATED( A % PrecValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN 
          CALL Info('ScaleLinearSystem','Scaling PrecValues',Level=20)
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % PrecValues(j) = A % PrecValues(j) * &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      IF ( ASSOCIATED( A % MassValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
          CALL Info('ScaleLinearSystem','Scaling MassValues',Level=20)
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % MassValues(j) = A % MassValues(j) * &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      IF ( ASSOCIATED( A % DampValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
          CALL Info('ScaleLinearSystem','Scaling DampValues',Level=20)
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % DampValues(j) = A % DampValues(j) * &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      !$OMP END PARALLEL

      DoCM = .FALSE.
      IF(PRESENT(ConstraintScaling)) DoCm=ConstraintScaling

      IF(doCM) THEN
        CM => A % ConstraintMatrix
        IF (ASSOCIATED(CM)) THEN
          CALL Info('ScaleLinearSystem','Scaling Constraints',Level=20)
          !$OMP PARALLEL DO &
          !$OMP SHARED(Diag, CM) &
          !$OMP PRIVATE(i,j) &
          !$OMP DEFAULT(NONE)
          DO i=1,CM % NumberOFRows
            DO j=CM % Rows(i), CM % Rows(i+1)-1
              CM % Values(j) = CM % Values(j) * Diag(CM % Cols(j))
            END DO
          END DO
          !$OMP END PARALLEL DO
        END IF
      END IF

      ! Scale r.h.s. and initial guess
      !--------------------------------
      A % RhsScaling=1._dp
      ! TODO: Add threading
      IF( PRESENT( b ) ) THEN
        CALL Info('ScaleLinearSystem','Scaling Rhs vector',Level=20)

        b(1:n) = b(1:n) * Diag(1:n)
        DoRHS = .TRUE.
        IF (PRESENT(RhsScaling)) DoRHS = RhsScaling
        IF (DoRHS) THEN
          IF( Parallel ) THEN
            bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))
          ELSE
            bnorm = SQRT(SUM(b(1:n)**2))
          END IF

          IF( bnorm < SQRT( TINY( bnorm ) ) ) THEN
            CALL Info('ScaleLinearSystem','Rhs vector is almost zero, skipping rhs scaling!',Level=20)
          ELSE
            A % RhsScaling = bnorm
            b(1:n) = b(1:n) / bnorm
          END IF
        END IF
      END IF
      
      IF( PRESENT(x) ) THEN
        x(1:n) = x(1:n) / (Diag(1:n) * A % RhsScaling) 
      END IF
      
    END SUBROUTINE ScaleLinearSystemDiagonal


    !-------------------------------------------------------------
    !>  Scale system Ax = b as (A/Ascl)y = b, with y = Ascl x.
    !>  By default, if the scaled system is rewritten as A'y = b',
    !>  this subroutine also perform an additional scaling, so that
    !>  the final form is given by A'(y/bscl) = b'/bscl, i.e.
    !>  (A/Ascl)(Ascl*x/bscl) = (b/bscl). Whether the last step is
    !>  taken depends on the optional argument RhsScaling.
    !-------------------------------------------------------------    
    SUBROUTINE ScaleLinearSystemConstant()

      INTEGER :: i,j,nSum
      REAL(KIND=dp) :: bnorm,s
      COMPLEX(KIND=dp) :: DiagC
      LOGICAL :: DoRHS, DoCM
      REAL(KIND=dp) :: Ascl, bscl, Xscl, bsum, DiagSum
      TYPE(Matrix_t), POINTER :: CM
           
      IF( Parallel ) THEN
        CALL Info('ScaleLinearSystem','Scaling matrix entries by constant in parallel',Level=10)
      ELSE
        CALL Info('ScaleLinearSystem','Scaling matrix entries by constant in serial',Level=10)
      END IF
      
      CALL Info('ScaleLinearSystem','Computing > AveScaling < constant',Level=12)

      A % ScalingMethod = 3
           
      DiagSum = 0.0_dp
      nSum = n

      DO i=1,n
        j = A % Diag(i)
        IF (j>0) THEN
          DiagSum = DiagSum + ABS(A % Values(j))
        END IF
      END DO

      IF ( Parallel ) THEN
        DiagSum = ParallelReduction( DiagSum )
        nSum = ParallelReduction( nSum )
      END IF

      Ascl = DiagSum / nSum
      A % AveScaling = Ascl

      WRITE( Message,'(A,ES12.3)') 'Average diagonal entry: ', Ascl
      CALL Info( 'ScaleLinearSystemConstant', Message, Level=8 )        

      ! Optionally we may just create the diag and leave the scaling undone
      !--------------------------------------------------------------------
      IF( PRESENT( ApplyScaling ) ) THEN
        IF(.NOT. ApplyScaling ) RETURN
      END IF
      
      A % Values = A % Values / Ascl
      IF ( ASSOCIATED( A % PrecValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
          A % PrecValues = A % PrecValues / Ascl
        END IF
      END IF
      IF ( ASSOCIATED( A % MassValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
          A % MassValues = A % MassValues / Ascl
        END IF
      END IF        
      IF ( ASSOCIATED( A % DampValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
          A % DampValues = A % DampValues / Ascl
        END IF
      END IF

      DoCM = .FALSE.
      IF(PRESENT(ConstraintScaling)) DoCm=ConstraintScaling
      IF(doCM) THEN
        CM => A % ConstraintMatrix
        IF (ASSOCIATED(CM)) THEN
          CM % Values = CM % Values / Ascl
        END IF
      END IF

      ! Scale r.h.s. and initial guess
      !--------------------------------
      A % RhsScaling = 1._dp

      IF( PRESENT( b ) ) THEN       
        DoRHS = .TRUE.
        IF (PRESENT(RhsScaling)) DoRHS = RhsScaling
        IF (DoRHS) THEN          
          bsum = SUM( ABS( b(1:n) ) ) 
          nSum = n
          
          IF ( Parallel ) THEN
            bSum = ParallelReduction( bSum )
            nSum = ParallelReduction( nSum )
          END IF

          bscl = bsum / nSum           
          b = b / bscl
                   
          A % RhsScaling = bscl
          
          WRITE( Message,'(A,ES12.3)') 'Average rhs entry: ', bscl
          CALL Info( 'ScaleLinearSystemConstant', Message, Level=7 )        
        END IF
      END IF
      
      IF( PRESENT(x) ) THEN
        Xscl = A % RhsScaling / Ascl
        x(1:n) = x(1:n) / Xscl 
      END IF

    END SUBROUTINE ScaleLinearSystemConstant

!-----------------------------------------------------------------------------
  END SUBROUTINE ScaleLinearSystem
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!>   Equilibrate the rows of the coefficient matrix A to
!>   minimize the condition number. The associated rhs vector f is also scaled.
!------------------------------------------------------------------------------
  SUBROUTINE RowEquilibration( Solver, A, f, Parallel, ApplyScaling )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), TARGET :: A
    REAL(KIND=dp), OPTIONAL :: f(:)
    LOGICAL :: Parallel
    LOGICAL, OPTIONAL :: ApplyScaling 
!-----------------------------------------------------------------------------
    LOGICAL :: ComplexMatrix, Found
    INTEGER :: i, j, n 
    REAL(kind=dp) :: norm, tmp
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:), Diag(:)
!-------------------------------------------------------------------------

    CALL Info('RowEquilibration','Scaling system such that abs rowsum is unity',Level=15)

    n = A % NumberOfRows
    ComplexMatrix = A % COMPLEX

    IF (ComplexMatrix) CALL Info('RowEquilibration','Using complex matrix norm',Level=15)

    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    IF( .NOT. ASSOCIATED(A % DiagScaling) ) THEN
      ALLOCATE( A % DiagScaling(n) ) 
    END IF
    Diag => A % DiagScaling    
    
    Diag = 0.0d0
    norm = 0.0d0

    A % ScalingMethod = 2
    
    !---------------------------------------------
    ! Compute 1-norm of each row
    !---------------------------------------------
    IF (ComplexMatrix) THEN
      DO i=1,n,2
        tmp = 0.0d0
        DO j=Rows(i),Rows(i+1)-1,2
          tmp = tmp + ABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
        END DO
        Diag(i) = tmp
        Diag(i+1) = tmp
      END DO
      IF (Parallel) THEN
        CALL ParallelSUMVector(A,Diag)
      END IF
    ELSE
      IF (Parallel) THEN
 
        IF(ListGetLogical(Solver % Values, 'Row Equilibration Use M-V',Found)) THEN
          BLOCK
             REAL(KIND=dp), ALLOCATABLE :: x(:),y(:),z(:)
             TYPE(Matrix_t), POINTER :: ap
             ALLOCATE(x(n),y(n),z(n))
             ap => a
             x = 1; y=0; z=0
             CALL ParallelInitSolve( ap, x, y, z )
             CALL ParallelMatrixVector( ap, x, Diag, Update=.TRUE.,UseABS=.TRUE. )
          END BLOCK
        ELSE
          DO i=1,n
            tmp = 0.0_dp
            DO j=Rows(i),Rows(i+1)-1        
              tmp = tmp + ABS(Values(j))          
            END DO
            Diag(i)  = tmp       
          END DO
        END IF
        CALL ParallelSUMVector(A,Diag)
      ELSE
        DO i=1,n
          tmp = 0.0d0
          DO j=Rows(i),Rows(i+1)-1        
            tmp = tmp + ABS(Values(j))          
          END DO
          Diag(i)  = tmp       
        END DO
      END IF
    END IF

    norm = MAXVAL(Diag(1:n))
    IF( Parallel ) THEN
      norm = ParallelReduction(norm,2)
    END IF
    
    WRITE( Message, * ) 'Unscaled matrix norm: ', norm    
    CALL Info( 'RowEquilibration', Message, Level=5 )
    
    !--------------------------------------------------
    ! Now, define the scaling matrix by inversion and 
    ! perform the actual scaling of the linear system
    !--------------------------------------------------
    IF (ComplexMatrix) THEN    
      DO i=1,n,2
        IF (Diag(i) > TINY(norm) ) THEN
          Diag(i) = 1.0_dp / Diag(i)
        ELSE
          Diag(i) = 1.0_dp
        END IF
        Diag(i+1) = Diag(i)
      END DO
    ELSE
      DO i=1,n      
        IF (Diag(i) > TINY(norm)) THEN
          Diag(i) = 1.0_dp / Diag(i)
        ELSE
          Diag(i) = 1.0_dp
        END IF
      END DO
    END IF

    IF( PRESENT( ApplyScaling ) ) THEN
      IF(.NOT. ApplyScaling) RETURN
    END IF
    
    DO i=1,n    
      DO j=Rows(i),Rows(i+1)-1
        Values(j) = Values(j) * Diag(i)
      END DO
    END DO

    IF (PRESENT(f)) THEN
      DO i=1,n    
        f(i) = Diag(i) * f(i)
      END DO
    END IF
    
    IF ( ASSOCIATED( A % PrecValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % PrecValues(j) = A % PrecValues(j) * Diag(i) 
          END DO
        END DO
      END IF
    END IF
    
!------------------------------------------------------------------------------
  END SUBROUTINE RowEquilibration
!------------------------------------------------------------------------------


  
!--------------------------------------------------------------
!>  Scale the system back to original.
!--------------------------------------------------------------
  SUBROUTINE BackScaleLinearSystem( Solver,A,b,x,DiagScaling,&
      ConstraintScaling, EigenScaling ) 

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), OPTIONAL :: b(:),x(:)
    LOGICAL, OPTIONAL :: ConstraintScaling, EigenScaling
    REAL(KIND=dp), OPTIONAL, TARGET :: DiagScaling(:)

    CHARACTER(:), ALLOCATABLE :: str
    LOGICAL :: Found
    INTEGER :: n

    CALL Info('BackScaleLinearSystem','Scaling back to original scale',Level=14)

    n = A % NumberOfRows
    
    SELECT CASE( A % ScalingMethod )
      
    CASE( 1 ) 
      CALL BackScaleLinearSystemDiagonal()
      
    CASE( 2 ) 
      CALL ReverseRowEquilibration( A, b )
      
    CASE( 3 ) 
      CALL BackScaleLinearSystemConstant()

    CASE DEFAULT
      CALL Info('BackScaleLinearSystem','No scaling was applied!',Level=20)
      
    END SELECT
    
  CONTAINS

    SUBROUTINE BackScaleLinearSystemDiagonal()

      REAL(KIND=dp), POINTER :: Diag(:)
      INTEGER :: i,j
      LOGICAL :: doCM      
      TYPE(Matrix_t), POINTER :: CM
      
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

      ! TODO: Add threading
      ! 
      !      Solve x:  INV(D)x = y
      !      -------------------------------------------
      IF( PRESENT( x ) ) THEN
        x(1:n) = x(1:n) * Diag(1:n) * A % RhsScaling
      END IF
      
      IF( PRESENT( b ) ) THEN
        b(1:n) = b(1:n) / Diag(1:n) * A % RhsScaling 
      END IF

      IF( PRESENT( EigenScaling ) ) THEN
        IF( EigenScaling ) THEN
          ! TODO: Add threading
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
        END IF
      END IF

      !$OMP PARALLEL &
      !$OMP SHARED(Diag, A, N) &
      !$OMP PRIVATE(i, j) &
      !$OMP DEFAULT(NONE)

      !$OMP DO
      DO i=1,n
        DO j=A % Rows(i), A % Rows(i+1)-1
          A % Values(j) = A % Values(j) / (Diag(i) * Diag(A % Cols(j)))
        END DO
      END DO
      !$OMP END DO NOWAIT

      IF ( ASSOCIATED( A % PrecValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % PrecValues(j) = A % PrecValues(j) / &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      IF ( ASSOCIATED( A % MassValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % MassValues(j) = A % MassValues(j) / &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      IF ( ASSOCIATED( A % DampValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % DampValues(j) = A % DampValues(j) / &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      !$OMP END PARALLEL

      ! TODO: Add threading
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

      IF(.NOT. PRESENT(DiagScaling) ) THEN
        DEALLOCATE(A % DiagScaling)
        A % DiagScaling=>NULL()
      END IF
        
    END SUBROUTINE BackScaleLinearSystemDiagonal


    SUBROUTINE BackScaleLinearSystemConstant()
      
      REAL(KIND=dp) :: Ascl, Bscl, Xscl
      LOGICAL :: doCM      
      TYPE(Matrix_t), POINTER :: CM
            
      Ascl = A % AveScaling
      Bscl = A % RhsScaling
      Xscl = Bscl / Ascl

      IF( PRESENT( x ) ) THEN
        x(1:n) = x(1:n) * Xscl 
      END IF

      IF( PRESENT( b ) ) THEN
        b(1:n) = Bscl * b(1:n) 
      END IF

      IF( PRESENT( EigenScaling ) ) THEN
        IF( EigenScaling ) THEN
          Solver % Variable % EigenVectors = &
              Solver % Variable % EigenVectors * Xscl
        END IF
      END IF
                 
      IF( Ascl > 0.0 .AND. ABS(Ascl-1.0_dp) > EPSILON(Ascl) ) THEN
        A % Values = Ascl * A % Values
        IF ( ASSOCIATED( A % PrecValues ) ) THEN
          IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
            A % PrecValues = AScl * A % PrecValues
          END IF
        END IF
        IF ( ASSOCIATED( A % MassValues ) ) THEN
          IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
            A % MassValues = Ascl * A % MassValues
          END IF
        END IF
        IF ( ASSOCIATED( A % DampValues ) ) THEN
          IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
            A % DampValues = Ascl * A % DampValues
          END IF
        END IF
        
        doCM = .FALSE.
        IF(PRESENT(ConstraintScaling)) doCM=ConstraintScaling
        IF(doCM) THEN
          CM => A % ConstraintMatrix
          IF( ASSOCIATED(CM ) ) THEN
            CM % Values = Ascl * CM % Values
          END IF
        END IF
      END IF

      A % RhsScaling = 1._dp

    END SUBROUTINE BackScaleLinearSystemConstant

    
  END SUBROUTINE BackScaleLinearSystem
      

!------------------------------------------------------------------------------
!> Scale the linear system back to original when the linear
!> system scaling has been done by row equilibration.
!------------------------------------------------------------------------------
  SUBROUTINE ReverseRowEquilibration( A, f )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), OPTIONAL :: f(:)
!-----------------------------------------------------------------------------
    INTEGER :: i, j, n
    INTEGER, POINTER :: Rows(:)
    REAL(KIND=dp), POINTER :: Values(:), Diag(:)
!-----------------------------------------------------------------------------
    CALL Info('ReverseRowEquilibration','Scaling back to original scale',Level=14)

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

    IF (PRESENT(f)) f(1:n) = f(1:n) / Diag(1:n)
    DO i=1,n    
      DO j = Rows(i), Rows(i+1)-1
        Values(j) = Values(j) / Diag(i)
      END DO
    END DO

    IF ( ASSOCIATED( A % PrecValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % PrecValues(j) = A % PrecValues(j) / Diag(i) 
          END DO
        END DO
      END IF
    END IF
    
    DEALLOCATE(A % DiagScaling)
    A % DiagScaling => NULL()

!------------------------------------------------------------------------------
  END SUBROUTINE ReverseRowEquilibration
!------------------------------------------------------------------------------


  SUBROUTINE CalculateLoads( Solver, Aaid, x, DOFs, UseBulkValues, NodalLoads, NodalValues ) 

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), POINTER  :: Aaid
    REAL(KIND=dp) CONTIG :: x(:)
    INTEGER :: DOFs
    LOGICAL :: UseBulkValues
    TYPE(Variable_t), POINTER, OPTIONAL :: NodalLoads
    REAL(KIND=dp), POINTER, OPTIONAL :: NodalValues(:)
    
    INTEGER :: i,j,k,l,m,ii,This,DOF
    REAL(KIND=dp), POINTER :: TempRHS(:), TempVector(:), Rhs(:), TempX(:)
    REAL(KIND=dp), POINTER CONTIG :: SaveValues(:)
    REAL(KIND=dp) :: Energy, Energy_im
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL :: Found, Rotated
    REAL(KIND=dp), ALLOCATABLE :: BoundarySum(:), BufReal(:)
    INTEGER, ALLOCATABLE :: BoundaryShared(:),BoundaryActive(:),DofSummed(:),BufInteg(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: bc, ind, NoBoundaryActive, NoBCs, ierr
    LOGICAL :: OnlyGivenBCs, IgnorePeriodic
    LOGICAL :: UseVar, Parallel


    Parallel = Solver % Parallel
      
    UseVar = .FALSE.
    IF(PRESENT( NodalLoads ) ) THEN
      UseVar = ASSOCIATED( NodalLoads )
      IF(.NOT. UseVar ) THEN
        CALL Warn('CalculateLoads','Load variable not associated!')
        RETURN
      END IF
    ELSE IF( PRESENT( NodalValues ) ) THEN
      IF(.NOT. ASSOCIATED( NodalValues ) ) THEN
        CALL Warn('CalculateLoads','Load values not associated!')
        RETURN
      END IF
    ELSE
      CALL Fatal('CalculateLoads','Give either loads variable or values as parameter!')
    END IF
    
    ALLOCATE( TempVector(Aaid % NumberOfRows) )

    IF( UseBulkValues ) THEN
      IF(.NOT. ASSOCIATED(Aaid % BulkValues)) THEN
        CALL Fatal('CalculateLoads','"BulkValues" are not associated!')
      END IF
      SaveValues => Aaid % Values      
      Aaid % Values => Aaid % BulkValues
      Rhs => Aaid % BulkRHS
    ELSE
      Rhs => Aaid % Rhs
    END IF
    
    IF ( Parallel ) THEN
      IF( ASSOCIATED( Rhs ) ) THEN
        ALLOCATE(TempRHS(SIZE(Rhs)))
        TempRHS = Rhs 
        CALL ParallelInitSolve( Aaid, x, TempRHS, Tempvector )
      END IF
      CALL ParallelMatrixVector( Aaid, x, TempVector, .TRUE. )
    ELSE
      CALL MatrixVectorMultiply( Aaid, x, TempVector )
    END IF

    IF(ListGetLogical(Solver % Values, 'Calculate Energy Inner Product', Found) .OR. &
        ListGetLogical(Solver % Values, 'Calculate Energy Norm', Found) ) THEN
      Energy = 0._dp
      IF( ListGetLogical(Solver % Values, 'Linear System Complex', Found) ) THEN
        Energy_im = 0._dp
        DO i = 1, (Aaid % NumberOfRows / 2)
          IF ( Parallel ) THEN
            IF ( Aaid% ParMatrix % ParallelInfo % &
                NeighbourList(2*(i-1)+1) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
          END IF
          Energy    = Energy    + x(2*(i-1)+1) * TempVector(2*(i-1)+1) + x(2*(i-1)+2) * TempVector(2*(i-1)+2)
          Energy_im = Energy_im + x(2*(i-1)+1) * TempVector(2*(i-1)+2) - x(2*(i-1)+2) * TempVector(2*(i-1)+1) 
        END DO
        Energy    = ParallelReduction(Energy)
        Energy_im = ParallelReduction(Energy_im)

        CALL ListAddConstReal( Solver % Values, 'Energy inner product', Energy)
        CALL ListAddConstReal( Solver % Values, 'Energy inner product im', Energy_im)

        WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy inner product'
        CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy )

        WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy inner product im'
        CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy_im )

        WRITE( Message, * ) 'Energy inner product: ', Energy, Energy_im
        CALL Info( 'CalculateLoads', Message, Level=5)
      ELSE 
        DO i=1,Aaid % NumberOfRows
          IF ( Parallel ) THEN
            IF ( Aaid % ParMatrix % ParallelInfo % &
                NeighbourList(i) % Neighbours(1) /= Parenv % MyPE ) CYCLE
          END IF
          Energy = Energy + x(i)*TempVector(i)
        END DO
        Energy = ParallelReduction(Energy)
        CALL ListAddConstReal( Solver % Values, 'Energy inner product', Energy )

        WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy inner product'
        CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy )

        WRITE( Message, * ) 'Energy inner product: ', Energy
        CALL Info( 'CalculateLoads', Message, Level=5)
      END IF
    END IF

    IF( ASSOCIATED( Rhs ) ) THEN
      IF ( Parallel ) THEN
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
    END IF
          
    IgnorePeriodic = ListGetLogical( Solver % Values,'Calculate Loads Ignore Periodic',Found )
    
    NoBCs = CurrentModel % NumberOfBCs

    IF( IgnorePeriodic ) THEN
      DO This=1,NoBCs
        Projector => CurrentModel  % BCs(This) % PMatrix
        IF (ASSOCIATED(Projector)) THEN
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
    END IF
      
    IF( UseVar ) THEN
      DO i=1,SIZE( NodalLoads % Perm )
        j = NodalLoads % Perm(i)
        k = Solver % Variable % Perm(i)
        IF ( j>0 .AND. k>0 ) THEN
          DO dof=1,DOFs
            NodalLoads % Values(DOFs*(j-1)+dof) = TempVector(DOFs*(k-1)+dof)
          END DO
        END IF
      END DO
    ELSE
      NodalValues = TempVector
    END IF
      
    IF( ListGetLogical( Solver % Values,'Calculate Boundary Fluxes',Found ) ) THEN
      CALL Info('CalculateLoads','Computing boundary fluxes from nodal loads',Level=6)

      IF( Solver % Mesh % MaxEdgeDofs > 1 .OR. Solver % Mesh % MaxFaceDOFs > 1 ) THEN
        CALL Warn('CalculateLoads','Boundary flux computation implemented only for nodes for now!')
      END IF

      IF(.NOT. UseVar ) THEN
        CALL Fatal('CalculateLoads','Boundary flux computation needs the variable parameter!')        
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
          IF ( Parallel ) THEN
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
      IF( Parallel ) THEN
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
          CALL Warn('CalculateLoads','Boundary '//I2S(i)//' includes inseparable dofs!')
        END IF
        NoBoundaryActive = NoBoundaryActive + 1

        DO j=1,DOFs
          IF( Dofs == 1 ) THEN
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' Flux over BC '//I2S(i)
          ELSE
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' '//I2S(j)//' Flux over BC '//I2S(i)
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
                ' '//I2S(j)//' Flux over all BCs'
          END IF
          WRITE( Message,'(A,ES12.5)') TRIM(Message)//': ',SUM(BoundarySum(j::DOFs))
          CALL Info('CalculateLoads',Message,Level=6)
        END DO
      END IF
      
      DEALLOCATE( DofSummed, BoundaryShared, BoundaryActive, BoundarySum )      


      IF( ListGetLogical( Solver % Values,'Calculate Boundary Weights', Found ) ) THEN
        BLOCK
          CHARACTER(MAX_NAME_LEN) :: Name   
          TYPE(Variable_t), POINTER :: WVar, FVar
          REAL(KIND=dp) :: eps

          CALL Info('CalculateLoads','Computing fluxes on boundaries!',Level=10)
          
          Name = GetVarName(Solver % Variable) // ' Boundary Weights'
          WVar => VariableGet( Solver % Mesh % Variables, Name )

          IF(.NOT. ASSOCIATED(WVar) ) THEN
            CALL Fatal('CalculateLoads','Weight variable is not available!')
          END IF

          Name = GetVarName(Solver % Variable) // ' Boundary Flux'        
          FVar => VariableGet( Solver % Mesh % Variables, Name ) 
          IF(.NOT. ASSOCIATED(FVar) ) THEN
            CALL VariableAddVector( Solver % Mesh % Variables,&
                Solver % Mesh, Solver, Name, Solver % Variable % DOFs, Secondary = .TRUE., &
                Perm = WVar % Perm )
            FVar => VariableGet( Solver % Mesh % Variables, Name ) 
          END IF
          IF(.NOT. ASSOCIATED(FVar) ) THEN
            CALL Fatal('CalculateLoads','Flux variable is not available!')
          END IF

          eps = EPSILON(eps)
          FVar % Values = 0.0_dp
          
          DO i=1,SIZE( WVar % Perm )
            j = WVar % Perm(i)
            IF(j==0) CYCLE
            k = Solver % Variable % Perm(i)
            IF(k==0) CYCLE

            DO dof=1,DOFs
              IF( WVar % Values(j) < eps ) CYCLE 
              FVar % Values(DOFs*(j-1)+dof) = TempVector(DOFs*(k-1)+dof) / WVar % Values(j)
            END DO
          END DO

        END BLOCK
      END IF
    END IF
    
    DEALLOCATE( TempVector )
    
    IF( UseBulkValues ) THEN
      Aaid % Values => SaveValues
    END IF

  END SUBROUTINE CalculateLoads





  ! Create a boundary matrix and at calculate step compute the boundary loads
  ! for one given body. This is not called by default but the user needs to
  ! include it in the code, both at assembly and after solution.
  !-----------------------------------------------------------------------------
  SUBROUTINE BCLoadsAssembly( Solver, Element, LocalMatrix, LocalRhs )

    TYPE(Solver_t) :: Solver
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: LocalMatrix(:,:)
    REAL(KIND=dp) :: LocalRhs(:)

    LOGICAL :: FirstStep
    INTEGER :: i,j,k,l,n,Row,Col,Dofs,ElemNo,TargetBody=-1
    TYPE(Matrix_t), POINTER :: BCMat
    REAL(KIND=dp) :: Val
    LOGICAL :: Found
    INTEGER, POINTER :: Perm(:), BCPerm(:)
    CHARACTER(:), ALLOCATABLE :: Name   
    TYPE(Variable_t), POINTER :: BCVar


    SAVE :: BCMat, TargetBody, BCPerm, Perm, Dofs


    FirstStep = ( Solver % ActiveElements(1) == Element % ElementIndex )

    IF( FirstStep ) THEN
      CALL Info('BCLoadsAssembly','Visiting first element',Level=6)
 
      BCMat => Solver % Matrix % EMatrix
      IF(.NOT. ASSOCIATED( BCMat ) ) THEN
        TargetBody = ListGetInteger( Solver % Values,'Boundary Loads Target Body',Found )
        IF( Found ) THEN
          CALL Info('BCLoadsAssembly','Target body set to: '//I2S(TargetBody),Level=6)       
        ELSE
          TargetBody = -1
          RETURN
        END IF

        CALL Info('BCLoadsAssembly','Allocating structures for load computation',Level=8)
        IF ( ParEnv % PEs > 1 ) THEN
          CALL Warn('BCLoadsAssembly','Not implemented in parallel')
        END IF

        ! Mark the active nodes
        ALLOCATE( BCPerm( Solver % Mesh % NumberOfNodes ) )
        BCPerm = 0

        ElemNo = 0
        k = Solver % Mesh % NumberOfBulkElements
        DO i = k+1,k + Solver % Mesh % NumberOfBoundaryElements
          Element => Solver % Mesh % Elements(i)
          Found = .FALSE.             
          IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
            Found = ( Element % BoundaryInfo % Left % BodyId == TargetBody )
          END IF
          IF(.NOT. Found ) THEN
            IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
              Found = ( Element % BoundaryInfo % Right % BodyId == TargetBody )
            END IF
          END IF
          IF( Found ) THEN
            ElemNo = ElemNo + 1
            BCPerm( Element % NodeIndexes ) = 1
          END IF
        END DO

        CALL Info('BCLoadsAssembly','Number of related boundary elements: '//I2S(ElemNo),Level=8)

        n = 0
        DO i=1,Solver % Mesh % NumberOfNodes
          IF( BCPerm(i) > 0 ) THEN
            n = n + 1
            BCPerm(i) = n
          END IF
        END DO
        CALL Info('BCLoadsAssembly','Number of active nodes: '//I2S(n),Level=8)

        ! Create the list matrix 
        BCMat => AllocateMatrix()
        BCMat % Format = MATRIX_LIST           
        CALL AddToMatrixElement( BCMat, n, n, 0.0_dp )
        Solver % Matrix % EMatrix => BCMat

        ALLOCATE( BCMat % Rhs(n) )
        BCMat % Rhs = 0.0_dp
      END IF

      ! When visiting the routine after the 1st iteration the matrix for is already CRS 
      IF( BCMat % Format == MATRIX_CRS ) THEN
        BCMat % Values = 0.0_dp
        BCMat % Rhs = 0.0_dp
      END IF

      Dofs = Solver % Variable % Dofs
      Perm => Solver % Variable % Perm

      Name = TRIM(Solver % Variable % Name)//' BCLoads'
      BCVar => VariableGet( Solver % Mesh % Variables, TRIM( Name ) )
      IF(.NOT. ASSOCIATED( BCVar ) ) THEN
        CALL Info('CalculateBCLoads','Creating variable: '//TRIM(Name),Level=7)
        CALL VariableAddVector( Solver % Mesh % Variables,&
            Solver % Mesh, Solver, Name, DOFs, Perm = BCPerm )
      END IF
      
    END IF

    IF( Element % BodyId == TargetBody ) THEN       
      n = Element % TYPE % NumberOfNodes
      DO i=1,n
        IF ( BCPerm( Element % NodeIndexes(i) ) == 0 ) CYCLE
        DO k=0,Dofs-1
          Row = Dofs * BCPerm( Element % NodeIndexes(i) ) - k
          BCMat % rhs(Row) = BCMat % rhs(Row) + LocalRhs(Dofs*i-k)
          DO j=1,n
            DO l=0,Dofs-1
              Col = Dofs * Perm( Element % NodeIndexes(j) ) - l
              Val = LocalMatrix(Dofs*i-k,Dofs*j-l)
              CALL AddToMatrixElement(BCMat,Row,Col,Val)
            END DO
          END DO
        END DO
      END DO
    END IF


  END SUBROUTINE BCLoadsAssembly


  ! Calculate the boundary loads resulting from the action of boundary matrix.
  !-----------------------------------------------------------------------------
  SUBROUTINE BCLoadsComputation( Solver )

    TYPE(Solver_t) :: Solver

    TYPE(Matrix_t), POINTER :: BCMat
    CHARACTER(:), ALLOCATABLE :: Name   
    TYPE(Variable_t), POINTER :: BCVar


    BCMat => Solver % Matrix % EMatrix
    IF(.NOT. ASSOCIATED( BCMat ) ) THEN
      CALL Fatal('BCLoadsComputation','We should have the boundary matrix!')
    END IF
        
    CALL Info('BCLoadsComputation','Computing boundary loads',Level=6)
    IF( BCMat % FORMAT == MATRIX_LIST ) THEN
      CALL List_ToCRSMatrix( BCMat )
      CALL Info('BCLoadsComputation','Matrix format changed to CRS',Level=8)
    END IF

    Name = TRIM(Solver % Variable % Name)//' BCLoads'
    BCVar => VariableGet( Solver % Mesh % Variables, TRIM( Name ) )
    IF(.NOT. ASSOCIATED( BCVar ) ) THEN
      CALL Fatal('BCLoadsComputation','Variable not present: '//TRIM(Name))
    END IF
    
    CALL MatrixVectorMultiply( BCMat, Solver % Variable % Values, BCVar % Values )
    BCVar % Values = BCVar % Values - BCMat % rhs

    CALL Info('BCLoadsComputation','All done',Level=12)

  END SUBROUTINE BCLoadsComputation


    
!------------------------------------------------------------------------------
!> Prints the values of the CRS matrix to standard output.
!------------------------------------------------------------------------------
  SUBROUTINE PrintMatrix( A, Parallel, CNumbering,SaveMass, SaveDamp, SaveStiff, SkipZeros )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A            !< Structure holding matrix
    LOGICAL :: Parallel    !< are we in parallel mode?
    LOGICAL :: CNumbering  !< Continuous numbering ?
    LOGICAL, OPTIONAL :: SaveMass  !< Should we save the mass matrix
    LOGICAL, OPTIONAL :: SaveDamp  !< Should we save the damping matrix
    LOGICAL, OPTIONAL :: SaveStiff !< Should we save the stiffness matrix
    LOGICAL, OPTIONAL :: SkipZeros !< Should we write zeros or not
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,IndMass,IndDamp,IndStiff,IndMax,row,col
    LOGICAL :: DoMass, DoDamp, DoStiff, Found, Skip0
    REAL(KIND=dp) :: Vals(3), val
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

    Skip0 = .FALSE.
    IF(PRESENT(SkipZeros)) Skip0 = SkipZeros
        
    IF(.NOT. (DoStiff .OR. DoDamp .OR. DoMass ) ) THEN
      CALL Warn('CRS_PrintMatrix','Saving just the topology!')
    END IF
    
    IndStiff = 0
    IndDamp = 0
    IndMass = 0

    IF( DoStiff ) IndStiff = 1
    IF( DoDamp ) IndDamp = IndStiff + 1
    IF( DoMass ) IndMass = MAX( IndStiff, IndDamp ) + 1
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

        IF( DoStiff ) THEN
          Vals(IndStiff) = A % Values(j)
        END IF
        IF( DoDamp ) THEN
          Vals(IndDamp) = A % DampValues(j)
        END IF
        IF( DoMass ) THEN
          Vals(IndMass) = A % MassValues(j)
        END IF

        IF( Skip0 ) THEN
          IF(SUM(ABS(Vals(1:IndMax))) < EPSILON(val)) CYCLE
        END IF
          
        WRITE(1,'(I0,A,I0,A)',ADVANCE='NO') row,' ',col,' '
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
  SUBROUTINE PrintRHS( A, Parallel, CNumbering, SaveSum )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A  !< Structure holding matrix
    LOGICAL :: Parallel, CNumbering, SaveSum
!------------------------------------------------------------------------------
    INTEGER :: i, j, row
    REAL(KIND=dp) :: val, asum, rsum
    LOGICAL :: SaveRhs

    SaveRhs = ASSOCIATED(A % rhs)
    
    DO i=1,A % NumberOfRows
      row = i
      IF(Parallel) THEN
        IF(Cnumbering) THEN
          row = A % Gorder(i)
        ELSE 
          row = A % ParallelInfo % GlobalDOFs(i)
        END IF
      END IF

      IF(SaveRhs) val = A % rhs(i)

      IF(SaveSum) THEN
        rsum = 0.0_dp
        asum = 0.0_dp
        DO j = A % Rows(i),A % Rows(i+1)-1
          rsum = rsum + A % Values(j)
          asum = asum + ABS(A % Values(j))
        END DO
      END IF
      
      WRITE(1,'(I0,A)',ADVANCE='NO') row,' '
      IF( SaveRhs .AND. SaveSum ) THEN
        WRITE(1,*) Val, rsum, asum
      ELSE IF( .NOT. SaveRhs .AND. SaveSum ) THEN
        WRITE(1,*) rsum, asum
      ELSE IF( ABS( Val ) <= TINY( val ) ) THEN
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
    USE EigenSolve

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
               SkipZeroRhs, SkipLoads, ComplexSystem, ComputeChangeScaled, &
               RecursiveAnalysis, CalcLoads
    INTEGER :: n,i,j,k,l,ii,m,DOF,istat,this,mn, AllocStat
    CHARACTER(:), ALLOCATABLE :: Method, Prec, SaveSlot
    INTEGER(KIND=AddrInt) :: Proc
    REAL(KIND=dp), ALLOCATABLE, TARGET :: Px(:), &
                TempRHS(:), NonlinVals(:)
    REAL(KIND=dp), POINTER :: Diag(:)
    REAL(KIND=dp) :: s,Relaxation,Beta,Gamma,bnorm,Energy,xn,bn
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: Aaid, Projector, MP
    REAL(KIND=dp), POINTER :: mx(:), mb(:), mr(:)
    TYPE(Variable_t), POINTER :: IterV
    LOGICAL :: NormalizeToUnity, AndersonAcc, AndersonScaled, NoSolve, Found
    REAL(KIND=dp), POINTER :: pv(:)
    CHARACTER(*), PARAMETER :: Caller = 'SolveLinearSystem'
    CHARACTER(LEN=MAX_NAME_LEN) :: str

    
    TARGET b, x 
    
    INTERFACE 
       SUBROUTINE VankaCreate(A,Solver)
          USE Types
          TYPE(Matrix_t) :: A
          TYPE(Solver_t) :: Solver
       END SUBROUTINE VankaCreate

       SUBROUTINE CircuitPrecCreate(A,Solver)
          USE Types
          TYPE(Matrix_t), TARGET :: A
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

    Params => Solver % Values
     
    IF( ListGetLogical( Params,'Linear System Skip Complex',GotIt ) ) THEN
      CALL Info(Caller,'This time skipping complex treatment',Level=20)
      A % COMPLEX = .FALSE.
      ComplexSystem = .FALSE.
    ELSE
      ComplexSystem = ListGetLogical( Params, 'Linear System Complex', GotIt )
      IF ( GotIt ) A % COMPLEX = ComplexSystem
    END IF

    IF( ListGetLogical( Params,'Linear System Skip Scaling',GotIt ) ) THEN     
      CALL Info(Caller,'This time skipping scaling',Level=20)
      ScaleSystem = .FALSE.
    ELSE
      ScaleSystem = ListGetLogical( Params, 'Linear System Scaling', GotIt )
      IF ( .NOT. GotIt  ) ScaleSystem = .TRUE.
    END IF

    SkipLoads = ListGetLogical( Params,'Linear System Skip Loads',GotIt)
    
    
    IF( A % COMPLEX ) THEN
      CALL Info(Caller,'Assuming complex valued linear system',Level=6)
    ELSE
      CALL Info(Caller,'Assuming real valued linear system',Level=8)
    END IF

    Parallel = Solver % Parallel
    
!------------------------------------------------------------------------------
!   If parallel execution, check for parallel matrix initializations
!------------------------------------------------------------------------------
    IF ( Parallel  ) THEN
      IF( .NOT. ASSOCIATED(A % ParMatrix) ) THEN
        CALL Info(Caller,'Creating parallel matrix structures',Level=8)
        CALL ParallelInitMatrix( Solver, A )
      ELSE
        CALL Info(Caller,'Using previously created parallel matrix structures!',Level=15)
      END IF      
      Parallel = ASSOCIATED(A % ParMatrix)       
    END IF

    IF( Parallel ) THEN
      CALL Info(Caller,'Assuming parallel linear system',Level=8)
    ELSE
      CALL Info(Caller,'Assuming serial linear system',Level=8)
    END IF  
        
    IF ( ListGetLogical( Solver % Values, 'Linear System Save',GotIt )) THEN
      saveslot = ListGetString( Solver % Values,'Linear System Save Slot', GotIt )
      IF(SaveSlot == 'linear solve') CALL SaveLinearSystem( Solver, A )
    END IF

!------------------------------------------------------------------------------

    n = A % NumberOfRows

    BackRotation = ListGetLogical(Params,'Back Rotate N-T Solution',GotIt)
    IF (.NOT.GotIt) BackRotation=.TRUE.
    BackRotation = BackRotation .AND. ASSOCIATED(Solver % Variable % Perm)

    IF ( Solver % Matrix % Lumped .AND. Solver % TimeOrder == 1 ) THEN
       Method = ListGetString( Params, 'Timestepping Method', GotIt)
       IF (  Method == 'runge-kutta' .OR. Method == 'explicit euler' ) THEN
         ALLOCATE(Diag(n), TempRHS(n))

         TempRHS= b(1:n)
         Diag = A % Values(A % Diag)

         IF( Parallel ) THEN
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

    EigenAnalysis = .FALSE.
    IF(.NOT. ListGetLogical( Params,'Skip Eigen Analysis',GotIt)) THEN
      EigenAnalysis = Solver % NOFEigenValues > 0 .AND. &
          ListGetLogical( Params, 'Eigen Analysis',GotIt )
    END IF
          
    HarmonicAnalysis = ( Solver % NOFEigenValues > 0 ) .AND. &
        ListGetLogical( Params, 'Harmonic Analysis',GotIt )
    
    ! These analyses types may require recursive strategies and may also have zero rhs
    RecursiveAnalysis = HarmonicAnalysis .OR. EigenAnalysis     

    ApplyLimiter = ListGetLogical( Params,'Apply Limiter',GotIt ) 
    SkipZeroRhs = ListGetLogical( Params,'Skip Zero Rhs Test',GotIt )
#ifdef HAVE_FETI4I
    IF ( C_ASSOCIATED(A % PermonMatrix) ) THEN
      ScaleSystem = .FALSE.
      SkipZeroRhs = .TRUE.
    END IF
#endif

    IF ( .NOT. ( RecursiveAnalysis .OR. ApplyLimiter .OR. SkipZeroRhs ) ) THEN
      bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))      
      IF ( bnorm <= TINY( bnorm) ) THEN
        CALL Info(Caller,'Solution trivially zero!',Level=5)
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
      
      IF ( ScaleSystem ) CALL BackScaleLinearSystem( Solver, A, EigenScaling = .TRUE. ) 
      IF ( BackRotation ) CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )

      Norm = ComputeNorm(Solver,n,x)
      Solver % Variable % Norm = Norm
      
      NormalizeToUnity = ListGetLogical( Solver % Values, &
          'Eigen System Normalize To Unity',Stat)         

      IF(NormalizeToUnity .OR. ListGetLogical( Solver % Values,  &
          'Eigen System Mass Normalize', Stat ) ) THEN

        CALL ScaleEigenVectors( A, Solver % Variable % EigenVectors, &
            SIZE(Solver % Variable % EigenValues), NormalizeToUnity ) 
      END IF

      CALL Info(Caller, 'Repointing '//I2S(Solver % Variable % DOFs)//&
          ' eigenvalue components for: '//TRIM(Solver % Variable % Name))

      IF (A % Complex) THEN
        DO k=1,Solver % Variable % DOFs
          str = ComponentName( Solver % Variable % Name, k )
          Var => VariableGet( Solver % Mesh % Variables, str, .TRUE. )

          IF( ASSOCIATED( Var ) ) THEN
            IF (Solver % Variable % DOFs == 2) THEN
              CALL Info('SolveEigenSystem', 'Eigenvalue component ' &
                  //I2S(k)//': '//TRIM(str))
              Var % EigenValues => Solver % Variable % EigenValues
              IF (.NOT. ASSOCIATED(Var % EigenVectors)) THEN
                ALLOCATE(Var % EigenVectors(SIZE(Var % EigenValues), &
                    SIZE(Solver % Variable % EigenVectors,2)), STAT=AllocStat)
              END IF

              DO i=1,SIZE(Var % EigenVectors,2)
                SELECT CASE(k)
                CASE(1)
                  Var % EigenVectors(:,i) = &
                      CMPLX(REAL(Solver % Variable % EigenVectors(:,i)), 0.0_dp, kind=dp)
                CASE(2)
                  ! This is the imaginary component as a real-valued array: 
                  Var % EigenVectors(:,i) = &
                      CMPLX(AIMAG(Solver % Variable % EigenVectors(:,i)), 0.0_dp, kind=dp)
                END SELECT
              END DO
            ELSE
              CALL Warn(Caller, 'A complex-valued system should have 2 DOFs')
            END IF
          END IF
        END DO
      END IF
      
      ! This is temporal (?) fix for a glitch where the complex eigen vector
      ! is expanded to one where real and complex parts follow each other. 
      IF( ListGetLogical( Solver % Values,'Expand Eigen Vectors', Stat ) ) THEN
        CALL ExpandEigenVectors( A, Solver % Variable % EigenVectors, &
            Solver % NOFEigenValues, Solver % Variable % dofs )
      END IF
        
      CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
          Solver % Variable % Name )
    END IF
    
    ! We have solved {harmonic,eigen,constraint} system and no need to continue further
    IF( RecursiveAnalysis ) THEN
      IF( HarmonicAnalysis ) CALL ListAddLogical( Solver % Values,'Harmonic Analysis',.TRUE.)
      IF( EigenAnalysis ) CALL ListAddLogical( Solver % Values,'Eigen Analysis',.TRUE.)
      RETURN
    END IF


! Check whether b=0 since then equation Ax=b has only the trivial solution, x=0. 
! In case of a limiter one still may need to check the limiter for contact.
!-----------------------------------------------------------------------------
    IF( Parallel ) THEN
      bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))
    ELSE
      bnorm = SQRT(SUM(b(1:n)**2))
    END IF

    IF ( bnorm <= TINY( bnorm) .AND..NOT.SkipZeroRhs) THEN
      CALL Info(Caller,'Solution trivially zero!',Level=5)
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

    AndersonAcc = ListGetLogical( Params,'Nonlinear System Acceleration',GotIt ) 
    AndersonScaled = ListgetLogical( Params,'Nonlinear System Acceleration Scaled',GotIt ) 
    
    IF( AndersonAcc .AND. .NOT. AndersonScaled ) THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .TRUE., NoSolve )
      IF(NoSolve) GOTO 120
    END IF
    
!   Convert rhs & initial value to the scaled system:
!   -------------------------------------------------
    IF ( ScaleSystem ) THEN
      CALL ScaleLinearSystem(Solver, A, b, x, &
          RhsScaling = (bnorm/=0._dp), ConstraintScaling=.TRUE. )
    END IF

    ComputeChangeScaled = ListGetLogical(Params,&
        'Nonlinear System Compute Change in Scaled System',GotIt)
    IF(.NOT.GotIt) ComputeChangeScaled = .FALSE.

    IF(ComputeChangeScaled) THEN
       ALLOCATE(NonlinVals(SIZE(x)))
       NonlinVals = x
       IF (ASSOCIATED(Solver % Variable % Perm)) & 
           CALL RotateNTSystemAll(NonlinVals, Solver % Variable % Perm, DOFs)
    END IF

    IF( AndersonAcc .AND. AndersonScaled ) THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .TRUE., NoSolve )
      IF( NoSolve ) GOTO 110
    END IF
    
    IF( ListGetLogical( Params,'Linear System Normalize Guess',GotIt ) ) THEN
      CALL NormalizeInitialGuess() 
    END IF

    IF( ListGetLogical( Params,'Linear System Nullify Guess',GotIt ) ) THEN
      CALL Info(Caller,'Nullifying initial guess!',Level=30)
      x(1:n) = 0.0_dp
    END IF
    
    Method = ListGetString(Params,'Linear System Solver',GotIt)
    IF(.NOT. GotIt) THEN
      CALL Fatal(Caller,'Give "Linear System Solver", e.g. "iterative" or "direct"')
    END IF
    
    IF (Method=='multigrid' .OR. Method=='iterative' ) THEN
      Prec = ListGetString(Params,'Linear System Preconditioning',GotIt)
      IF( GotIt ) THEN
        CALL Info(Caller,'Linear System Preconditioning: '//TRIM(Prec),Level=8)
        CALL ResetTimer("Prec0-"//TRIM(Prec))
        IF( SEQL(Prec,'vanka') ) THEN
          IF(LEN(Prec)>=6) THEN
            i = ICHAR(Prec(6:6)) - ICHAR('0')
            CALL ListAddNewInteger( Params,'Vanka Mode',i) 
          END IF
          CALL VankaCreate(A,Solver)
        ELSE IF ( Prec=='circuit' ) THEN
          CALL CircuitPrecCreate(A,Solver)
#if 0 
          IF( ListGetLogical(Params,'Linear System Save', GotIt) ) THEN
            IF( ASSOCIATED( A % CircuitMatrix ) ) THEN
              CALL SaveLinearSystem(Solver, Ain = A % CircuitMatrix, LinSysName = "circuit")
            END IF
            IF( ASSOCIATED( Solver % Matrix % AddMatrix ) ) THEN
              CALL SaveLinearSystem(Solver, Ain = Solver % Matrix % AddMatrix, LinSysName = "addmatrix")
            END IF
          END IF
#endif
        END IF
        CALL CheckTimer("Prec0-"//TRIM(Prec),Level=8,Delete=.TRUE.)                  
      END IF
    END IF

    IF (Method=='iterative' ) THEN
      k = ListGetInteger(Params,'Linear System Dense Dofs',GotIt)
      IF(GotIt) A % ndeg = k

      IF( A % ndeg > 1 ) THEN
        IF( CRS_CheckStructuredDofs( A, A % ndeg) ) THEN
          CALL Fatal(Caller,'CRS matrix failed the dense test of size '//I2S(A % ndeg))
          A % ndeg = 0
        ELSE
          CALL Info(Caller,'CRS matrix passed the dense test of size '//I2S(A % ndeg),Level=12)
        END IF
      END IF
    END IF
     
    IF( InfoActive(30) ) THEN
      CALL VectorValuesRange(A % values,SIZE(A % values),'A')
      pv => b
      CALL VectorValuesRange(pv,SIZE(pv),'b')
    END IF


    IF(ListGetLogical(Params, 'Linear System Use Rocalution', Found)) &
      Method = 'rocalution'
    
    IF ( .NOT. Parallel ) THEN
      IF(ListGetLogical(Params, 'Linear System Use Hypre', Found)) Method = 'hypre'

      CALL Info(Caller,'Serial linear System Solver: '//TRIM(Method),Level=8)
      
      SELECT CASE(Method)
      CASE('multigrid')
        CALL MultiGridSolve( A, x, b, &
            DOFs, Solver, Solver % MultiGridLevel )
      CASE('iterative')
        CALL IterSolver( A, x, b, Solver )
      CASE('feti')
        CALL Fatal(Caller,'Feti solver available only in parallel.')
      CASE('block')
        CALL BlockSolveExt( A, x, b, Solver )
      CASE('amgx')
        CALL AMGXSolver( A, x, b, Solver )
      CASE('rocalution')
        CALL ROCSolver( A, x, b, Solver )
      CASE('hypre')
        CALL SolveHypre( A, x, b, Solver )
      CASE('direct')
        CALL DirectSolver( A, x, b, Solver )        
      CASE DEFAULT        
        CALL Fatal(Caller,'Unknown "Linear System Solver": '//TRIM(Method))
      END SELECT
    ELSE
      CALL Info(Caller,'Parallel linear System Solver: '//TRIM(Method),Level=8)

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
      CASE('amgx')
        CALL AMGXSolver( A, x, b, Solver )
      CASE('rocalution')
        CALL ROCSolver( A, x, b, Solver )
      CASE('direct')
        CALL DirectSolver( A, x, b, Solver )
      CASE DEFAULT        
        CALL Fatal(Caller,'Unknown "Linear System Solver": '//TRIM(Method))
      END SELECT
    END IF

    IF( InfoActive(30) ) THEN
      pv => x
      CALL VectorValuesRange(pv,SIZE(pv),'x')
    END IF
    
110 IF( AndersonAcc .AND. AndersonScaled )  THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .FALSE.)
    END IF

    IF(ComputeChangeScaled) THEN
      CALL ComputeChange(Solver,.FALSE.,n, x, NonlinVals, Matrix=A, RHS=b )
      DEALLOCATE(NonlinVals)
    END IF

    IF ( ScaleSystem ) THEN
      CALL BackScaleLinearSystem( Solver, A, b, x, ConstraintScaling=.TRUE. )
    END IF

120 IF( AndersonAcc .AND. .NOT. AndersonScaled )  THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .FALSE.)
    END IF
    
    Aaid => A
    IF (PRESENT(BulkMatrix)) THEN
      IF (ASSOCIATED(BulkMatrix) ) Aaid=>BulkMatrix
    END IF

    NodalLoads => NULL()
    IF(.NOT. SkipLoads ) THEN
      NodalLoads => VariableGet( Solver % Mesh % Variables, &
          GetVarName(Solver % Variable) // ' Loads' )
      IF( ASSOCIATED( NodalLoads ) ) THEN
        ! Nodal loads may be allocated but the user may have toggled
        ! the 'calculate loads' flag such that no load computation should be performed.
        CalcLoads = ListGetLogical( Solver % Values,'Calculate Loads',GotIt )
        IF( .NOT. GotIt ) CalcLoads = .TRUE.
        IF( CalcLoads ) THEN
          CALL Info(Caller,'Calculating nodal loads for: '//&
              GetVarName(Solver % Variable),Level=6)
          CALL CalculateLoads( Solver, Aaid, x, Dofs, .TRUE., NodalLoads ) 
        END IF
      END IF
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

    IF(.NOT. SkipLoads ) THEN
      NodalLoads => VariableGet( Solver % Mesh % Variables, &
          GetVarName(Solver % Variable) // ' Residual' )
      IF( ASSOCIATED( NodalLoads ) ) THEN
        CalcLoads = ListGetLogical( Solver % Values,'Calculate Residual',GotIt )
        IF( .NOT. GotIt ) CalcLoads = .TRUE.
        IF( CalcLoads ) THEN
          CALL Info(Caller,'Calculating nodal residual',Level=6)
          CALL CalculateLoads( Solver, Aaid, x, Dofs, .FALSE., NodalLoads ) 
        END IF
      END IF
    END IF
    
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
! In order to be able to change the preconditioners or solvers the old matrix structures
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

 CONTAINS


   ! Sometimes the r.h.s. may abruptly diminish in value resulting to significant 
   ! convergence issues or it may be that the system scales linearly with the source. 
   ! This flag tries to improve on the initial guess of the linear solvers, and may 
   ! sometimes even result to the exact solution.
   !--------------------------------------------------------------------------------
   SUBROUTINE NormalizeInitialGuess()
     REAL(KIND=dp) :: xn, bn
     REAL(KIND=dp), ALLOCATABLE, TARGET :: TempVector(:)

     
     CALL Info(Caller,'Normalizing initial guess!',Level=30)

     ALLOCATE( TempVector(A % NumberOfRows) )

     IF ( Parallel ) THEN
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
       CALL Info(Caller,Message,Level=6) 
     END IF
     DEALLOCATE( TempVector )

   END SUBROUTINE NormalizeInitialGuess


   
!------------------------------------------------------------------------------
  END SUBROUTINE SolveLinearSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AMGXMatrixVectorMultiply( A, u, v, Solver )
!------------------------------------------------------------------------------
    USE iso_c_binding, only: C_INTPTR_T, C_CHAR, C_NULL_CHAR

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: u(*), v(*)

#ifdef HAVE_AMGX
    INTERFACE
      SUBROUTINE AMGXmv(AMGX, n, rows, cols, vals, b, x, nonlin_update ) BIND(C, Name="AMGXmv")

         USE Types
         USE ISO_C_BINDING, ONLY: C_CHAR, C_INTPTR_T

         IMPLICIT NONE

         INTEGER(KIND=C_INTPTR_T) :: AMGX
         REAL(KIND=dp) :: vals(*), b(*), x(*)
         INTEGER :: rows(*), cols(*), nonlin_update, n
      END SUBROUTINE AMGXmv
    END INTERFACE

    LOGICAL :: Found
    INTEGER :: nonlin_update, i

    nonlin_update = 0
    CALL AMGXmv( A % AMGXMV, A % NumberOfRows, A % Rows-1, A % Cols-1, &
                  A % Values, u, v, nonlin_update )
#else
    CALL Fatal('AMGXSolver', "AMGX doesn't seem to be included.")
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE AMGXMatrixVectorMultiply
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AMGXSolver( A, x, b, Solver )
!------------------------------------------------------------------------------
    USE iso_c_binding, only: C_INTPTR_T, C_CHAR, C_NULL_CHAR

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: x(:), b(:)

#ifdef HAVE_AMGX
    INTERFACE
      SUBROUTINE AMGXSolve(AMGX, n, rows, cols, vals, b, x, &
              nonlin_update, config, comm, ng, part_vec, bnrm ) BIND(C, Name="AMGXSolve")

         USE Types
         USE ISO_C_BINDING, ONLY: C_CHAR, C_INTPTR_T

         IMPLICIT NONE

         INTEGER(KIND=C_INTPTR_T) :: AMGX
         REAL(KIND=dp) :: vals(*), b(*), x(*), bnrm
         CHARACTER(KIND=C_CHAR) :: config(*)
         INTEGER :: rows(*), cols(*), nonlin_update, n, comm, ng, part_vec(*)
      END SUBROUTINE AMGXSolve
    END INTERFACE


    CHARACTER(KIND=C_CHAR) :: cfg(MAX_NAME_LEN)
    CHARACTER(LEN=MAX_NAME_LEN) :: config
    LOGICAL :: found, isparallel , nlin
    INTEGER :: nonlin_update, i, j, n, lrow, me, you
    REAL(KIND=dp)  :: bnrm

    TYPE(Matrix_t), POINTER :: Bm
    INTEGER, SAVE :: ng
    INTEGER, ALLOCATABLE, SAVE :: Owner(:), APerm(:), part_vec(:), iLPerm(:)

    nlin = ListGetLogical( Solver % Values, 'Linear System Refactorize', Found )
    IF ( nlin .OR. .NOT. Found ) THEN
      nonlin_update = 1
    ELSE
      nonlin_update = 0
    END IF

    config = ListGetString( Solver % Values, 'AMGX Config')
    DO i=1,LEN_TRIM(config)
      cfg(i) = config(i:i)
    END DO
    cfg(i) = C_NULL_CHAR

    me =  Parenv % MyPe
    isParallel = Parenv % PEs>1

    IF(isParallel) THEN

      BLOCK
        INTEGER, ALLOCATABLE :: part_vec_tmp(:)

        INTEGER, ALLOCATABLE, SAVE :: GlobalToLocal(:),rRows(:),rSize(:),cBuf(:),SendTo(:)
        REAL(KIND=dp), ALLOCATABLE :: vBuf(:)
        INTEGER :: status(MPI_STATUS_SIZE),ierr,i,j,k,l,lrow,you,rcnt,proc

        TYPE SendStuff_t
          INTEGER, ALLOCATABLE :: Size(:), Rows(:)
        END TYPE SendStuff_t
        TYPE(SendStuff_t), ALLOCATABLE :: SendStuff(:)
 
        IF(.NOT. ASSOCIATED(A % CollectionMatrix)) THEN
          IF(ALLOCATED(Owner))THEN
            DEALLOCATE(Owner,APerm,ILperm,part_vec, SendTo, GlobalToLocal)
          END IF

          n = SIZE(A % ParallelInfo % GlobalDOFs)
          ALLOCATE( Owner(n), Aperm(n) )
          CALL ContinuousNumbering(A % ParallelInfo,A % Perm,APerm,Owner)

          ng = 0
          DO i=1, A % NumberOfRows
            IF(A % ParallelInfo % NeighbourList(i) % Neighbours(1) == me) ng=ng+1
          END DO
          ng = NINT(ParallelReduction(1._dp*ng))

          ALLOCATE(part_vec(ng), part_vec_tmp(ng));
          part_vec_tmp=-1; part_vec=-1
          DO i=1,A % NumberOfRows
            part_vec_tmp(APerm(i)) = A % ParallelInfo % NeighbourList(i) % Neighbours(1)
          END DO
          CALL MPI_ALLREDUCE( part_vec_tmp, part_vec, ng, MPI_INTEGER, MPI_MAX, &
                            ELMER_COMM_WORLD, ierr )
          DEALLOCATE(part_vec_tmp)

          Bm => AllocateMatrix(); Bm % Format = MATRIX_LIST
          ALLOCATE(SendTo(ParEnv % Pes),GlobalToLocal(ng),iLPerm(A % NumberOfRows))
        ELSE
          Bm => A % CollectionMatrix
        END IF

        IF( Bm % Format == MATRIX_LIST .OR. nonlin_update==1 ) THEN
          iLPerm = 0
          LRow = 0
          SendTo = 0
          GlobalToLocal = 0
 
          IF (Bm % Format == MATRIX_CRS  ) Bm % Values = 0._dp
          DO i=1,A % NumberofRows
            you = A % ParallelInfo % NeighbourList(i) % Neighbours(1)
            IF ( you == me ) THEN
              LRow = LRow+1
              iLPerm(LRow) = i
              DO j=A % Rows(i+1)-1, A % Rows(i),-1
                CALL AddToMatrixElement(Bm,LRow,APerm(A  % Cols(j)), A % Values(j))
              END DO
              GlobalToLocal(APerm(i)) = LRow
            ELSE
              SendTo(you+1) = SendTo(you+1)+1
            END IF
          END DO

          ALLOCATE(SendStuff(ParEnv % Pes))
          DO i=1,ParEnv % PEs
            IF( i-1==me ) CYCLE
            IF(.NOT.ParEnv % IsNeighbour(i))  CYCLE
  
            ALLOCATE( SendStuff(i) % Rows(SendTo(i)) )
            ALLOCATE( SendStuff(i) % Size(SendTo(i)) )
          END DO
 
          SendTo = 0
          DO i=1,a % NumberOfRows
            you = A % ParallelInfo % NeighbourList(i) % Neighbours(1)
            IF ( you /= me ) THEN
              SendTo(you+1) = SendTo(you+1)+1
              SendStuff(you+1) % Size(Sendto(you+1))  = A % Rows(i+1)-A % Rows(i)
              SendStuff(you+1) % Rows(Sendto(you+1))  = i
            END IF
          END DO

          DO i=1,ParEnV % PEs
            IF(i-1==me .OR. .NOT. ParEnv % IsNeighbour(i)) CYCLE

            CALL MPI_BSEND(SendTo(i),1,MPI_INTEGER,i-1,1200,ELMER_COMM_WORLD, ierr)
            IF(Sendto(i)==0) CYCLE

            CALL MPI_BSEND(APerm(SendStuff(i) % Rows),SendTo(i),MPI_INTEGER,i-1, &
                          1201,ELMER_COMM_WORLD,ierr )
  
            CALL MPI_BSEND( SendStuff(i) % Size,SendTo(i),MPI_INTEGER,i-1, &
                          1202,ELMER_COMM_WORLD,ierr )
            DO j=1,SendTo(i)
              k = SendStuff(i) % Rows(j)
              CALL MPI_BSEND(APerm(A % Cols(A % Rows(k):A % Rows(k+1)-1)),SendStuff(i) % Size(j), &
                         MPI_INTEGER,i-1, 1203,ELMER_COMM_WORLD, ierr )

              CALL MPI_BSEND(A % Values(A % Rows(k):A % Rows(k+1)-1),SendStuff(i) % Size(j), &
                      MPI_DOUBLE_PRECISION,i-1,1204,ELMER_COMM_WORLD, ierr )
            END DO
          END DO

          DO i=1,ParEnV % PEs
            IF(i-1==me .OR. .NOT. ParEnv % IsNeighbour(i)) CYCLE

            CALL MPI_RECV(rcnt,1,MPI_INTEGER,MPI_ANY_SOURCE,1200,ELMER_COMM_WORLD,status,ierr)
            IF(rcnt==0) CYCLE

            proc = status(MPI_SOURCE)
            ALLOCATE( rRows(rcnt), rSize(rcnt) )

            CALL MPI_RECV(rRows,rcnt,MPI_INTEGER,proc,1201,ELMER_COMM_WORLD,status,ierr)
            CALL MPI_RECV(rSize,rcnt,MPI_INTEGER,proc,1202,ELMER_COMM_WORLD,status,ierr)
            DO j=1,rcnt
              k = GlobalToLocal(rRows(j))
              IF ( k==0 ) THEN
                PRINT*,Parenv % MyPE,proc, 'not mine then ?', rRows(j)
                CYCLE
              END IF
              ALLOCATE(cBuf(rSize(j)), vBuf(rSize(j)))
              CALL MPI_RECV(cBuf,rSize(j),MPI_INTEGER,proc,1203, &
                        ELMER_COMM_WORLD,status,ierr )

              CALL MPI_RECV(vBuf,rSize(j),MPI_DOUBLE_PRECISION,proc, &
                      1204,ELMER_COMM_WORLD, status,ierr)

              DO l=1,rSize(j)
                CAll AddToMatrixElement(Bm,k,cBuf(l),vBuf(l))
              END DO
              DEALLOCATE(cbuf, vbuf)
            END DO

            DEALLOCATE( rRows, rSize )
          END DO
          CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)
 
          IF( Bm % Format == MATRIX_LIST ) THEN
            CALL List_toCRSMatrix(Bm)
            A % CollectionMatrix => Bm
          END IF
        END IF
      END BLOCK

      BLOCK
        REAL(KIND=dp), ALLOCATABLE :: bb(:),xb(:), r(:)

        n = Bm % NumberOfRows
        j = A  % NumberOfRows
        ALLOCATE(bb(n), xb(n), r(j) )

        r = b(1:j)
        CALL ParallelSumVector(A, r)

        DO i=1,n
          bb(i) = r(iLPerm(i))
          xb(i) = x(iLPerm(i))
        END DO
        bnrm = SQRT(ParallelReduction(SUM(bb**2)))

        CALL AMGXSolve( A % AMGX, n, Bm % Rows-1, Bm % Cols-1, Bm % Values,  &
          bb, xb, nonlin_update, cfg, ELMER_COMM_WORLD, ng, part_vec, bnrm )

        x = 0
        DO i=1,n
          x(iLPerm(i)) = xb(i)
        END DO
        CALL ParallelSumVector(A, x)
      END BLOCK
    ELSE
      n = A % NumberOfRows
      bnrm = SQRT(SUM(b**2))

      CALL AMGXSolve( A % AMGX, n, A % Rows-1, A % Cols-1, &
        A % Values, b, x, nonlin_update, cfg, ELMER_COMM_WORLD, n, &
           A % Diag, bnrm ) ! <--- a % diag  for dummy
    END IF
#else
    CALL Fatal('AMGXSolver', "AMGX doesn't seem to be included.")
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE AMGXSolver
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE ROCSolver( A, x, b, Solver )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: x(:), b(:)

#ifdef HAVE_ROCALUTION

    ! Interfaces to c++ code calling ROCalution library
    ! -------------------------------------------------
    INTERFACE
      SUBROUTINE ROCSerialSolve(n, rows, cols, vals, b, x, nonlin_update, &
                  imethod, prec, maxiter, tol) BIND(C, Name="ROCSerialSolve")
        USE Types
        USE ISO_C_BINDING, ONLY: C_CHAR, C_INTPTR_T

        IMPLICIT NONE
        REAL(KIND=dp) :: vals(*), b(*), x(*), tol
        INTEGER :: rows(*), cols(*), nonlin_update, n, imethod, prec, maxiter
      END SUBROUTINE ROCSerialSolve


      SUBROUTINE ROCParallelSolve(gn, n, rows, cols, vals, b, x, bnrm, goffset, &
         fcomm, imethod, prec, maxiter, tol) BIND(C, Name="ROCParallelSolve")
        USE Types
        USE ISO_C_BINDING, ONLY: C_CHAR, C_INTPTR_T

        IMPLICIT NONE
        REAL(KIND=dp) :: vals(*), b(*), x(*), bnrm, tol
        INTEGER :: gn, n, rows(*), cols(*), goffset(*), fcomm, imethod, prec, maxiter
      END SUBROUTINE ROCParallelSolve
    END INTERFACE

    ! local variables:
    ! ----------------
    LOGICAL :: found, isParallel 
    INTEGER :: nonlin_update, i, j, k,l,n, gn, me

    TYPE(Matrix_t), POINTER ::Rmatrix

    INTEGER, POINTER ::  aPerm(:), iLperm(:), gOffset(:)
    REAL(KIND=dp), ALLOCATABLE :: dBuf(:)
    INTEGER, ALLOCATABLE :: Owner(:), SendTo(:), iBuf(:), tOffset(:), rRows(:), rSize(:)

    INTEGER :: status(MPI_STATUS_SIZE),ierr,lrow,you,rcnt,proc

    INTEGER :: buf_size, procs

    ! Define these to get somewhat shorter MPI subroutine calls:
    ! -----------------------------------------------------------
    INTEGER :: xmpi_comm
    INTEGER :: xmpi_src  = MPI_ANY_SOURCE
    INTEGER :: xmpi_max  = MPI_MAX
    INTEGER :: xmpi_int  = MPI_INTEGER
    INTEGER :: xmpi_dbl  = MPI_DOUBLE_PRECISION


    REAL(KIND=dp) :: TOL
    INTEGER :: Imethod, Prec, MaxIter, ILULevel, own_n

    TYPE(ValueList_t), POINTER :: Params

    TYPE SendStuff_t
      INTEGER, ALLOCATABLE :: Size(:), Rows(:)
    END TYPE SendStuff_t
    TYPE(SendStuff_t), ALLOCATABLE :: SendStuff(:)

    Params => Solver % Values

    ! Extract some controls to ROCalution from the simulation control info:
    ! ---------------------------------------------------------------------
    nonlin_update = 1
    IF ( .NOT. ListGetLogical( Params, 'Linear System Refactorize', Found ) ) &
      nonlin_update = 0;

    SELECT CASE(ListGetString(Params,'Linear System Iterative Method',Found))
      CASE('cg')
         Imethod = 0;
      CASE('bicgstab')
         Imethod = 1;
      CASE('bicgstabl')
         Imethod = 2;
      CASE('gmres')
         Imethod = 3;
      CASE('fgmres')
         Imethod = 4;
      CASE DEFAULT
         Imethod = 0;
    END SELECT

    SELECT CASE(ListGetString(Params,'Linear System Preconditioning',Found))
      CASE('jacobi')
        Prec = 0;
      CASE('sgs')
        Prec = 1;
      CASE('ilu')
        Prec = 2; ILULevel = 0
      CASE('ilu0')
       Prec = 2; ILULevel = 0
      CASE('ilu1')
       Prec = 2; ILULevel = 1
      CASE('ilu2')
       Prec = 2; ILULevel = 2
      CASE DEFAULT
       Prec = 0;
    END SELECT

    MaxIter = ListGetInteger(Params,'Linear System Max Iterations')
    TOL = ListGetCReal(Params,'Linear System Convergence Tolerance')

    ! ---------------------------------------------------------------------

    n = A % NumberOfRows

    procs =  Parenv % PEs
    isParallel = procs>1

    IF(isParallel) THEN
      me    =  Parenv % MyPe
      xmpi_comm = ELMER_COMM_WORLD

      IF (.NOT.ASSOCIATED(A % ParMatrix)) CALL ParallelInitMatrix(Solver,A)

      ParEnv => A % ParMatrix % ParEnv
      ParEnv % ActiveComm = A % Comm

      ! Enforce continuous ascending numbering as required by ROCalution 
      ! ----------------------------------------------------------------
      IF (.NOT. ASSOCIATED(A % RocParams % Rmatrix)) THEN
        n = SIZE(A % ParallelInfo % GlobalDOFs)
        ALLOCATE( Owner(n), Aperm(n), tOffset(0:procs), gOffset(0:procs) )

        tOffset = -1
        CALL ContinuousNumbering(A % ParallelInfo,A % Perm,aPerm,Owner,gStart=tOffset(me))

        gOffset = -1
        CALL MPI_ALLREDUCE(tOffset,gOffset,procs,xmpi_int,xmpi_max,xmpi_comm,ierr)

        own_n = SUM(Owner)
        gOffset(procs) = ParallelReduction(own_n)

        ALLOCATE(iLPerm(A % NumberOfRows))

        Rmatrix => AllocateMatrix();
        Rmatrix % Format = MATRIX_LIST
        Rmatrix % ListMatrix => List_AllocateMatrix(own_n)

        A % RocParams % Rmatrix => Rmatrix
        A % RocParams % CntPerm => aPerm
        A % RocParams % LocPerm => iLperm
        A % RocParams % gOffset => gOffset
      ELSE
        Rmatrix => A % RocParams % Rmatrix
        aPerm   => A % RocParams % CntPerm
        iLPerm  => A % RocParams % LocPerm
        gOffset => A % RocParams % gOffset
      END IF

      ! Complete the matrix rows such that each partition has full rows of the 'owned' dofs
      ! -----------------------------------------------------------------------------------
      IF( Rmatrix % Format == MATRIX_LIST .OR. nonlin_update==1 ) THEN
 
        IF (Rmatrix % Format == MATRIX_CRS  ) Rmatrix % Values = 0._dp

        ! Create inside matrix + count rows with values to send for each neighbour
        ! -------------------------------------------------------------------------
        ALLOCATE(SendTo(procs))
        iLPerm = 0
        LRow = 0
        SendTo = 0
        DO i=1,A % NumberofRows
          you = A % ParallelInfo % NeighbourList(i) % Neighbours(1)
          IF ( you == me ) THEN
            lRow = lRow + 1
            iLPerm(lRow) = i
            DO j=A % Rows(i+1)-1, A % Rows(i),-1
              CALL AddToMatrixElement(Rmatrix, lRow, aPerm(A  % Cols(j)), A % Values(j))
            END DO
          ELSE
            SendTo(you+1) = SendTo(you+1)+1
          END IF
        END DO

        ALLOCATE(SendStuff(ParEnv % Pes))
        DO i=1,ParEnv % PEs
          IF( i-1==me ) CYCLE
          IF(.NOT.ParEnv % IsNeighbour(i))  CYCLE

          ALLOCATE( SendStuff(i) % Rows(SendTo(i)) )
          ALLOCATE( SendStuff(i) % Size(SendTo(i)) )
        END DO
 
        ! Count number of columns of each neighbour's rows to be sent
        ! -----------------------------------------------------------
        SendTo   = 0
        buf_size = 0
        DO i=1,a % NumberOfRows
          you = A % ParallelInfo % NeighbourList(i) % Neighbours(1)
          IF ( you /= me ) THEN
            SendTo(you+1) = SendTo(you+1)+1
            SendStuff(you+1) % Size(Sendto(you+1))  = A % Rows(i+1)-A % Rows(i)
            SendStuff(you+1) % Rows(Sendto(you+1))  = i
            buf_size = buf_size + A % Rows(i+1) - A % Rows(i)
          END IF
        END DO
        CALL CheckBuffer(ParEnv % PEs*(4+4*MPI_BSEND_OVERHEAD) + 4*buf_size)

        ! Send data to neighbours
        ! -----------------------
        DO i=1,ParEnV % PEs
          IF(i-1==me .OR. .NOT. ParEnv % IsNeighbour(i)) CYCLE

          CALL MPI_BSEND(SendTo(i),1,xmpi_int,i-1,1200,xmpi_comm,ierr)
          IF(Sendto(i)==0) CYCLE

          ibuf = aPerm(SendStuff(i) % Rows)
          CALL MPI_BSEND(ibuf,SendTo(i),xmpi_int,i-1,1201,xmpi_comm,ierr)

          ibuf = SendStuff(i) % Size
          CALL MPI_BSEND(ibuf,SendTo(i),xmpi_int,i-1,1202,xmpi_comm,ierr)

          DO j=1,SendTo(i)
            k = SendStuff(i) % Rows(j)
            l = SendStuff(i) % Size(j)
            dBuf =  A % Values(A % Rows(k):A % Rows(k+1)-1)
            iBuf =  aPerm(A % Cols(A % Rows(k):A % Rows(k+1)-1))
            CALL MPI_BSEND(iBuf,l,xmpi_int,i-1,1203,xmpi_comm,ierr)
            CALL MPI_BSEND(dBuf,l,xmpi_dbl,i-1,1204,xmpi_comm,ierr)
          END DO
        END DO

        ! receive data from neighbours
        ! ----------------------------
        ALLOCATE( rRows(100), rSize(100) )
        DO i=0,procs-1
          IF(i==me .OR. .NOT. ParEnv % IsNeighbour(i+1)) CYCLE

          CALL MPI_RECV(rcnt,1,xmpi_int,xmpi_src,1200,xmpi_comm,status,ierr)
          IF(rcnt==0) CYCLE

          IF(rcnt > SIZE(rRows)) THEN
            DEALLOCATE(rRows,Rsize)
            ALLOCATE( rRows(rcnt), rSize(rcnt) )
          END IF

          proc = status(MPI_SOURCE)
          CALL MPI_RECV(rRows,rcnt,xmpi_int,proc,1201,xmpi_comm,status,ierr)
          CALL MPI_RECV(rSize,rcnt,xmpi_int,proc,1202,xmpi_comm,status,ierr)
          DO j=1,rcnt
            k = rRows(j)

            IF ( k<= gOffset(me) .OR. k> gOffset(me+1) ) THEN
              PRINT*,Parenv % MyPE,proc, 'not mine then ?', rRows(j), gOffset(me), gOffset(me+1)
              CYCLE
            END IF

            IF(rSize(j) > SIZE(iBuf)) THEN
              DEALLOCATE(iBuf,dBuf)
              ALLOCATE( iBuf(rSize(j)), dBuf(rSize(j)) )
            END IF

            CALL MPI_RECV(iBuf,rSize(j),xmpi_int,proc,1203,xmpi_comm,status,ierr)
            CALL MPI_RECV(dBuf,rSize(j),xmpi_dbl,proc,1204,xmpi_comm,status,ierr)

            DO l=1,rSize(j)
              CAll AddToMatrixElement(Rmatrix,k-gOffset(me),iBuf(l),dBuf(l))
            END DO
          END DO
        END DO
        ! ----------

        CALL MPI_BARRIER(A % Comm,ierr)
 
        IF(Rmatrix % Format == MATRIX_LIST) CALL List_toCRSMatrix(Rmatrix)
        n = Rmatrix % NumberOfRows
        gn = ParallelReduction(n);
      END IF


      !  the linear solver
      ! ----------------------
      BLOCK
        REAL(KIND=dp), ALLOCATABLE :: pb(:),px(:), r(:)
        REAL(KIND=dp) :: bnrm

        n = Rmatrix % NumberOfRows
        j = A  % NumberOfRows
        ALLOCATE(pb(n), px(n), r(j) )

        ! complete the RHS of the partition dofs
        !---------------------------------------
        r = b(1:j)
        CALL ParallelSumVector(A, r)

        ! Extract 'owned' dofs
        !---------------------
        DO i=1,n
          pb(i) = r(iLPerm(i))
          px(i) = x(iLPerm(i))
        END DO

        bnrm = SQRT(ParallelReduction(SUM(pb**2)))
        IF(bnrm <AEPS) bnrm=1;

        !  call ROCalution
        ! ----------------
        CALL ROCParallelSolve( gn, n, Rmatrix % Rows-1, Rmatrix % Cols-1, &
             Rmatrix % Values, pb, px, bnrm, gOffset, A % comm, imethod, prec, maxiter, tol )

        ! distribute the result such that all dofs present (even shared and not 'owned')
        ! in a partition have consistent values
        ! ------------------------------------------------------------------------------
        x = 0
        DO i=1,n
          x(iLPerm(i)) = px(i)
        END DO
        CALL ParallelSumVector(A, x)
      END BLOCK

      ! Cleanup, remains to be reconsidered for optimizations
      ! -----------------------------------------------------
      CALL FreeMatrix(Rmatrix);
      DEALLOCATE(APerm,ILperm,gOffset)

      A % RocParams % Rmatrix => Null()
      A % RocParams % CntPerm => Null()
      A % RocParams % LocPerm => Null()
      A % RocParams % gOffset => Null()
    ELSE
      ! Serial case: call the linear solver
      ! -----------------------------------
      CALL ROCSerialSolve( n, A % Rows-1, A % Cols-1, A % Values, b, x, &
              nonlin_update, imethod, prec, maxiter, tol )
    END IF
#else
    CALL Fatal('ROCSolver', "Rocalution doesn't seem to be included.")
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE ROCSolver
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Given a linear system Ax=b make a change of variables such that we will 
!> be solving for the residual Adx=b-Ax0 where dx=x-x0.
!------------------------------------------------------------------------------
  SUBROUTINE LinearSystemResidual( A, b, x, r )
!------------------------------------------------------------------------------

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
!------------------------------------------------------------------------------
  END SUBROUTINE LinearSystemResidual
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Given a linear system Ax=b make a change of variables such that we will 
!> be solving for the residual Adx=b-Ax0 where dx=x-x0.
!------------------------------------------------------------------------------
  FUNCTION LinearSystemMaskedResidualNorm( A, b, x, ActiveRow, ActiveCol ) RESULT ( Nrm )

    REAL(KIND=dp) CONTIG :: b(:)   
    REAL(KIND=dp) CONTIG :: x(:)   
    TYPE(Matrix_t), POINTER :: A
    LOGICAL, DIMENSION(:) :: ActiveRow(:), ActiveCol(:)
    REAL(KIND=dp) :: Nrm
    
    REAL(KIND=dp), ALLOCATABLE :: r(:)
    INTEGER :: i,n,totn
    REAL(KIND=dp) :: r2sum

    n = A % NumberOfRows 

    ALLOCATE(r(n))
   
    IF (Parenv % Pes>1) THEN
      CALL Fatal('LinearSystemMaskedResidualNorm','Not implemented in parallel yet!')
!      CALL ParallelMatrixVector(A, x, r, .TRUE.)
    ELSE
      CALL MaskedMatrixVectorMultiply( A, x, r, ActiveRow, ActiveCol )
    END IF

    DO i=1,n
      IF( ActiveRow(i) ) THEN
        r(i) = b(i) - r(i)
      END IF
    END DO

    totn = ParallelReduction(n)

    r2sum = SUM( r**2 )
    Nrm = SQRT( ParallelReduction(r2sum) / totn )

    DEALLOCATE( r ) 
    
  END FUNCTION LinearSystemMaskedResidualNorm



  FUNCTION HaveRestrictionMatrix( A ) RESULT( HaveConstraint ) 

    TYPE(Matrix_t), POINTER :: A
    LOGICAL :: HaveConstraint

    INTEGER :: n

    IF( .NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal('HaveRestrictionMatrix','Matrix A not associated!')
    END IF

    n = 0
    IF ( ASSOCIATED(A % ConstraintMatrix) )  THEN
      IF ( A % ConstraintMatrix % NumberOFRows > 0 ) n = n + 1 
    END IF
         
    IF ( ASSOCIATED(A % AddMatrix) )  THEN
      IF ( A % AddMatrix % NumberOFRows > 0 ) n = n + 1
    END IF
    
    n = ParallelReduction(n)
    HaveConstraint = ( n > 0 ) 
    
  END FUNCTION HaveRestrictionMatrix

  
  
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
    LOGICAL :: Relax, Found, NeedPrevSol, Timing, ResidualMode, &
        RestrictionMode, BlockMode, GloNum, FirstLoop
    INTEGER :: n,i,j,k,l,m,istat,nrows,ncols,colsj,rowoffset
    CHARACTER(:), ALLOCATABLE :: Method, VariableName
    INTEGER(KIND=AddrInt) :: Proc
    REAL(KIND=dp) :: Relaxation,Alpha,Beta,Gamma
    REAL(KIND=dp), ALLOCATABLE :: Diag(:), TempVector(:)
    REAL(KIND=dp), POINTER :: bb(:),Res(:)
    REAL(KIND=dp) :: t0,rt0,rst,st,ct
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NMode, LinModes
    CHARACTER(*), PARAMETER :: Caller = 'SolveSystem'
    
    INTERFACE
      SUBROUTINE BlockSolveExt(A,x,b,Solver)
        USE Types
        TYPE(Matrix_t), POINTER :: A
        TYPE(Solver_t) :: Solver
        REAL(KIND=dp) :: x(:), b(:)
      END SUBROUTINE BlockSolveExt
    END INTERFACE   


!------------------------------------------------------------------------------
    Params => Solver % Values

    CALL Info(Caller,'Solving linear system',Level=10)

    Timing = ListCheckPrefix(Params,'Linear System Timing')
    IF( Timing ) THEN
      t0 = CPUTime(); rt0 = RealTime()
    END IF

    n = A % NumberOfRows

    RestrictionMode = HaveRestrictionMatrix( A ) 

    ResidualMode = ListGetLogical( Params,'Linear System Residual Mode',Found )      
      
    BlockMode = ListGetLogical( Params,'Linear System Block Mode',Found ) 
    
!------------------------------------------------------------------------------
! The allocation of previous values has to be here in order to 
! work properly with the Dirichlet elimination.
!------------------------------------------------------------------------------
    NeedPrevSol = ResidualMode
    
    IF(.NOT. NeedPrevSol ) THEN
      Relaxation = ListGetCReal( Params, &
          'Nonlinear System Relaxation Factor', Found )
      IF( Found ) NeedPrevSol = (Relaxation /= 1.0_dp)
    END IF

    IF(.NOT. NeedPrevSol ) THEN
      Method = ListGetString( Params, &
        'Nonlinear System Convergence Measure', Found ) 
      NeedPrevSol = ( Method == 'residual' .OR. Method == 'solution' )
    END IF

    IF( NeedPrevSol ) THEN
      CALL Info(Caller,'Previous solution must be stored before system is solved',Level=10)
      Found = ASSOCIATED(Solver % Variable % NonlinValues)
      IF( Found ) THEN
        IF ( SIZE(Solver % Variable % NonlinValues) /= n) THEN
          DEALLOCATE(Solver % Variable % NonlinValues)
          Found = .FALSE.
        END IF
      END IF
      IF(.NOT. Found) THEN
        ALLOCATE( Solver % Variable % NonlinValues(n), STAT=istat ) 
        IF ( istat /= 0 ) CALL Fatal( Caller, 'Memory allocation error.' )
      END IF
      Solver % Variable % NonlinValues = x(1:n)
    END IF

    IF ( Solver % LinBeforeProc /= 0 ) THEN
      CALL Info(Caller,'Calling procedure before solving system',Level=7)
      istat = ExecLinSolveProcs( Solver % LinBeforeProc,CurrentModel,Solver, &
                       A, b, x, n, DOFs, Norm )
       IF ( istat /= 0 ) GOTO 10
    END IF

    ! If residual mode is requested make change of variables:
    ! Ax=b -> Adx = b-Ax0 = r
    IF( ResidualMode ) THEN
      CALL Info(Caller,'Changing the equation to residual based mode',Level=10)
      ALLOCATE( Res(n) ) 

      ! If needed move the current solution to N-T coordinate system
      ! before computing the residual.
      IF (ASSOCIATED(Solver % Variable % Perm)) &
          CALL RotateNTSystemAll(x, Solver % Variable % Perm, DOFs)

      CALL LinearSystemResidual( A, b, x, res )
      bb => res
      ! Set the initial guess for the residual system to zero
      x = 0.0_dp
    ELSE
      bb => b
    END IF

    FirstLoop = .TRUE.
    Nmode = 0
20  CALL ConstraintModesDriver( A, x, b, Solver, .TRUE., Nmode, LinModes, FirstLoop = FirstLoop )  
    
    
    IF( BlockMode ) THEN
      CALL Info(Caller,'Solving linear system with block strategy',Level=10)
      ! Here activate constraint solve only if constraints are not treated as blocks
      IF( RestrictionMode .AND. &
          ListGetLogical( Params, 'Eliminate Linear Constraints', Found) ) THEN
        BLOCK 
          TYPE(Matrix_t), POINTER :: Acoll      
          Acoll  => AllocateMatrix()
          Acoll % FORMAT = MATRIX_LIST        
          CALL Info(Caller,'Eliminating constraints before going into block matrix!')
          CALL EliminateLinearRestriction( A, bb, A % ConstraintMatrix, Acoll, Solver, .TRUE. )
          CALL List_ToCRSMatrix(Acoll)

          Acoll % Comm = A % Comm 
          Acoll % AddMatrix => A % AddMatrix
          CALL ParallelInitMatrix(Solver, Acoll)
          
          CALL BlockSolveExt( Acoll, x, Acoll % rhs, Solver )

          CALL Info(Caller,'Freeing collection matrix after solution',Level=10)
          NULLIFY( Acoll % AddMatrix )         

          CALL FreeMatrix(Acoll)
          ParEnv => A % ParMatrix % ParEnv

          Acoll => NULL()
        END BLOCK
      ELSE
        CALL BlockSolveExt( A, x, bb, Solver )
      END IF
    ELSE IF ( RestrictionMode ) THEN
      CALL Info(Caller,'Solving linear system with linear restrictions!',Level=10)
      IF( ListGetLogical( Params,'Save Constraint Matrix',Found ) ) THEN
        GloNum = ListGetLogical( Params,'Save Constraint Matrix Global Numbering',Found )
        CALL SaveProjector(A % ConstraintMatrix,.TRUE.,'cm',Parallel=GloNum)
      END IF
      CALL SolveWithLinearRestriction( A,bb,x,Norm,DOFs,Solver )
    ELSE ! standard mode
      CALL Info(Caller,'Solving linear system in standard way',Level=12)
      CALL SolveLinearSystem( A,bb,x,Norm,DOFs,Solver )
    END IF
    CALL Info(Caller,'System solved',Level=12)

    
    IF( LinModes > 0 .OR. Nmode > 0 ) THEN
      CALL ConstraintModesDriver( A, x, b, Solver, .FALSE., FirstLoop = FirstLoop ) 
      FirstLoop = .FALSE.
      IF( Nmode < LinModes ) GOTO 20
    END IF
    
    ! Even in the residual mode the system is reverted back to complete vectors 
    ! and we may forget about the residual.
    IF( ResidualMode ) DEALLOCATE( Res ) 

!------------------------------------------------------------------------------

10  CONTINUE

    IF ( Solver % LinAfterProc /= 0 ) THEN
      CALL Info(Caller,'Calling procedure after solving system',Level=7)
      istat = ExecLinSolveProcs( Solver % LinAfterProc, CurrentModel, Solver, &
              A, b, x, n, DOFs, Norm )
    END IF

    IF ( Solver % TimeOrder == 2 ) THEN
      CALL Info(Caller,'Setting up PrevValues for 2nd order transient equations',Level=12)

      IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        CALL Update2ndOrder(n,Solver % dt,x, &
            Solver % Variable % PrevValues,Solver % Alpha,Solver % Beta)
      END IF
    END IF

    IF( Timing ) THEN
      st  = CPUTime() - t0;
      rst = RealTime() - rt0

      WRITE(Message,'(a,f8.2,f8.2,a)') 'Linear system time (CPU,REAL) for '&
          //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
      CALL Info(Caller,Message,Level=4)    
      
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

        WRITE(Message,'(a,f8.2,f8.2,a)') 'Linear system time cumulative (CPU,REAL) for '&
            //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
        CALL Info(Caller,Message,Level=7)    
        
        CALL ListAddConstReal(CurrentModel % Simulation,'res: cum linsys cpu time '&
            //GetVarName(Solver % Variable),st)
        CALL ListAddConstReal(CurrentModel % Simulation,'res: cum linsys real time '&
            //GetVarName(Solver % Variable),rst)
      END IF

    END IF

    CALL Info(Caller,'Finished solving the system',Level=12)

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
    COMPLEX(KIND=dp), TARGET :: EigenValues(:),EigenVectors(:,:)
    REAL(KIND=dp) :: Norm
    TYPE(Matrix_t), POINTER :: StiffMatrix
    INTEGER :: NOFEigen
    TYPE(Solver_t) :: Solver
    !------------------------------------------------------------------------------
    COMPLEX(KIND=dp), POINTER :: dvecs(:,:), evecs(:,:)
    REAL(KIND=dp), POINTER :: p(:)
    INTEGER :: i,j,k,n, AllocStat
    TYPE(Matrix_t), POINTER :: A
    TYPE(Variable_t), POINTER :: Var 
    LOGICAL :: Damped, Direct, Found
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    !------------------------------------------------------------------------------
    n = StiffMatrix % NumberOfRows
    EVecs => EigenVectors

    A => StiffMatrix

    Damped = ListGetLogical( Solver % Values, 'Eigen System Damped', Found)
    Direct = ListGetString( Solver % Values, 'Linear System Solver', Found) == 'direct'

    IF( Damped .AND. Direct ) THEN
      A => GenerateStateEquationSystem(Solver,A,n)
      ALLOCATE(EVecs(NOFEigen,n)); EVecs=0
      CALL ListAddLogical( Solver % Values,'Eigen System Damped', .FALSE. )
    END IF

    IF ( .NOT. A % Complex ) THEN
      CALL Info('SolveEigenSystem','Solving real valued eigen system of size: '//I2S(n),Level=8)

      IF ( ParEnv % PEs <= 1 ) THEN
        CALL ArpackEigenSolve( Solver, A, n, NOFEigen, EigenValues, evecs )
      ELSE
        CALL ParallelArpackEigenSolve( Solver, A, n, NOFEigen, EigenValues, evecs )
      END IF
    ELSE
      CALL Info('SolveEigenSystem','Solving complex valued eigen system of size: '//I2S(n/2),Level=8)
      IF ( ParEnv % PEs <= 1 ) THEN
        CALL ArpackEigenSolveComplex( Solver, A, n/2, NOFEigen, EigenValues, evecs )
      ELSE
        CALL ParallelArpackEigenSolveComplex( Solver, A, n/2, NOFEigen, EigenValues, EVecs )
      END IF
    END IF

    IF (Damped.AND.Direct) THEN
      n = n/2
      DO i=1,NOFEigen
        DO j=1,n
          EigenVectors(i,j) = EVecs(i,j+n)
        END DO
      END DO
      DEALLOCATE(EVecs)
      CALL ListAddLogical( Solver % Values,'Eigen System Damped', .TRUE. )
    END IF

 CONTAINS

     FUNCTION GenerateStateEquationSystem(Solver,A,n) RESULT(B)

       TYPE(Matrix_t), POINTER :: A,B
       TYPE(Solver_t):: Solver

       INTEGER :: i,j,k,l,n, max_gdofs

! Make the system matrices for the damped eigensystem

!      (M 0)           lambda *  ( C K)
!      (0 M) [y x]^T =           (-M 0) [y x]^T

       B => AllocateMatrix()
       B % Format = MATRIX_CRS
       B % Complex = A % Complex

       B % NumberOfRows = 2*n
       l = 4*SIZE(A % Cols)
       ALLOCATE(B % Cols(l), B % Rows(2*n+1), B % Values(l), B % MassValues(l))
       B % Values = 0; B % MassValues = 0

       B % Rows(1) = 1
       l = 0
       DO i=1,n
         DO k=A % Rows(i), A % Rows(i+1)-1
           l = l + 1
           B % Cols(l) = A % Cols(k)
           B % Values(l) = A % DampValues(k)
           B % MassValues(l) = A % MassValues(k)
         END DO
         DO k=A % Rows(i), A % Rows(i+1)-1
           l = l + 1
           B % Cols(l) = A % Cols(k)+n
           B % Values(l) = A % Values(k)
         END DO
         B % Rows(i+1) = l+1
       END DO

       DO i=1,n
         DO k=A % Rows(i), A % Rows(i+1)-1
           l = l + 1
           B % Cols(l) = A % Cols(k)
           B % Values(l) = -A % MassValues(k)
         END DO
 
         DO k=A % Rows(i), A % Rows(i+1)-1
           l = l + 1
           B % Cols(l) = A % Cols(k)+n
           B % MassValues(l) = A % MassValues(k)
         END DO
         B % Rows(n+i+1) = l+1
       END DO

       IF(ParEnv % PEs>1) THEN
         IF(.NOT.ASSOCIATED(A % ParMatrix)) CALL ParallelInitMatrix(Solver,A)

         ALLOCATE(B % ParallelInfo)
         ALLOCATE(B % Perm(2*n))
         ALLOCATE(B % ParallelInfo % GInterface(2*n))
         ALLOCATE(B % ParallelInfo % GlobalDOFs(2*n))
         ALLOCATE(B % ParallelInfo % NeighbourList(2*n))
 
         max_gdofs = ParallelReduction( MAXVAL(A % ParallelInfo % GlobalDOFs), 2 )
         DO i=1,n
           B % ParallelInfo % NeighbourList(i) % Neighbours => A % ParallelInfo % NeighbourList(i) % Neighbours
           B % ParallelInfo % NeighbourList(i+n) % Neighbours => A % ParallelInfo % NeighbourList(i) % Neighbours

           B % ParallelInfo % Ginterface(i) = A % ParallelInfo % GInterface(i)
           B % ParallelInfo % Ginterface(i+n) = A % ParallelInfo % GInterface(i)

           B % ParallelInfo % GlobalDOFs(i) = A % ParallelInfo % GlobalDOFs(i)
           B % ParallelInfo % GlobalDOFs(i+n) = A % ParallelInfo % GlobalDOFs(i)+max_gdofs

           B % Perm(i) = A % Perm(i)
           B % Perm(i+n) = A % Perm(i)+n
         END DO
         B % Parmatrix => ParInitMatrix(B, B % ParallelInfo)
       END IF
       n = 2*n
     END FUNCTION GenerateStateEquationSystem
    
!------------------------------------------------------------------------------
END SUBROUTINE SolveEigenSystem
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Compute lumped fluxes, for example for capacitance or impedance matrices. 
!------------------------------------------------------------------------------
SUBROUTINE StoreLumpedFluxes( Solver, NoModes, iMode, FluxesRow, FluxesRowIm, &
    FluxesRhs, FluxesRhsIm, ImpRe, ImpIm ) 
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  INTEGER :: NoModes, iMode
  REAL(KIND=dp) :: FluxesRow(:)
  REAL(KIND=dp), OPTIONAL :: FluxesRowIm(:)
  REAL(KIND=dp), OPTIONAL :: FluxesRhs, FluxesRhsIm
  REAL(KIND=dp), OPTIONAL :: ImpRe, ImpIm
!------------------------------------------------------------------------------
  REAL(KIND=dp), POINTER :: FluxesMatrix(:,:), FluxesMatrixIm(:,:)
  INTEGER :: i,j,k,n,Nmode,Mmode
  REAL(KIND=dp), POINTER :: PValues(:)
  TYPE(LumpedModel_t), POINTER :: Lumped
  CHARACTER(*), PARAMETER :: Caller = 'StoreLumpedFluxes'
      
  CALL Info(Caller,'Storing lumped fluxes',Level=10)

  Lumped => Solver % Lumped
  IF(.NOT. ASSOCIATED( Lumped ) ) THEN
    CALL Info(Caller,'Allocating lumped model of size: '//I2S(NoModes), Level=5)
    ALLOCATE( Solver % Lumped )
    Lumped => Solver % Lumped

    Lumped % NoModes = NoModes
    ALLOCATE( Lumped % CMatrix( NoModes, NoModes ) )
    Lumped % CMatrix = 0.0_dp

    IF( PRESENT(FluxesRowIm) ) THEN
      CALL Info(Caller,'Storing lumped fluxes imaginary',Level=20)
      Lumped % IsComplex = .TRUE.
      ALLOCATE( Lumped % CMatrixIm( NoModes, NoModes ) )
      Lumped % CMatrixIm = 0.0_dp
    END IF
    IF( PRESENT(FluxesRhs) ) THEN
      ALLOCATE( Lumped % CRhs( NoModes ) )
      Lumped % Crhs = 0.0_dp
    END IF
    IF( PRESENT(FluxesRhsIm) ) THEN
      Lumped % IsComplex = .TRUE.
      ALLOCATE( Lumped % CRhsIm( NoModes ) )
      Lumped % CrhsIm = 0.0_dp
    END IF
    IF( PRESENT(ImpRe) ) THEN
      ALLOCATE( Lumped % ImpRe( NoModes ) )
      Lumped % ImpRe = 0.0_dp
    END IF
    IF( PRESENT(ImpIm) ) THEN
      ALLOCATE( Lumped % ImpIm( NoModes ) )
      Lumped % ImpIm = 0.0_dp
    END IF
  END IF

  IF( Lumped % IsComplex ) FluxesMatrixIm => Solver % Lumped % CMatrixIm

  Lumped % CMatrix(iMode,1:NoModes) = FluxesRow(1:NoModes)
  IF(PRESENT(FluxesRowIm) ) THEN
    Lumped % CMatrixIm(iMode,1:NoModes) = FluxesRowIm(1:NoModes)
  END IF
  IF(PRESENT(FluxesRhs) ) THEN
    Lumped % CRhs(iMode) = FluxesRhs
  END IF
  IF(PRESENT(FluxesRhsIm) ) THEN
    Lumped % CRhsIm(iMode) = FluxesRhsIm
  END IF
  IF( PRESENT(ImpRe) ) THEN
    Lumped % ImpRe(iMode) = ImpRe
  END IF
  IF( PRESENT(ImpIm) ) THEN
    Lumped % ImpIm(iMode) = ImpIm
  END IF
        
  Lumped % CntModes = iMode

  BLOCK
    INTEGER :: ADepth
    IF( ListGetLogicalAnySolver( CurrentModel,'Adaptive Mesh Refinement') ) THEN      
      Adepth = CurrentModel % Mesh % AdaptiveDepth
      WRITE(Message,*) 'Row'//I2S(iMode)//' Depth'//I2S(Adepth)//':',FluxesRow(1:NoModes)
      CALL Info('StoreLumpedFluxes',Message,Level=4)
    END IF
  END BLOCK
  
END SUBROUTINE StoreLumpedFluxes


SUBROUTINE FinalizeLumpedMatrix( Solver )

  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
  REAL(KIND=dp), POINTER :: FluxesMatrix(:,:), FluxesMatrixIm(:,:)
  INTEGER :: i,j,k,NoModes
  LOGICAL :: Symmetric, IsComplex, EmWaveMode, CoilMode, Found, DoIt
  REAL(KIND=dp) :: nrm
  CHARACTER(:), ALLOCATABLE :: MatrixFile
  CHARACTER(*), PARAMETER :: Caller = 'FinalizeLumpedMatrix'
  COMPLEX(KIND=dp) :: cx, cb
  REAL(KIND=dp), ALLOCATABLE :: Checksum(:)
  TYPE(LumpedModel_t), POINTER :: Lumped

  CALL Info(Caller,'Finalizing lumped matrix',Level=8)
  
  IF(.NOT. ASSOCIATED( Solver % Lumped ) ) THEN
    CALL Fatal(Caller,'We should not be here without the lumped system!')
  END IF

  EmWaveMode = ListGetLogical( Solver % Values,'Constraint Modes EM Wave', Found )
  CoilMode = ListGetLogical( Solver % Values,'Constraint Modes Coils', Found )
  
  Lumped => Solver % Lumped  
  NoModes = Lumped % NoModes
  IF( Lumped % CntModes /= NoModes ) THEN
    CALL Fatal(Caller,'Trying to deduce '//I2S(NoMOdes)//&
        ' rows with '//I2S(Lumped % CntModes)//' lines of data!')
  END IF
  
  IsComplex = Lumped % IsComplex

  FluxesMatrix => Lumped % CMatrix
  IF( IsComplex ) FluxesMatrixIm => Lumped % CMatrixIm

  ! Energies and Impedances are related by factor two when currents are "1".   
  IF(CoilMode) THEN
    BLOCK
      REAL(KIND=dp) :: DesiredCurr(NoModes)  
      ! Normalize the inductance in case it has been computed with non-unity currents!
      DO i=1,CurrentModel % NumberOfComponents
        j = ListGetInteger( CurrentModel % Components(i) % Values,'Constraint Mode',Found )
        IF(j > 0) THEN
          DesiredCurr(j) = ListGetCReal( CurrentModel % Components(i) % Values,'Desired Coil Current',Found )
          IF(.NOT. Found) DesiredCurr(j) = 1.0_dp
        END IF
      END DO
      DO i=1,NoModes
        DO j=1,NoModes
          FluxesMatrix(i,j) = FluxesMatrix(i,j) / (DesiredCurr(i)*DesiredCurr(j)) 
        END DO
      END DO
    END BLOCK
  END IF

  ! Normalize the S-parameter matrix
  IF(EmWaveMode) THEN    
    ! Normalize by the source    
    BLOCK
      
      LOGICAL :: FixIt      
      FixIt =  ListGetLogical( Solver % Values,'Enforce Unity rowsum',Found )

      IF( InfoActive(20) ) THEN        
        CALL Info( Caller,'Showing matrix before normalization!')
        DO i=1,NoModes
          DO j=1,NoModes
            WRITE( Message, '(I3,I3,2ES17.9)' ) i,j,FluxesMatrix(i,j),FluxesMatrixIm(i,j)
            CALL Info( Caller, Message )
          END DO
        END DO
        DO i=1,NoModes
          WRITE( Message,*) 'Normalization vector '//I2S(i)//':',Lumped % Crhs(i)
          CALL Info( Caller, Message )
        END DO
      END IF      
   
      DO i=1,NoModes
        DO j=1,NoModes         
          nrm = SQRT(Lumped % Crhs(j) * Lumped % Crhs(i))                               
          FluxesMatrix(i,j) = FluxesMatrix(i,j) / nrm
          FluxesMatrixIm(i,j) = FluxesMatrixIm(i,j) / nrm
        END DO
      END DO
        
      IF( FixIt ) THEN
        DO i=1,NoModes
          nrm = 2*FluxesMatrix(i,i) / SUM(FluxesMatrix(i,:)**2 + FluxesMatrixIm(i,:)**2) 
          FluxesMatrix(i,:) = FluxesMatrix(i,:) * nrm
          FluxesMatrixIm(i,:) = FluxesMatrixIm(i,:) * nrm
          
          WRITE(Message,*) 'Fixing Multiplier '//I2S(i)//':',nrm
          CALL Info( Caller, Message, Level=6)          
        END DO
      END IF

      DO i=1,NoModes
        FluxesMatrix(i,i) = FluxesMatrix(i,i) - 1.0_dp
      END DO
    END BLOCK    
  END IF
  
  Symmetric = ListGetLogical( Solver % Values,&
      'Constraint Modes Fluxes Symmetric', Found ) 
  IF(.NOT. Found) Symmetric = ListGetLogical( Solver % Values,&
      'Constraint Modes Matrix Symmetric', Found ) 
  IF( Symmetric ) THEN        
    CALL Info(Caller,'Enforcing symmetry of reduced system!',Level=8)
    
    IF( InfoActive(10) ) THEN
      CALL Info( Caller,'Showing asymmetry of matrix before enforced symmetry!')
      DO i=1,NoModes
        DO j=i+1,NoModes
          IF( IsComplex ) THEN
            WRITE( Message, '(I3,I3,2ES17.9)' ) i,j,&
                FluxesMatrix(i,j)-FluxesMatrix(j,i),FluxesMatrixIm(i,j)-FluxesMatrixIm(j,i)
          ELSE                
            WRITE( Message, '(I3,I3,ES17.9)' ) i,j,FluxesMatrix(i,j)-FluxesMatrix(j,i)
          END IF
          CALL Info( Caller, Message )
        END DO
      END DO
    END IF
    
    FluxesMatrix = 0.5_dp * ( FluxesMatrix + TRANSPOSE( FluxesMatrix ) )
    IF( IsComplex ) THEN
      FluxesMatrixIm = 0.5_dp * ( FluxesMatrixIm + TRANSPOSE( FluxesMatrixIm ) )
    END IF
  END IF
  
  CALL Info( Caller,'Lumped Matrix', Level=5 )
  DO i=1,NoModes
    DO j=1,NoModes
      IF( Symmetric .AND. j < i ) CYCLE
      IF( IsComplex ) THEN
        WRITE( Message, '(I3,I3,2ES17.9)' ) i,j,FluxesMatrix(i,j),FluxesMatrixIm(i,j)
      ELSE
        WRITE( Message, '(I3,I3,ES17.9)' ) i,j,FluxesMatrix(i,j)
      END IF
      CALL Info( Caller, Message, Level=5 )
    END DO
  END DO

  MatrixFile = ListGetString(Solver % Values,'Constraint Modes Fluxes Filename',Found )
  IF(.NOT. Found) MatrixFile = ListGetString(Solver % Values,'Constraint Modes Matrix Filename',Found )
  IF( Found ) THEN
    ! Find the lowest active partition working here.
    k = ParallelReduction(ParEnv % MyPe,1)     
    IF( k == ParEnv % MyPe ) THEN
      OPEN(10, FILE=MatrixFile)
      DO i=1,NoModes
        WRITE (10,*) FluxesMatrix(i,:)
      END DO
      CLOSE(10)

      IF( IsComplex ) THEN
        OPEN( 11, FILE=TRIM(MatrixFile)//'_im')
        DO i=1,NoModes
          WRITE (11,*) FluxesMatrixIm(i,:) 
        END DO
        CLOSE(11)
        
        ALLOCATE(CheckSum(NoModes))
        OPEN( 11, FILE=TRIM(MatrixFile)//'_abs')
        DO i=1,NoModes
          WRITE (11,*) SQRT(FluxesMatrix(i,:)**2+FluxesMatrixIm(i,:)**2)           
          CheckSum(i) = SUM(FluxesMatrix(i,:)**2+FluxesMatrixIm(i,:)**2) 
        END DO
        CLOSE(11)

        OPEN( 11, FILE=TRIM(MatrixFile)//'_angle')
        DO i=1,NoModes
          WRITE (11,*) ( 180.0_dp / PI ) * ATAN2(FluxesMatrixIm(i,:),FluxesMatrix(i,:))           
        END DO
        CLOSE(11)

        IF( ASSOCIATED( Lumped % ImpRe ) ) THEN
          OPEN( 11, FILE=TRIM(MatrixFile)//'_Z')
          DO i=1,NoModes
            WRITE (11,*) Lumped % ImpRe(i), Lumped % ImpIm(i) 
          END DO
          CLOSE(11)
        END IF

        WRITE(Message,*) 'Normalization checksum: ',CheckSum
        CALL Info(Caller,Message) 
      END IF
      CALL Info( Caller,'Constraint modes fluxes was saved to file '//TRIM(MatrixFile),Level=5)
    END IF
  END IF
  
  nrm = 0.0_dp
  DO i=1,NoModes
    DO j=1,NoModes
      nrm = nrm + ABS(FluxesMatrix(i,j))
      IF(IsComplex) nrm = nrm + ABS(FluxesMatrixIm(i,j))
    END DO
  END DO

  WRITE(Message,'(A,ES12.5)') TRIM(Solver % Variable % Name)//' lumped matrix norm: ',Nrm
  CALL Info(Caller, Message, Level=6)

  DoIt = ListGetLogical( Solver % Values,'Constraint Modes Fluxes Norm', Found )
  IF(.NOT. Found) DoIt = ListGetLogical( Solver % Values,'Constraint Modes Matrix Norm', Found )
  IF(DoIt) CALL ListAddConstReal( CurrentModel % Simulation,&
      'res: '//TRIM(Solver % Variable % Name)//' lumped matrix norm',nrm)
  
  IF( ListGetLogical( Solver % Values,'Constraint Modes Matrix Results', Found ) ) THEN
    CALL Info(Caller,'Adding Constraint Modes Fluxes with "res:" to list',Level=5)
    DO i=1,NoModes
      DO j=1,NoModes
        CALL ListAddConstReal( CurrentModel % Simulation,'res: CMF '//I2S(10*i+j),FluxesMatrix(i,j))
      END DO
    END DO
    IF( IsComplex ) THEN
      DO i=1,NoModes
        DO j=1,NoModes
          CALL ListAddConstReal( CurrentModel % Simulation,'res: CMF Im '//I2S(10*i+j),FluxesMatrixIm(i,j))
        END DO
      END DO
    END IF
  END IF

END SUBROUTINE FinalizeLumpedMatrix


! Given an edge element and associated field variable computes the circulation over the
! element. The consistency of the circulation is checked 1st ensuring the each edge is
! accounted for correctly and the element is consistently accounted for with the reference
! normal.
!------------------------------------------------------------------------------------------
SUBROUTINE EdgeElementCirculation( Circ, Mesh, Element, Nodes, Solver, RefNormal, Var ) 
  REAL(KIND=dp) :: Circ(:)
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  TYPE(Solver_t), POINTER :: Solver
  TYPE(Variable_t), POINTER, OPTIONAL :: Var
  REAL(KIND=dp), OPTIONAL :: RefNormal(3)

  REAL(KIND=dp) :: EdgeVector(3), EdgeVector1(3), Normal(3), DotProd, s
  INTEGER :: dofs, dofIndexes(20), i, j, k, n, i1, i2, j1, j2, nb
  INTEGER :: ownParent, totParent
  TYPE(Variable_t), POINTER :: pVar

  Circ = 0.0_dp
  nb = mGetElementDOFs( dofIndexes, Element, Solver )
  n = Element % Type % NumberOfEdges
  
  IF(PRESENT(Var)) THEN
    pVar => Var
  ELSE
    pVar => Solver % Variable 
  END IF

  IF( ParEnv % PEs > 1 ) THEN
    TotParent = 0
    OwnParent = 0   
    IF(ASSOCIATED(Element % BoundaryInfo % Left) ) THEN
      TotParent = TotParent + 1
      IF(Element % Boundaryinfo % Left % PartIndex == ParEnv % MyPe ) ownParent = ownParent + 1
    END IF
    IF(ASSOCIATED(Element % BoundaryInfo % Right) ) THEN
      TotParent = TotParent + 1
      IF(Element % Boundaryinfo % Right % PartIndex == ParEnv % MyPe ) ownParent = ownParent + 1
    END IF
    IF(ownParent == 0) RETURN
  END IF
  
  DO i=1,nb
    i1 = i
    i2 = MODULO(i,n)+1

    ! Vector in the direction of the edge
    EdgeVector(1) = Nodes % x(i2) - Nodes % x(i1)
    EdgeVector(2) = Nodes % y(i2) - Nodes % y(i1)
    EdgeVector(3) = Nodes % z(i2) - Nodes % z(i1)

    ! Among different elements define a consistent direction for circulation
    ! Store 1st tangent vector and compute the normal from the cross product. 
    IF( PRESENT(RefNormal) ) THEN
      IF(i==1) THEN
        Edgevector1 = -EdgeVector
      ELSE IF(i==2) THEN
        Normal = CrossProduct( EdgeVector1, EdgeVector )
        DotProd = SUM(Normal * RefNormal) 
      END IF
    END IF

    ! Among edges for one element define consistent direction for circulation
    ! Edges are (1,2), (2,3), (3,1) for triangle and likewise for quad.
    j1 = Element % NodeIndexes(i1)
    j2 = Element % NodeIndexes(i2) 
    IF( ParEnv % PEs > 1 ) THEN                            
      j1 = Mesh % ParallelInfo % GlobalDOFs(i1)             
      j2 = Mesh % ParallelInfo % GlobalDOFs(i2)             
    END IF

    ! Integration length and direction
    s = SQRT(SUM(EdgeVector**2))
    IF( j1 < j2) s = -s

    j = pVar % Perm(DofIndexes(i))
    DO k=1,pVar % Dofs
      Circ(k) = Circ(k) + s * pVar % Values(pVar % Dofs*(j-1)+k)
    END DO
  END DO

  ! If the element normal is pointing to different direction reverse the flux.
  IF( PRESENT(RefNormal) ) THEN
    IF( DotProd < 0.0_dp ) Circ = -Circ
  END IF

  ! If we are at interface of partitions the face element may appear twice. Then take just half. 
  IF( ParEnv % PEs > 1 ) THEN
    Circ = Circ * ownParent / totParent    
  END IF
  
END SUBROUTINE EdgeElementCirculation


SUBROUTINE BoundaryCirculation( BCCirc, BcInd, Solver, Var ) 
  REAL(KIND=dp) :: BCCirc(:)
  INTEGER :: BcInd 
  TYPE(Solver_t), POINTER :: Solver
  TYPE(Variable_t), POINTER :: Var
  
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t), SAVE :: Nodes
  REAL(KIND=dp) :: Normal0(3), ElemCirc(3), parCoeff
  INTEGER :: i,t,k,n,nb,elem,dofs,ownParent,TotParent
  LOGICAL :: NormalSet
  TYPE(Solver_t), POINTER :: pSolver
  
  Mesh => Solver % Mesh
  IF( Mesh % MeshDim < 3 ) THEN
    CALL Fatal('BoundaryCirculation','Currently only available in 3D!')
  END IF
  pSolver => Solver
  dofs = Var % dofs
  BCCirc = 0.0_dp
  NormalSet = .FALSE.
  
  DO elem=Mesh % NumberOfBulkElements + 1, &
      Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements        
    
    Element => Mesh % Elements(elem)    

    IF(.NOT. ASSOCIATED(Element % BoundaryInfo) ) CYCLE
    IF(Element % BoundaryInfo % Constraint /= CurrentModel % BCs(bcInd) % Tag) CYCLE
    
    n = Element % Type % NumberOfNodes
    CALL CopyElementNodesFromMesh( Nodes, Solver % Mesh, n, Element % NodeIndexes)

    ! Compute a reference normal to which all other face normals are compared against.
    IF( .NOT. NormalSet ) THEN
      Normal0 = NormalVector(Element,Nodes)
      NormalSet = .TRUE.
    END IF
      
    CALL EdgeElementCirculation( ElemCirc, Mesh, Element, Nodes, pSolver, Normal0, Var )     
    BcCirc(1:dofs) = BcCirc(1:dofs) + ElemCirc(1:dofs)
  END DO

  IF( ParEnv % PEs > 1 ) THEN
    DO i=1,dofs
      BcCirc(i) = ParallelReduction(BCCirc(i))
    END DO
  END IF
    
END SUBROUTINE BoundaryCirculation



!------------------------------------------------------------------------------
!> Solve a linear system with permutated constraints.
!------------------------------------------------------------------------------
SUBROUTINE ConstraintModesDriver( A, x, b, Solver, PreSolve, ThisMode, LinSysModes, FirstLoop )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t), TARGET :: Solver
    REAL(KIND=dp) CONTIG :: x(:),b(:)
    LOGICAL :: PreSolve
    INTEGER, OPTIONAL :: ThisMode, LinSysModes
    LOGICAL, OPTIONAL :: FirstLoop
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i,j,k,n,NoModes,Mmode,ierr
    LOGICAL :: PrecRecompute, Stat, Found, ComputeFluxes, ComputeLinkage, Symmetric, &
        IsComplex, Parallel, ConsiderP, RhsMode, EmWaveMode, CoilMode, GaussLaw
    REAL(KIND=dp), ALLOCATABLE :: Fluxes(:), b0(:), A0(:), TempRhs(:)
    REAL(KIND=dp), ALLOCATABLE :: FluxesRow(:), FluxesRowIm(:)
    REAL(KIND=dp) :: FluxesRhs, FluxesRhsIm, ImpRe, ImpIm
    LOGICAL, ALLOCATABLE :: ConstrainedDOF0(:)
    REAL(KIND=dp) :: flux
    CHARACTER(:), ALLOCATABLE :: MatrixFile, BCName
    INTEGER :: NMode = 0, dof
    TYPE(Variable_t), POINTER :: pVar
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: LinsysMode, EigenMode, GotBC, LumpedMode 
    CHARACTER(*), PARAMETER :: Caller = 'ConstraintModesDriver'

    SAVE FluxesRow, FluxesRowIm, Fluxes, TempRhs, A0, b0, ConstrainedDOF0, LinsysMode, NMode

    
    !------------------------------------------------------------------------------
    IF(PRESENT(LinSysModes)) LinSysModes = 0
    NoModes = Solver % NumberOfConstraintModes 

    IF(NoModes <= 0) RETURN    
    Params => Solver % Values

    ! We can also have a combination of standard analysis + constraint modes
    ! analysis of the frozen state. Then the default solution slot should really do
    ! the standard analysis.
    IF( ListGetLogical( Params,'Constraint Modes Analysis Frozen',Found ) ) THEN
      Solver % Variable % FrozenMode = .TRUE.
      RETURN
    END IF

    EigenMode = ListgetLogical( Solver % Values, 'Eigen Analysis', Found ) 
    
    Var => Solver % Variable
    n = A % NumberOfRows        
    Parallel = Solver % Parallel
    
    IsComplex = ListGetLogical( Params,'Linear System Complex',Found)


    ! If the mode is nodal it is not lumped
    ! If it relates to whole boundary it is. 
    LumpedMode = ListGetLogical( Params,'Constraint Modes Lumped',Found )
    
    ! This is to my understanding not needed. To estimate the fluxes we
    ! basically integrate over basis functions that estimate unity.
    ! For p-elements this means using the linear nodal basis only, not any
    ! of the fake fluxes associated to p-degrees of freedom. Anyways,
    ! we leave this option here for testing etc. 
    ConsiderP = ListGetLogical( Params,'Consider P Fluxes',Found ) 

    ComputeFluxes = ListGetLogical( Params,'Constraint Modes Fluxes',Found) 
    ComputeLinkage = ListGetLogical( Params,'Constraint Modes Linkage',Found )
    RhsMode = ListGetLogical( Params,'Constraint Modes rhs',Found )
    EmWaveMode = ListGetLogical( Params,'Constraint Modes EM Wave', Found ) 
    CoilMode = ListGetLogical( Params,'Constraint Modes Coils',Found ) 

    IF(EmWaveMode) IsComplex = .TRUE.
    GaussLaw = .FALSE.
    IF(EMWaveMode) GaussLaw = ListGetLogical( Params,'Use Gauss Law',Found ) 
    
    IF( EmWaveMode .OR. CoilMode ) THEN
      ComputeFluxes = .TRUE.
      RhsMode = .TRUE.
    END IF

    
    IF( PreSolve ) THEN
      CALL Info(Caller,'Number of constraint modes is: '//I2S(NoModes),Level=8)
          
      ! We loop over the mode if it is not given in some external loop.
      !---------------------------------------------------------------------
      pVar => NULL()
      IF( ListGetLogical( Params,'Nonlinear System Constraint Modes', Found ) ) &
          pVar => VariableGet( Solver % Mesh % Variables,'nonlin iter')    
      IF( ListGetLogical( Params,'Steady State Constraint Modes', Found ) ) &
          pVar => VariableGet( Solver % Mesh % Variables,'coupled iter')
      IF( ListGetLogical( Params,'Run Control Constraint Modes', Found ) .OR. &
          ListGetLogical( CurrentModel % Control,'Constraint Modes Analysis', Found ) ) &
          pVar => VariableGet( Solver % Mesh % Variables,'run')    
      LinSysMode = .NOT. ASSOCIATED(pVar)
      
      IF(LinSysMode) THEN
        LinSysModes = NoModes
        ! If we combined eigen analysis base + constraint modes base do an empty
        ! cycle that does the eigenmodes first without updating the counter. 
        IF( PRESENT( FirstLoop ) ) THEN
          IF( FirstLoop .AND. EigenMode ) RETURN
        END IF        
        Nmode = ThisMode + 1
      ELSE
        Nmode = NINT( pVar % Values(1) ) 
        LinSysModes = 0
      END IF
      ThisMode = Nmode
      
      IF(EigenMode) CALL ListAddLogical( Solver % Values,'Skip Eigen Analysis',.TRUE.)

      IF( SIZE(x) /= n ) THEN
        CALL Fatal(Caller,'Conflicting sizes for matrix and variable! ('//I2S(SIZE(x))//','//I2S(n)//')')
      END IF

      IF( ComputeFluxes .OR. ComputeLinkage ) THEN
        ALLOCATE( FluxesRow(NoModes), Fluxes(n) )       
        IF( IsComplex ) ALLOCATE( FluxesRowIm(NoModes) )         
        
        IF( Parallel ) THEN
          ALLOCATE(TempRHS(SIZE(A % BulkRhs)))
          TempRhs = 0.0_dp
        END IF
        
        IF( IsComplex .OR. CoilMode) THEN
          ALLOCATE( A0(SIZE(A % Values)), b0(n), ConstrainedDOF0(n) )
          A0 = A % Values
          b0 = A % Rhs
          ConstrainedDOF0 = A % ConstrainedDOF
        END IF
      END IF
      
      IF(LinSysMode .AND. NMode == 2 ) THEN
        CALL ListAddLogical( Params,'No Precondition Recompute',.TRUE.)
      END IF

      CALL Info(Caller,'Setting up constrained mode: '//I2S(NMode),Level=6)
      i = Nmode


      ! By default constraint modes are set to 0/1.
      ! However, we can also set the BC's in some other way using prefix "mode 1:" etc.  
      GotBC = .FALSE.
      IF( LumpedMode ) THEN
        DO dof=1,Var % dofs
          BcName = 'mode '//I2S(Nmode)//': '//ComponentName(Var % name,dof)
          IF(ListCheckPresentAnyBC(CurrentModel, BcName ) ) GotBC = .TRUE.
        END DO
      END IF
        
      ! The matrix has been manipulated already before. This ensures
      ! that the system has values 1 at the constraint mode i.
      IF( CoilMode ) THEN                
        b = b0
        WHERE( Var % ConstraintModesIndeces > 0 .AND. Var % ConstraintModesIndeces /= Nmode ) 
          b = 0.0_dp
        END WHERE
      ELSE IF( RhsMode ) THEN
        IF( IsComplex ) THEN       
          A % Values = A0
          WHERE( Var % ConstraintModesIndeces > 0 )
            b = 0.0_dp
          END WHERE
          WHERE( Var % ConstraintModesIndeces == 2*Nmode-1 .OR. &
              Var % ConstraintModesIndeces == 2*Nmode )
            ! Revert rhs to previous
            b = b0
          END WHERE
        ELSE       
          IF( Nmode > 1 .AND. LinSysMode ) THEN
            WHERE( Var % ConstraintModesIndeces == Nmode-1 ) 
              b = 0.0_dp
            END WHERE
          END IF
          WHERE( Var % ConstraintModesIndeces == Nmode ) 
            b = Var % ConstraintModesWeights
          END WHERE
        END IF
      ELSE        
        IF( IsComplex ) THEN
          ! Quick and a little dirty fix for complex capacitance matrix.
          IF( ListGetLogical( Solver % Values,'Calculate Capacitance Matrix',Found ) ) THEN
            WHERE( Var % ConstraintModesIndeces /= 0 ) 
              A % ConstrainedDOF = .TRUE.
              A % DValues = 0.0_dp
            END WHERE
            WHERE( Var % ConstraintModesIndeces == 2*Nmode-1 )
              A % DValues = 1.0_dp
            END WHERE
          ELSE
            A % Values = A0
            b = b0
            A % ConstrainedDOF = ConstrainedDOF0

            WHERE( Var % ConstraintModesIndeces == 2*Nmode-1 )
              A % ConstrainedDOF = .TRUE.
              A % DValues = 1.0_dp
            END WHERE
            WHERE( Var % ConstraintModesIndeces == 2*Nmode ) 
              A % ConstrainedDOF = .TRUE.
              A % DValues = 0.0_dp
            END WHERE
          END IF
          CALL EnforceDirichletConditions( Solver, A, b )

        ELSE IF( GotBC ) THEN
          
          IF( Nmode > 1 .AND. LinSysMode ) THEN
            DO dof=1,Var % dofs
              WHERE( Var % ConstraintModesIndeces == Var % Dofs*(Nmode-2)+dof ) 
                A % DValues = 0.0_dp
              END WHERE
            END DO
          END IF
          
          DO dof=1,Var % dofs
            BcName = 'mode '//I2S(Nmode)//': '//ComponentName(Var % name,dof)
            IF(ListCheckPresentAnyBC(CurrentModel, BcName ) ) THEN            
              CALL Info(Caller,"Setting constraint for: "//TRIM(BCName),Level=7)
              CALL SetDirichletBoundaries( CurrentModel, A, b, &
                  BcName, dof, Var % DOFs, Var % Perm )
            END IF
          END DO

          CALL EnforceDirichletConditions( Solver, A, b )
          
        ELSE
          
          IF( Nmode > 1 .AND. LinSysMode ) THEN
            WHERE( Var % ConstraintModesIndeces == Nmode-1 ) 
              A % DValues = 0.0_dp
            END WHERE
          END IF

          WHERE( Var % ConstraintModesIndeces == Nmode ) 
            A % DValues = 1.0_dp
          END WHERE

          CALL EnforceDirichletConditions( Solver, A, b )                    
        END IF
      END IF
      CALL ListAddLogical( Params,'Skip Zero Rhs Test',.TRUE. )
    END IF
      
    IF( .NOT. PreSolve ) THEN 
      IF( PRESENT( FirstLoop ) ) THEN
        IF( FirstLoop .AND. EigenMode ) RETURN
      END IF

      CALL Info(Caller,'Mode '//I2S(NMode)//' computed, doing some postprocessing',Level=10)

      IF( .NOT. ( IsComplex .OR. CoilMode ) ) THEN
        WHERE( Var % ConstraintModesIndeces == Nmode ) b = 0.0_dp
      END IF
      
      IF( NMode <= Var % NumberOfConstraintModes ) THEN
        Var % ConstraintModes(NMode,:) = x
      END IF

      IF( ComputeFluxes .OR. ComputeLinkage ) THEN
        CALL Info(Caller,'Computing lumped fluxes',Level=8)
        
        IF( CoilMode ) THEN
          CALL MagneticEnergies()
        ELSE IF(EmWaveMode ) THEN
          CALL EMWaveFluxes()
        ELSE IF( ComputeFluxes ) THEN
          CALL ConstraintModesFluxes(EmWaveMode)
        ELSE IF( ComputeLinkage ) THEN
          CALL ConstraintModesLinkage()
        END IF
        
        ! Do parallel communication here at one sweep, not before!
        IF(.NOT. EMWaveMode ) THEN
          CALL CommunicateConstraintModesFluxes()
        END IF
        
        IF( IsComplex ) THEN
          IF(EmWaveMode ) THEN
            CALL StoreLumpedFluxes(Solver, NoModes, NMode, FluxesRow, FluxesRowIm, FluxesRhs, FluxesRhsIm, ImpRe, ImpIm)
          ELSE
            CALL StoreLumpedFluxes(Solver, NoModes, NMode, FluxesRow, FluxesRowIm, FluxesRhs, FluxesRhsIm)
          END IF
        ELSE
          CALL StoreLumpedFluxes(Solver, NoModes, NMode, FluxesRow ) 
        END IF
      END IF
            
      IF(LinSysMode .AND. NMode == NoModes ) THEN
        IF( ComputeFluxes .OR. ComputeLinkage .OR. CoilMode ) THEN
          CALL FinalizeLumpedMatrix( Solver )            
        END IF
        CALL ListAddLogical( Params,'No Precondition Recompute',.FALSE.)
      END IF

      IF( ComputeFluxes .OR. ComputeLinkage ) THEN
        DEALLOCATE( Fluxes, FluxesRow )
        IF( IsComplex ) DEALLOCATE( FluxesRowIm )
        IF( Parallel ) DEALLOCATE(TempRHS)
        IF( IsComplex .OR. CoilMode) THEN
          A % Values = A0
          A % Rhs = b0
          A % ConstrainedDOF = ConstrainedDOF0           
          DEALLOCATE( A0, b0, ConstrainedDOF0 )
        END IF
      END IF

      CALL ListAddLogical( Params,'Skip Zero Rhs Test',.FALSE. )

      IF(EigenMode) THEN
        CALL ListAddLogical( Solver % Values,'Skip Eigen Analysis',.FALSE.)
      END IF


    END IF


    
  CONTAINS

    SUBROUTINE MagneticEnergies()
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Variable_t), POINTER :: CVar, AVar
      TYPE(ValueList_t), POINTER :: Vlist
      INTEGER, POINTER :: MasterEntities(:)
      REAL(KIND=dp) :: Nrm, FL
      
      Mesh => Solver % Mesh
      CVar => VariableGet( Mesh % Variables,'CoilCurrent e',ThisOnly=.TRUE.)
      IF(.NOT. ASSOCIATED(CVar)) THEN
        CVar => VariableGet( Mesh % Variables,'CoilCurrent',ThisOnly=.TRUE.)
      END IF
      AVar => Solver % Variable

      FluxesRow = 0.0_dp      
      DO j=1,CurrentModel % NumberOfComponents
        Vlist => CurrentModel % Components(j) % Values       
        k = ListGetInteger( Vlist,'Constraint Mode', Found )
        IF(.NOT. Found) CYCLE

        IF( ListGetLogical( Vlist,'Flux Linkage', Found ) ) THEN
          FL = ComponentStokesTheorem(CurrentModel, Mesh, Vlist, AVar, Surf = .TRUE.)
          WRITE(Message,'(A,ES12.3)') 'Flux Linkage '//I2S(NMode)//' '//I2S(k)//':', FL
          CALL Info(Caller,Message)
        END IF

        MasterEntities => ListGetIntegerArray( Vlist,'Master Bodies',Found )
        FluxesRow(k) = ComponentCoilEnergy(CurrentModel, Mesh, MasterEntities, avar, cvar )

        WRITE(Message,'(A,ES12.3)') 'Coil Inductance '//I2S(Nmode)//' '//I2S(k)//':', FluxesRow(k)
        CALL Info(Caller,Message)
      END DO
            
      Nrm = SUM(FluxesRow)
      WRITE(Message,'(A,ES12.3)') 'Energy norm of Impedance matrix row '//I2S(NMode)//': ',Nrm
      CALL Info(Caller,Message,Level=8)

    END SUBROUTINE MagneticEnergies


    SUBROUTINE EMWaveFluxes()
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Variable_t), POINTER :: AVar
      TYPE(ValueList_t), POINTER :: Vlist
      INTEGER, POINTER :: MasterEntities(:)
      COMPLEX(KIND=dp) :: OutFlux,InFlux,PortFlux,InImp,PortImp
      INTEGER :: i,j,k,n,port,alloc     
      
      CALL Info(Caller,'Using <Ej,Ej> for lumping',Level=10)
      
      Mesh => Solver % Mesh
      AVar => Solver % Variable      
      FluxesRow = 0.0_dp
      FluxesRowIm = 0.0_dp

      FluxesRhs = 1.0_dp 
      FluxesRhsIm = 0.0_dp

      DO port = 1, NoModes
        DO alloc = 0, 1
          n = 0
          DO j=1,CurrentModel % NumberOfBCs        
            Vlist => CurrentModel % BCs(j) % Values       
            k = ListGetInteger( Vlist,'Constraint Mode', Found )
            IF(k == port) THEN
              n = n+1
              IF(alloc == 1 ) MasterEntities(n) = j
            END IF
          END DO
          IF(alloc == 1 ) EXIT
          ALLOCATE(MasterEntities(n))
        END DO

        OutFlux = BoundaryWaveFlux(CurrentModel, Mesh, MasterEntities, Avar, PortFlux, PortImp, port==NMode )
          
        ! Memorize the coefficient for normalization: <Ec,Ej>/<Ei,Ei>                
        ! Real and imag part of: <Ec,Ej>
        FluxesRow(port) = REAL(OutFlux)
        FluxesRowIm(port) = AIMAG(OutFlux) 

        ! Memorize diagonal entry <Ej,Ej*> for future normalization
        IF(port==NMode) THEN
          InFlux = PortFlux
          InImp = PortImp
        END IF
               
        DEALLOCATE(MasterEntities)                
      END DO

      FluxesRhs = REAL(InFlux)      
      FluxesRhsIm = AIMAG(InFlux)
      ImpRe = REAL(InImp)
      ImpIm = AIMAG(InImp)
      
    END SUBROUTINE EMWaveFluxes

      
    
    SUBROUTINE ConstraintModesFluxes(EmWaveMode)
      LOGICAL :: EmWaveMode

      REAL(KIND=dp), POINTER CONTIG :: PValues(:)
      INTEGER :: poffset
      REAL(KIND=dp) :: flux, w
      COMPLEX(KIND=dp) :: cflux, cx, cmult, cb, crhs
      
      ! Use the initial bulk values that do not include Dirichlet conditions
      PValues => A % Values
      IF( .NOT. ASSOCIATED( A % BulkValues ) ) THEN
        CALL Fatal(Caller,'BulkValues not associated!')
      END IF
      A % Values => A % BulkValues
            
      Fluxes = 0.0_dp
      IF ( Parallel ) THEN
        TempRHS = A % BulkRhs 
        CALL ParallelInitSolve( A, x, TempRHS, Fluxes )
        CALL ParallelMatrixVector( A, x, Fluxes, .TRUE. )
      ELSE
        CALL MatrixVectorMultiply( A, x, Fluxes )
      END IF

      ! Revert pointer back 
      A % Values => PValues

      poffset = 2*(NoModes + 1)
      FluxesRow = 0.0_dp
      IF(IsComplex) THEN
        FluxesRowIm = 0.0_dp
        FluxesRhs = 0.0_dp
        FluxesRhsIm = 0.0_dp
      END IF
      
      IF( EmWaveMode ) THEN
        w = ListGetAngularFrequency( Found = Found )
        IF(.NOT. Found) CALL Fatal(Caller,'Energy mode requires "Angular Frequency"!')
        cmult = 1.0/(2*w*CMPLX(0.0_dp,1.0_dp)) 
      END IF
      
      DO j=1,n
        k = Var % ConstraintModesIndeces(j)
        
        IF( ConsiderP ) THEN
          ! P dofs are associated with negative index as they are not included in ConstraintModesAnalysis.
          IF( k < -1 ) k = k + poffset
        END IF

        IF( k > 0 ) THEN          
          IF( IsComplex ) THEN
            Mmode = (k+1)/2
            IF( MOD(k,2) == 1 ) THEN                
              IF( EmWaveMode ) THEN
                cx = CMPLX(x(j),x(j+1))
                cflux = cmult * cx * CONJG(CMPLX(Fluxes(j),Fluxes(j+1)))                
                crhs = cmult * cx * CONJG(CMPLX(b(j),b(j+1)))
              ELSE
                cflux = CMPLX(Fluxes(j),Fluxes(j+1))                            
                crhs = CMPLX(b(j),b(j+1))
              END IF
              IF( Nmode /= Mmode ) THEN
                FluxesRow(Mmode) = FluxesRow(Mmode) - REAL(cflux)
                FluxesRowIm(Mmode) = FluxesRowIm(Mmode) - AIMAG(cflux)
              ELSE
                FluxesRow(Nmode) = FluxesRow(Nmode) + REAL(cflux)
                FluxesRowIm(Nmode) = FluxesRowIm(Nmode) + AIMAG(cflux)
                FluxesRhs = FluxesRhs + REAL(crhs)
                FluxesRhsIm = FluxesRhsIm + AIMAG(crhs)
              END IF
            END IF
          ELSE
            flux = Fluxes(j)
            Mmode = k 
            IF( Nmode /= Mmode ) THEN                
              FluxesRow(Mmode) = FluxesRow(Mmode) - flux
            END IF
            FluxesRow(Nmode) = FluxesRow(Nmode) + flux
          END IF
        END IF
      END DO

      IF( EmWaveMode ) THEN
        IF(InfoActive(20)) THEN
          PRINT *,'ConstrainModesFluxes:',FluxesRow
          PRINT *,'ConstrainModesFluxesIm:',FluxesRowIm
          PRINT *,'ConstrainModesFluxes:',FluxesRhs,FluxesRhsIm        
        END IF
      END IF
        
    END SUBROUTINE ConstraintModesFluxes

   
    SUBROUTINE ConstraintModesCirculation()
      TYPE(Element_t), POINTER :: Element
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Solver_t), POINTER :: pSolver
      INTEGER, POINTER :: dofIndexes(:), nodeIndexes(:), Perm(:)
      REAL(KIND=dp), POINTER :: Values(:)
      REAL(KIND=dp) :: Normal(3), Normal0(3), Circ(3)
      TYPE(Nodes_t), SAVE :: Nodes
      INTEGER :: i,t,k,n,nb,elem
      LOGICAL :: NormalSet
            
      Mesh => Solver % Mesh
      IF( Mesh % MeshDim < 3 ) THEN
        CALL Fatal('ConstraintModesCirculation','Currently only available in 2D!')
      END IF
      pSolver => Solver

      FluxesRow = 0.0_dp
      IF( IsComplex ) FluxesRowIm = 0.0_dp
      
      DO elem=Mesh % NumberOfBulkElements + 1, &
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements        

        Element => Mesh % Elements(elem)
        n = Element % Type % NumberOfNodes
        nodeIndexes => Element % NodeIndexes

        CALL CopyElementNodesFromMesh( Nodes, Solver % Mesh, n, Element % NodeIndexes)
                      
        nb = mGetElementDOFs( dofIndexes, Element, Solver )

        k = Solver % Variable % ConstraintModesIndeces(Perm(dofIndexes(1)))
        IF(k==0) CYCLE
              
        IF(ANY(Var % ConstraintModesIndeces(Perm(dofIndexes(1:nb))) /= k)) CYCLE        
        
        ! Compute a reference normal to which all other face normals are compared against.
        Normal0 = NormalVector(Element,Nodes)

        CALL EdgeElementCirculation( Circ, Mesh, Element, Nodes, pSolver, Normal0, Var ) 
         
        ! Store elemental fluxes to the total sums. 
        IF( IsComplex ) THEN            
          Mmode = (k+1)/2
          IF( MOD(k,2) == 1 ) THEN                
            FluxesRow(Mmode) = FluxesRow(Mmode) + Circ(1)
          ELSE
            FluxesRowIm(Mmode) = FluxesRowIm(Mmode) + Circ(2)
          END IF
        ELSE
          Mmode = k 
          FluxesRow(Mmode) = FluxesRow(Mmode) + Circ(1)
        END IF
      END DO

    END SUBROUTINE ConstraintModesCirculation

    
    SUBROUTINE ConstraintModesLinkage()
      TYPE(Element_t), POINTER :: Element
      TYPE(Mesh_t), POINTER :: Mesh
      INTEGER, POINTER :: Indexes(:), Perm(:)
      TYPE(GaussIntegrationPoints_t) :: IP
      TYPE(Nodes_t), SAVE :: Nodes
      REAL(KIND=dp) :: detJ, flux, s, PotAtIp, Basis(12), pot(12)
      INTEGER :: t,k,n,elem

      
      Mesh => Solver % Mesh
      IF( Mesh % MeshDim == 3 ) THEN
        CALL Fatal('ConstraintModesLinkage','Currently only available in 2D!')
      END IF

      FluxesRow = 0.0_dp
      IF( IsComplex ) FluxesRowIm = 0.0_dp
      Perm => Solver % Variable % Perm
      
      DO elem=Mesh % NumberOfBulkElements + 1, &
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements        

        Element => Mesh % Elements(elem)
        n = Element % Type % NumberOfNodes
        Indexes => Element % NodeIndexes        
        k = Var % ConstraintModesIndeces(Perm(Indexes(1)))

        IF(k==0) CYCLE
        
        IF(ANY(Var % ConstraintModesIndeces(Perm(Indexes(1:n))) /= k)) CYCLE        
        pot(1:n) = Solver % Variable % Values(Solver % Variable % Perm(Indexes(1:n)))

        CALL CopyElementNodesFromMesh( Nodes, Solver % Mesh, n, Indexes)
        IP = GaussPoints( Element )
        
        DO t=1,IP % n
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )          
          s = IP % s(t) * DetJ
          PotAtIp = SUM(Basis(1:n)*pot(1:n))

          flux = s * PotAtIp 
          
          IF( IsComplex ) THEN
            Mmode = (k+1)/2
            IF( MOD(k,2) == 1 ) THEN                
              FluxesRow(Mmode) = FluxesRow(Mmode) + flux
            ELSE
              FluxesRowIm(Mmode) = FluxesRowIm(Mmode) + flux
            END IF
          ELSE
            Mmode = k 
            FluxesRow(Mmode) = FluxesRow(Mmode) + flux
          END IF
        END DO        
      END DO

    END SUBROUTINE ConstraintModesLinkage
      

    SUBROUTINE CommunicateConstraintModesFluxes()

      REAL(KIND=dp), ALLOCATABLE :: tmpFluxesRow(:)

      IF( Parallel ) THEN
        ALLOCATE(tmpFluxesRow(NoModes))
        tmpFluxesRow = FluxesRow
        CALL MPI_ALLREDUCE(tmpFluxesRow, FluxesRow, NoModes, MPI_DOUBLE_PRECISION, &
            MPI_SUM, ELMER_COMM_WORLD, ierr)
        IF( IsComplex ) THEN
          tmpFluxesRow = FluxesRowIm
          CALL MPI_ALLREDUCE(tmpFluxesRow, FluxesRowIm, NoModes, MPI_DOUBLE_PRECISION, &
              MPI_SUM, ELMER_COMM_WORLD, ierr)        
          CALL MPI_ALLREDUCE(FluxesRhs, FluxesRhs, 1, MPI_DOUBLE_PRECISION, &
              MPI_SUM, ELMER_COMM_WORLD, ierr)        
          CALL MPI_ALLREDUCE(FluxesRhsIm, FluxesRhsIm, 1, MPI_DOUBLE_PRECISION, &
              MPI_SUM, ELMER_COMM_WORLD, ierr)        
        END IF
        DEALLOCATE(tmpFluxesRow)
      END IF

    END SUBROUTINE CommunicateConstraintModesFluxes
    
    
!------------------------------------------------------------------------------
  END SUBROUTINE ConstraintModesDriver
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> A parser of the variable name that returns the true variablename
!> where the inline options have been interpreted.
!------------------------------------------------------------------------------
  SUBROUTINE VariableNameParser(var_name, NoOutput, Global, Dofs, &
      IpVariable, ElemVariable, DgVariable, NodalVariable )

    CHARACTER(LEN=*)  :: var_name
    LOGICAL, OPTIONAL :: NoOutput, Global
    INTEGER, OPTIONAL :: Dofs
    LOGICAL, OPTIONAL :: IpVariable
    LOGICAL, OPTIONAL :: ElemVariable
    LOGICAL, OPTIONAL :: DgVariable
    LOGICAL, OPTIONAL :: NodalVariable

    INTEGER :: i,j,k,m

    IF(PRESENT(NoOutput)) NoOutput = .FALSE.
    IF(PRESENT(Global)) Global = .FALSE.
    IF(PRESENT(Dofs)) Dofs = 0
    IF(PRESENT(IpVariable)) IpVariable = .FALSE.
    IF(PRESENT(DgVariable)) DgVariable = .FALSE.      
    IF(PRESENT(ElemVariable)) ElemVariable = .FALSE.      
    IF(PRESENT(NodalVariable)) NodalVariable = .FALSE.      
    
    
    DO WHILE( var_name(1:1) == '-' )

      m = 0
      IF ( SEQL(var_name, '-nooutput ') ) THEN
        IF(PRESENT(NoOutput)) NoOutput = .TRUE.
        m = 10

      ELSE IF ( SEQL(var_name, '-global ') ) THEN
        IF(PRESENT(Global)) Global = .TRUE.
        m = 8

      ELSE IF ( SEQL(var_name, '-ip ') ) THEN
        IF(PRESENT(IpVariable)) IpVariable = .TRUE.      
        m = 4

      ELSE IF ( SEQL(var_name, '-dg ') ) THEN
        IF(PRESENT(DgVariable)) DgVariable = .TRUE.      
        m = 4

      ELSE IF ( SEQL(var_name, '-elem ') ) THEN
        IF(PRESENT(ElemVariable)) ElemVariable = .TRUE.      
        m = 6

      ELSE IF ( SEQL(var_name, '-nodal ') ) THEN
        IF(PRESENT(NodalVariable)) NodalVariable = .TRUE.      
        m = 7
      END IF

      IF( m > 0 ) THEN
        var_name(1:LEN(var_name)-m) = var_name(m+1:)
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


  !> Create permutation for fields on integration points, optionally with mask.
  !> The non-masked version is saved to Solver structure for reuse while the
  !> masked version may be unique to every variable. 
  !-----------------------------------------------------------------------------------
  SUBROUTINE CreateIpPerm( Solver, MaskPerm, MaskName, SecName, UpdateOnly )

    TYPE(Solver_t), POINTER :: Solver
    INTEGER, POINTER, OPTIONAL :: MaskPerm(:)
    LOGICAL, OPTIONAL :: UpdateOnly
    CHARACTER(*), OPTIONAL :: MaskName, SecName      

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Element_t), POINTER :: Element
    INTEGER :: t, n, IpCount , RelOrder, nIp
    LOGICAL :: Found, ActiveElem, ActiveElem2
    INTEGER, POINTER :: IpOffset(:) 
    TYPE(ValueList_t), POINTER :: BF
    LOGICAL :: UpdatePerm
    CHARACTER(:), ALLOCATABLE :: EquationName

    n = 0
    IF( PRESENT( MaskPerm ) ) n = n + 1
    IF( PRESENT( MaskName ) ) n = n + 1
    IF( PRESENT( SecName ) ) n = n + 1
    IF( PRESENT( UpdateOnly ) ) n = n + 1

    ! Currently a lazy check
    IF( n /= 0 .AND. n /= 3 .AND. n /= 2) THEN
      CALL Fatal('CreateIpPerm','Only some optional parameter combinations are possible')
    END IF

    UpdatePerm = .FALSE.
    IF( PRESENT( UpdateOnly ) ) UpdatePerm = UpdateOnly

    IF( UpdatePerm ) THEN
      CALL Info('CreateIpPerm','Updating IP permutation table',Level=8)       
    ELSE IF( PRESENT( MaskPerm ) ) THEN
      CALL Info('CreateIpPerm','Creating masked permutation for integration points',Level=8)
    ELSE       
      IF( ASSOCIATED( Solver % IpTable ) ) THEN
        CALL Info('CreateIpPerm','IpTable already allocated, returning',Level=8)
      END IF
      CALL Info('CreateIpPerm','Creating permutation for integration points',Level=8)
    END IF

    EquationName = ListGetString( Solver % Values, 'Equation', Found)
    IF( .NOT. Found ) THEN
      CALL Fatal('CreateIpPerm','Equation not present!')
    END IF

    Mesh => Solver % Mesh
    NULLIFY( IpOffset ) 

    n = Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements

    IF( UpdatePerm ) THEN
      IpOffset => MaskPerm
      ActiveElem = (IpOffset(2)-IpOffset(1) > 0 )
      IF( n >= 2 ) ActiveElem2 = (IpOffset(3)-IpOffset(2) > 0 )
    ELSE
      ALLOCATE( IpOffset( n + 1) )     
      IpOffset = 0
      IF( PRESENT( MaskPerm ) ) MaskPerm => IpOffset
    END IF
    IpCount = 0

    nIp = ListGetInteger( Solver % Values,'Gauss Points on Ip Variables', Found ) 

    DO t=1,Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements
      Element => Mesh % Elements(t)

      IF( .NOT. UpdatePerm ) THEN
        ActiveElem = .FALSE.
        IF( Element % PartIndex == ParEnv % myPE ) THEN
          IF ( CheckElementEquation( CurrentModel, Element, EquationName ) ) THEN             
            IF( PRESENT( MaskName ) ) THEN
              BF => ListGetSection( Element, SecName )
              ActiveElem = ListGetLogicalGen( BF, MaskName )
            ELSE
              ActiveElem = .TRUE.
            END IF
          END IF
        END IF
      END IF

      IF( ActiveElem ) THEN
        IF( nIp > 0 ) THEN
          IpCount = IpCount + nIp
        ELSE
          IP = GaussPointsAdapt( Element )
          IpCount = IpCount + Ip % n
        END IF
      END IF

      ! We are reusing the permutation table hence we must be one step ahead 
      IF( UpdatePerm .AND. n >= t+1) THEN
        ActiveElem = ActiveElem2
        ActiveElem2 = (IpOffset(t+2)-IpOffset(t+1) > 0 )
      END IF

      IpOffset(t+1) = IpCount
    END DO

    IF( .NOT. PRESENT( MaskPerm ) ) THEN
      ALLOCATE( Solver % IpTable ) 
      Solver % IpTable % IpOffset => IpOffset
      Solver % IpTable % IpCount = IpCount
    END IF

    IF( UpdatePerm ) THEN
      CALL Info('CreateIpPerm','Updated permutation for IP points: '//I2S(IpCount),Level=8)  
    ELSE       
      CALL Info('CreateIpPerm','Created permutation for IP points: '//I2S(IpCount),Level=8)  
    END IF

  END SUBROUTINE CreateIpPerm


  SUBROUTINE UpdateIpPerm( Solver, Perm )

    TYPE(Solver_t), POINTER :: Solver
    INTEGER, POINTER :: Perm(:)

    CALL CreateIpPerm( Solver, Perm, UpdateOnly = .TRUE.)

  END SUBROUTINE UpdateIpPerm



!------------------------------------------------------------------------------
!> Updates values for exported variables which are typically auxiliary variables derived
!> from the solution.
!------------------------------------------------------------------------------
  SUBROUTINE UpdateExportedVariables( Solver )  
!------------------------------------------------------------------------------
    TYPE(Solver_t), TARGET :: Solver

    INTEGER :: i,j,k,l,n,m,t,bf_id,dofs,nsize,i1,i2,NoGauss
    REAL(KIND=dp), POINTER :: Values(:), Solution(:), LocalSol(:), LocalCond(:)
    INTEGER, POINTER :: Indexes(:), VarIndexes(:), Perm(:)
    LOGICAL :: Found, Conditional, GotIt, Stat, StateVariable, DoIt, AllocationsDone = .FALSE.
    LOGICAL, POINTER :: ActivePart(:),ActiveCond(:),ActivePartBC(:),ActiveCondBC(:)
    TYPE(Variable_t), POINTER :: ExpVariable
    TYPE(ValueList_t), POINTER :: ValueList
    TYPE(Element_t),POINTER :: Element  
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    REAL(KIND=dp) :: detJ
    TYPE(ValueHandle_t) :: LocalSol_h
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: pSolver
    CHARACTER(LEN=MAX_NAME_LEN) :: var_name
    CHARACTER(:), ALLOCATABLE :: str, tmpname,condname
    CHARACTER(*), PARAMETER :: Caller = 'UpdateExportedVariables'


    SAVE LocalSol_h

    CALL Info(Caller,'Updating variables, if any!',Level=20)

    AllocationsDone = .FALSE.
    Mesh => Solver % Mesh

    l = 0
    DO WHILE( .TRUE. )
      l = l + 1

      str = ComponentName( 'exported variable', l )    

      var_name = ListGetString( Solver % Values, str, GotIt )    
      IF(.NOT. GotIt) EXIT

      CALL Info(Caller,'Trying to set values for variable: '//TRIM(Var_name),Level=20)

      CALL VariableNameParser( var_name ) 

      ExpVariable => VariableGet( Mesh % Variables, Var_name )
      IF( .NOT. ASSOCIATED(ExpVariable)) CYCLE

      CALL Info(Caller,'Setting values for variable: '//TRIM(Var_name),Level=20)

      IF(.NOT. AllocationsDone ) THEN
        m = CurrentModel % NumberOFBodyForces
        ALLOCATE( ActivePart(m), ActiveCond(m) )

        m = CurrentModel % NumberOFBCs
        ALLOCATE( ActivePartBC(m), ActiveCondBC(m) )

        m = Mesh % MaxElementDOFs
        ALLOCATE( LocalSol(m), LocalCond(m))

        m =  CurrentModel % MaxElementNodes
        ALLOCATE( Basis(m), Nodes % x(m), Nodes % y(m), Nodes % z(m) )

        AllocationsDone = .TRUE.
      END IF

      Dofs = ExpVariable % DOFs
      Values => ExpVariable % Values
      Perm => ExpVariable % Perm
      n = LEN_TRIM( var_name )

      StateVariable = ( SIZE( Values ) == DOFs ) .OR. ( ExpVariable % Type == Variable_Global ) 
      IF( StateVariable ) THEN
        CALL Info(Caller,'Updating state variable',Level=20)
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
            CALL Info(Caller,&
                'Found a proper definition for state variable',Level=6)
            Solution = ListGetCReal( CurrentModel % BodyForces(bf_id) % Values,TmpName)
            EXIT
          END IF
        END DO
        CYCLE
      END IF

      CALL Info(Caller,'Updating field variable with dofs: '//I2S(DOFs),Level=12)

      
      DO j=1,DOFs

100     Values => ExpVariable % Values
        IF( Dofs > 1 ) THEN
          tmpname = ComponentName( var_name(1:n), j )
          Solution => Values( j:: DOFs ) 
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
        ActivePartBC = .FALSE.
        ActiveCondBC = .FALSE.

        DO bf_id=1,CurrentModel % NumberOFBodyForces
          ActivePart(bf_id) = ListCheckPresent( &
              CurrentModel % BodyForces(bf_id) % Values,TmpName ) 
          ActiveCond(bf_id) = ListCheckPresent( &
              CurrentModel % BodyForces(bf_id) % Values,CondName )      
        END DO
        DO bf_id=1,CurrentModel % NumberOFBCs
          ActivePartBC(bf_id) = ListCheckPresent( &
              CurrentModel % BCs(bf_id) % Values,TmpName ) 
          ActiveCondBC(bf_id) = ListCheckPresent( &
              CurrentModel % BCs(bf_id) % Values,CondName )      
        END DO

        m = COUNT(ActivePart) + COUNT(ActivePartBC)
        IF (m == 0) CYCLE
        CALL Info(Caller,'Exported Variable '//I2S(l)//' defined in '//I2S(m)//' sections',Level=8)

        IF( ExpVariable % TYPE == Variable_on_gauss_points ) THEN 
          ! Initialize handle when doing values on Gauss points!
          CALL ListInitElementKeyword( LocalSol_h,'Body Force',TmpName )
        END IF

        DO t = 1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

          Element => Mesh % Elements(t) 

          DoIt = .FALSE.
          IF( Element % BodyId > 0 ) THEN
            bf_id = ListGetInteger( CurrentModel % Bodies(Element % BodyId) % Values,'Body Force',GotIt)
            IF( bf_id > 0 ) THEN
              ValueList => CurrentModel % BodyForces(bf_id) % Values          
              DoIt = ActivePart(bf_id)
              IF(DoIt) Conditional = ActiveCond(bf_id)
            END IF
          END IF
          IF( .NOT. DoIt .AND. t > Mesh % NumberOfBulkElements ) THEN
            ! If we don't have an active "body force" section check still the boundary section.
            IF(ASSOCIATED( Element % BoundaryInfo ) ) THEN
              DO bf_id=1,CurrentModel % NumberOfBCs
                IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bf_id) % Tag ) EXIT
              END DO
              IF ( bf_id <= CurrentModel % NumberOfBCs ) THEN            
                ValueList => CurrentModel % BCs(bf_id) % Values                         
                DoIt = ActivePartBC(bf_id)
                IF(DoIt) Conditional = ActiveCondBC(bf_id)
              END IF
            END IF
          END IF

          IF(.NOT. DoIt) CYCLE
          
          CurrentModel % CurrentElement => Element
          m = Element % TYPE % NumberOfNodes
          Indexes => Element % NodeIndexes

          IF( ExpVariable % TYPE == Variable_on_gauss_points ) THEN 

            i1 = Perm( Element % ElementIndex )
            i2 = Perm( Element % ElementIndex + 1 )
            NoGauss = i2 - i1

            ! This is not active here
            IF( NoGauss == 0 ) CYCLE

            IP = GaussPointsAdapt( Element, Solver )

            IF( NoGauss /= IP % n ) THEN

              CALL Info(Caller,&
                  'Number of Gauss points has changed, redoing permutations!',Level=8)

              pSolver => Solver
              CALL UpdateIpPerm( pSolver, Perm )
              nsize = MAXVAL( Perm )

              CALL Info(Caller,'Total number of new IP dofs: '//I2S(nsize),Level=7)

              IF( SIZE( ExpVariable % Values ) /= ExpVariable % Dofs * nsize ) THEN
                DEALLOCATE( ExpVariable % Values )
                ALLOCATE( ExpVariable % Values( nsize * ExpVariable % Dofs ) )
              END IF
              ExpVariable % Values = 0.0_dp
              GOTO 100 
            END IF

            Nodes % x(1:m) = Mesh % Nodes % x(Indexes)
            Nodes % y(1:m) = Mesh % Nodes % y(Indexes)
            Nodes % z(1:m) = Mesh % Nodes % z(Indexes)

            IF( Conditional ) THEN
              CALL Warn(Caller,'Elemental variable cannot be conditional!')
            END IF

            DO k=1,IP % n
              stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
                  IP % W(k), detJ, Basis )
              Solution(i1+k) = ListGetElementReal( LocalSol_h,Basis,Element,Found,GaussPoint=k) 
            END DO

          ELSE IF( ExpVariable % TYPE == Variable_on_elements ) THEN
            IF( Conditional ) THEN
              CALL Warn(Caller,'Elemental variables not conditional!')
            END IF
            LocalSol(1:m) = ListGetReal(ValueList, TmpName, m, Indexes(1:m) )
            i = Perm( Element % ElementIndex ) 
            IF( i > 0 ) Solution(i) = SUM( LocalSol(1:m) ) / m

          ELSE
            IF( ExpVariable % TYPE == Variable_on_nodes_on_elements ) THEN
              VarIndexes => Element % DGIndexes
            ELSE
              VarIndexes => Indexes
            END IF

            LocalSol(1:m) = ListGetReal(ValueList, TmpName, m, Indexes(1:m) )

            IF( Conditional ) THEN
              LocalCond(1:m) = ListGetReal(ValueList, CondName, m, Indexes(1:m) )
              DO i=1,m
                IF( LocalCond(i) > 0.0_dp ) THEN
                  IF( Perm(VarIndexes(i)) > 0 ) THEN
                    Solution( Perm(VarIndexes(i)) ) = LocalSol(i)
                  END IF
                END IF
              END DO
            ELSE
              IF( ALL( Perm(VarIndexes(1:m)) > 0 ) ) THEN
                Solution( Perm(VarIndexes(1:m)) ) = LocalSol(1:m)
              END IF
            END IF

          END IF
        END DO

      END DO
    END DO


    l = 0
    DO WHILE( .TRUE. )
      l = l + 1
      str = ComponentName( 'project variable', l )    

      var_name = ListGetString( Solver % Values, str, GotIt )    
      IF(.NOT. GotIt) EXIT
      CALL Info(Caller,'Variable requested for projection: '//TRIM(var_name),Level=20)
      
      CALL VariableNameParser( var_name )     
      ExpVariable => VariableGet( Mesh % Variables, Var_name )
      IF( .NOT. ASSOCIATED(ExpVariable)) THEN
        CALL Warn(Caller,'Could not find variable for projection: '//TRIM(Var_name))
        CYCLE
      END IF

      IF( ExpVariable % TYPE /= Variable_on_gauss_points ) THEN
        CALL Fatal(Caller,'Variable projection implemented only for IP variable!')
      END IF

      k = Variable_on_nodes

      str = ComponentName( 'project type', l )
      tmpname = ListGetString( Solver % Values, str, GotIt )    

      k = Variable_on_nodes
      IF( GotIt ) THEN
        IF( tmpname == 'elem' ) THEN
          k = Variable_on_elements
        ELSE IF( tmpname == 'dg' ) THEN
          k = Variable_on_nodes_on_elements
        END IF
      END IF

      CALL Info(Caller,'Variable type for projection: '//I2S(k),Level=20)

      CALL Ip2DgSwapper( Mesh, ExpVariable, ToType = k )

      CALL Info(Caller,'Finished projection of variable',Level=30)
    END DO

    IF( AllocationsDone ) THEN
      DEALLOCATE(ActivePart, ActiveCond, ActivePartBC, ActiveCondBC, &
          LocalSol, LocalCond, Basis, Nodes % x, Nodes % y, Nodes % z )
    END IF
      
  END SUBROUTINE UpdateExportedVariables


!------------------------------------------------------------------------------
!> Derivates values for exported variables to come up with velocity and
!> acceleration fields.
!------------------------------------------------------------------------------
  SUBROUTINE DerivateExportedVariables( Solver )  
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var, DerVar, dtVar
  INTEGER :: VarNo, Cnt
  LOGICAL :: Found, DoIt
  REAL(KIND=dp) :: dt
  CHARACTER(LEN=MAX_NAME_LEN) :: str, var_name
  
  
  CALL Info('DerivateExportedVariables','Derivating variables, if any!',Level=20)

  Mesh => Solver % Mesh
  Params => Solver % Values

  VarNo = 0
  Cnt = 0
  
  DO WHILE( .TRUE. )
    VarNo = VarNo + 1

    str = ComponentName( 'exported variable', VarNo )    
    
    var_name = ListGetString( Solver % Values, str, Found )    
    IF(.NOT. Found) EXIT
    
    CALL VariableNameParser( var_name ) 

    Var => VariableGet( Mesh % Variables, Var_name )
    IF( .NOT. ASSOCIATED(Var)) CYCLE
    IF( .NOT. ASSOCIATED(Var % PrevValues) ) CYCLE
    
    str = TRIM( ComponentName(Var_name) )//' Calculate Velocity'
    DoIt = ListGetLogical( Params, str, Found )        
    IF( DoIt ) THEN
      str = TRIM( ComponentName(var_name) ) // ' Velocity'
      DerVar => VariableGet( Solver % Mesh % Variables, str )        
      IF(.NOT. ASSOCIATED(DerVar)) THEN
        CALL Warn('DerivateExportedVariables','Variable does not exist:'//TRIM(str))
        CYCLE
      END IF

      dtVar => VariableGet( Solver % Mesh % Variables, 'timestep size' )
      dt = dtVar % Values(1) 
      
      CALL Info('DerivateExportedVariables','Computing numerical derivative for:'//TRIM(str),Level=8)     
      DerVar % Values = (Var % Values(:) - Var % PrevValues(:,1)) / dt
      Cnt = Cnt + 1
    END IF

    str = TRIM( ComponentName(Var_name) )//' Calculate Acceleration'
    DoIt = ListGetLogical( Params, str, Found )        
    IF( DoIt ) THEN
      str = TRIM( ComponentName(var_name) ) // ' Acceleration'
      DerVar => VariableGet( Solver % Mesh % Variables, str )        
      IF(.NOT. ASSOCIATED(DerVar)) THEN
        CALL Warn('DerivateExportedVariables','Variable does not exist:'//TRIM(str))
        CYCLE
      END IF

      dtVar => VariableGet( Solver % Mesh % Variables, 'timestep size' )
      dt = dtVar % Values(1) 

      CALL Info('DerivateExportedVariables','Computing numerical derivative for:'//TRIM(str),Level=8)     
      DerVar % Values = (Var % Values(:) - 2*Var % PrevValues(:,1) - Var % PrevValues(:,2)) / dt**2
      Cnt = Cnt + 1
    END IF

  END DO

  CALL Info('DerivateExportedVariables','Derivating done for variables: '//I2S(Cnt),Level=20)

END SUBROUTINE DerivateExportedVariables



!------------------------------------------------------------------------------
!> Solves a harmonic system.
!------------------------------------------------------------------------------
SUBROUTINE SolveHarmonicSystem( G, Solver )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), TARGET :: G
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: BMatrix, A => NULL()
    INTEGER :: i,j,k,n, kr, ki, DOFs, ne, niter
    LOGICAL :: stat, Found, OptimizeBW, Real_given,Imag_given
    REAL(KIND=dp) :: Omega, norm, s
    REAL(KIND=dp), POINTER :: Freqv(:,:)
    REAL(KIND=dp), ALLOCATABLE :: x(:)
    REAL(KIND=dp), POINTER :: b(:)
    REAL(KIND=dp) :: frequency
    INTEGER :: Nfrequency
    TYPE(ValueList_t), POINTER :: BC
    CHARACTER(:), ALLOCATABLE :: Name

    CALL Info( 'SolveHarmonicSystem', 'Solving initially transient style system as harmonic one', Level=5)
    
    n = Solver % Matrix % NumberofRows
    DOFs = Solver % Variable % DOFs * 2

    A => G
    DO WHILE( ASSOCIATED(A) )
      BMatrix => A
      A => A % EMatrix
      IF ( ASSOCIATED(A) ) THEN
        IF ( A % COMPLEX ) THEN
          CALL Info('SolveHarmonicSystem','Reusing existing harmonic system',Level=10)
          EXIT
        END IF
      END IF
    END DO

    IF ( .NOT. ASSOCIATED(A) ) THEN      
      CALL Info('SolveHarmonicSystem','Creating new matrix for harmonic system',Level=10)      

      OptimizeBW = ListGetLogical(Solver % Values, 'Optimize Bandwidth', Found)
      IF ( .NOT. Found ) OptimizeBW = .TRUE.
      
      A => CreateMatrix( CurrentModel, Solver, Solver % Mesh,   &
              Solver % Variable % Perm, DOFs, MATRIX_CRS, OptimizeBW, &
              ListGetString( Solver % Values, 'Equation') )
      A % COMPLEX = .TRUE.
      BMatrix % EMatrix => A
      ALLOCATE( A % rhs(2*n) )
      
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
    END IF

    b => A % rhs
    ALLOCATE( x(2*n) )
    x = 0
    
    b(1:2*n:2) = G % RHS(1:n)
    b(2:2*n:2) = G % RHS_im(1:n)

    
    Nfrequency = ListGetInteger( Solver % Values,'Harmonic System Values',Found )
    IF( Nfrequency > 1 ) THEN
      freqv => ListGetConstRealArray( Solver % Values, 'Frequency' )
    ELSE
      Frequency = ListGetAngularFrequency( Solver % Values, Found ) / (2*PI)
      IF( .NOT. Found ) THEN
        CALL Fatal( 'SolveHarmonicSystem', '> Frequency < must be given for harmonic analysis.' )
      END IF
      
      Nfrequency = 1
      ! Add the number of frequencies even for case of one for some postprocessing stuff to work 
      CALL ListAddInteger( Solver % Values,'Harmonic System Values',Nfrequency )
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
      CALL Info( 'SolveHarmonicSystem', Message, Level=4 )

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

      CALL EnforceDirichletConditions( Solver, A, b )
 
      
      CALL SolveLinearSystem( A, b, x, Norm, DOFs, Solver )
      
      DO j=1,n
        Solver % Variable % EigenVectors(i,j) = &
                 CMPLX( x(2*(j-1)+1),x(2*(j-1)+2),KIND=dp )
      END DO
    END DO

    Solver % NOFEigenValues = ne

    DEALLOCATE( x )
!------------------------------------------------------------------------------
 END SUBROUTINE SolveHarmonicSystem
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Just toggles the initial system to harmonic one and back
!------------------------------------------------------------------------------
SUBROUTINE MergeSlaveSolvers( Solver, PreSolve )
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  LOGICAL :: PreSolve
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Solver_t), POINTER :: Solver2
  TYPE(Matrix_t), POINTER :: A1, A2, A
  TYPE(Variable_t), POINTER :: betaVar
  REAL(KIND=dp), POINTER :: vals1(:), vals2(:), vals(:)
  REAL(KIND=dp) :: alpha, invAlpha
  INTEGER :: i,j,k,n,dofs1,dofs2
  INTEGER, POINTER :: SlaveSolverIndexes(:)
  INTEGER, POINTER :: perm1(:), perm2(:), perm(:)
  CHARACTER(:), ALLOCATABLE :: str
  LOGICAL :: Found

  SAVE :: A1, perm1, dofs1, vals1, &
      A2, perm2, dofs2, vals2, A, perm, vals, &
      alpha, invAlpha, betaVar

  CALL Info('MergeSlaveSolvers','Monolithic treatment of solvers')
  
  IF( .NOT. ASSOCIATED( Solver % Variable ) ) THEN
    CALL Fatal('MergeSlaveSolvers','Not applicable without a variable')
    RETURN    
  END IF
  IF( .NOT. ASSOCIATED( Solver % Matrix ) ) THEN
    CALL Fatal('MergeSlaveSolvers','Not applicable without a matrix')
    RETURN    
  END IF
  
  Params => Solver % Values  
  SlaveSolverIndexes => ListGetIntegerArray( Params,'Slave Solvers',Found )
  IF(.NOT. Found ) RETURN
  
  IF(SIZE(SlaveSolverIndexes) > 1 ) THEN
    CALL Warn('MergeSlaveSolvers','Cannot current only deal with one slave solver!')   
  END IF
  i = SlaveSolverIndexes(1)

  IF(PreSolve) THEN
    A1 => Solver % Matrix
    perm1 => Solver % Variable % Perm
    dofs1 = Solver % Variable % Dofs
    vals1 => Solver % Variable % Values

    Solver2 => CurrentModel % Solvers(i)
    A2 => Solver2 % Matrix
    perm2 => Solver2 % Variable % Perm
    dofs2 = Solver2 % Variable % Dofs
    vals2 => Solver2 % Variable % Values

    IF(dofs1 /= 1 .OR. dofs2 /= 1 ) THEN
      CALL Fatal('MergeSlaveSolver','Enabled for 1 dof only: '//I2S(dofs1)//','//I2S(dofs2))
    END IF

    alpha = ListGetCReal( Params,'Slave Field Scaling',Found)
    IF(.NOT. Found) alpha = 1.0_dp
    invAlpha = 1.0_dp / alpha

    betaVar => NULL()
    str = ListGetString( Params,'Slave Field Offset Variable',Found )
    IF( Found ) THEN
      betaVar => VariableGet(Solver % Mesh % Variables,str)
      IF(.NOT. ASSOCIATED(betaVar)) THEN
        CALL Fatal('MergeSlaveSolver','Could not find offset variable: '//TRIM(str))
      END IF
      IF(ANY(betaVar % Perm /= perm2 ) ) THEN
        CALL Fatal('MergeSlaveSolver','The offset field should have same permutation!')
      END IF      
    END IF
      
    A => AllocateMatrix()
    CALL CRS_MergeMatrix(A1, A2, C = A, PermA = perm1, PermB = perm2, PermC = perm) 

    ALLOCATE( A % Diag(A % NumberOfRows) )
    A % diag = 0
    CALL CRS_SortMatrix( A, .TRUE. )
    
    Solver % Matrix => A
    Solver % Variable % Perm => Perm

    ALLOCATE( vals( A % NumberOfRows) )
    Solver % Variable % Values => vals

    CALL MergeRhsAndSolutions()
  ELSE
    CALL UnmergeSolutions()

    
    Solver % Matrix => A1
    Solver % Variable % Values => vals1
    Solver % Variable % Perm => perm1

    DEALLOCATE(vals)
    DEALLOCATE(perm)
    CALL FreeMatrix(A)
  END IF

CONTAINS

  SUBROUTINE MergeRhsAndSolutions()

    INTEGER :: i,j,j1,j2
    REAL(KIND=dp), ALLOCATABLE :: rhsadd(:)
    REAL(KIND=dp) :: c1,x2
       
    CALL Info('MergeSlaveSolvers','Merging rhs and initial guess for monolithic solution!',Level=10)
    
    IF(.NOT. ASSOCIATED(A % rhs)) THEN
      ALLOCATE(A % rhs(A % NumberOfRows))
      A % RHS = 0._dp
    END IF

    IF(ASSOCIATED(betaVar)) THEN
      ALLOCATE(rhsadd(SIZE(betaVar % Values)))
      rhsadd = 0.0_dp
      CALL MatrixVectorMultiply( A2, betaVar % Values, rhsadd )
      rhsAdd = -InvAlpha * rhsAdd 
    END IF
      
    DO i=1,SIZE(perm)
      j = perm(i)
      IF(j==0) CYCLE

      j1 = perm1(i)
      j2 = perm2(i)

      IF(j1>0) A % rhs(j) = A1 % rhs(j1) 
      IF(j2>0) THEN
        A % rhs(j) = A % rhs(j) + InvAlpha * A2 % rhs(j2) 
        IF(ASSOCIATED(betaVar)) A % rhs(j) = A % rhs(j) + rhsAdd(j2)
      END IF

      vals(j) = 0.0_dp
      c1 = 0.0_dp      
      IF(j1>0) THEN
        c1 = 1.0_dp
        vals(j) = vals1(j1) 
      END IF      
      IF(j2>0) THEN
        x2 = vals2(j2)
        IF(ASSOCIATED(betaVar)) x2 = x2 - betaVar % Values(j2)
        vals(j) = (vals(j) + InvAlpha * x2)/(c1+1.0_dp) 
      END IF
    END DO    
    
  END SUBROUTINE MergeRhsAndSolutions

  
  SUBROUTINE UnmergeSolutions()

    INTEGER :: i,j,j1,j2

    CALL Info('MergeSlaveSolvers','Unmerging the solution back to composite solvers!',Level=10)
    
    DO i=1,SIZE(perm)
      j = perm(i)
      IF(j==0) CYCLE
      j1 = perm1(i)
      j2 = perm2(i)
      IF(j1>0) vals1(j1) = vals(j)
      IF(j2>0) THEN
        vals2(j2) = Alpha * vals(j)
      END IF
    END DO

    ! The permutation has been checked for these
    IF(ASSOCIATED(BetaVar)) vals2 = vals2 + betaVar % Values

    
  END SUBROUTINE UnmergeSolutions
    
END SUBROUTINE MergeSlaveSolvers

 

!------------------------------------------------------------------------------
!> Just toggles the initial system to harmonic one and back
!------------------------------------------------------------------------------
SUBROUTINE ChangeToHarmonicSystem( Solver, BackToReal )
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  LOGICAL, OPTIONAL :: BackToReal
  !------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER :: Are => NULL(), Aharm => NULL(), SaveMatrix 
  INTEGER :: i,j,k,n, kr, ki, DOFs
  LOGICAL :: stat, Found, OptimizeBW, Real_given, Imag_given
  CHARACTER(:), ALLOCATABLE :: Name
  REAL(KIND=dp) :: Omega, s, val
  REAL(KIND=dp), POINTER :: b(:), TmpVals(:)
  REAL(KIND=dp) :: frequency
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Variable_t), POINTER :: TmpVar, ReVar, HarmVar, SaveVar
  LOGICAL :: ToReal, ParseName, AnyDirichlet, Diagonal, HarmonicReal, EigenMode
  
  IF( .NOT. ASSOCIATED( Solver % Variable ) ) THEN
    CALL Warn('ChangeToHarmonicSystem','Not applicable without a variable')
    RETURN    
  END IF

  IF( .NOT. ASSOCIATED( Solver % Matrix ) ) THEN
    CALL Warn('ChangeToHarmonicSystem','Not applicable without a matrix')
    RETURN    
  END IF

  EigenMode = ListgetLogical( Solver % Values, 'Eigen Analysis', Found )

  ToReal = .FALSE.
  IF( PRESENT( BackToReal ) ) ToReal = BackToReal

  IF( ToReal ) THEN
    IF( ASSOCIATED( Solver % Variable % Evar ) ) THEN
      IF( Solver % Variable % Evar % Dofs < Solver % Variable % Dofs ) THEN
        CALL Info('ChangeToHarmonicSystem','Changing the harmonic results back to real system!',Level=6)

        SaveVar => Solver % Variable
        SaveMatrix => Solver % Matrix 

        Solver % Variable => Solver % Variable % Evar
        Solver % Variable % Evar => SaveVar

        Solver % Matrix => Solver % Matrix % EMatrix
        Solver % Matrix % Ematrix => SaveMatrix

        ! Eliminate cyclic dependence that is a bummer when deallocating stuff
        NULLIFY( Solver % Matrix % EMatrix % Ematrix )
      END IF
    END IF
    RETURN
  END IF


  CALL Info('ChangeToHarmonicSystem','Changing the real transient system to harmonic one!',Level=6)

  SaveMatrix => Solver % Matrix
  SaveVar => Solver % Variable     

  n = Solver % Matrix % NumberofRows
  DOFs = SaveVar % Dofs
  Are => Solver % Matrix

  CALL Info('ChangeToHarmonicSystem','Number of real system rows: '//I2S(n),Level=16)

  ! Obtain the frequency, it may depend on iteration step etc. 
  Omega = 0._dp
  IF (.NOT. EigenMode) THEN
    Frequency = ListGetAngularFrequency( Solver % Values, Found ) / (2*PI)
    IF( .NOT. Found ) THEN
      CALL Fatal( 'ChangeToHarmonicSystem', '> Frequency < must be given for harmonic analysis.' )
    END IF
    WRITE( Message, '(a,e12.3)' ) 'Frequency value: ', frequency
    CALL Info( 'ChangeToHarmonicSystem', Message, Level=5 )

     omega = 2 * PI * Frequency
     CALL ListAddConstReal( CurrentModel % Simulation, 'res: frequency', Frequency )
  END IF

  
  HarmonicReal = ListGetLogical( Solver % Values,'Harmonic Mode Real',Found ) 
  IF( HarmonicReal ) THEN
    CALL Info('ChangeToHarmonicSystem','Enforcing harmonic system to be real valued',Level=8)
    IF (ASSOCIATED(Are % MassValues)) THEN
      ARe % Values = Are % Values - omega**2* Are % MassValues
    ELSE
      CALL Fatal('ChangeToHarmonicSystem','Harmonic system requires mass!')
    END IF
    ! This is set outside so that it can be called more flexibilly
    CALL EnforceDirichletConditions( Solver, Are, Are % rhs  )
    RETURN
  END IF

 
  Diagonal = ListGetLogical( Solver % Values,'Harmonic Mode Block Diagonal',Found )  
  IF(.NOT. Found ) Diagonal = .NOT. ASSOCIATED(Are % DampValues)
  IF( Diagonal ) THEN
    CALL Info('ChangeToHarmonicSystem','Undamped system is assumed to be block diagonal',Level=8)
  END IF

  
  ! Find whether the matrix already exists
  Aharm => Are % EMatrix
  IF( ASSOCIATED( Aharm ) ) THEN
    CALL Info('ChangeToHarmonicSystem','Found existing harmonic system',Level=10)
    IF( ALLOCATED( Aharm % ConstrainedDOF ) ) Aharm % ConstrainedDOF = .FALSE.
  ELSE    
    ! Create the matrix if it does not
    
    Aharm => CreateChildMatrix( Are, Dofs, 2*Dofs, CreateRhs = .TRUE., Diagonal = Diagonal )

    Aharm % COMPLEX = ListGetLogical( Solver % Values,'Linear System Complex', Found ) 
    IF( .NOT. Found ) Aharm % COMPLEX = .NOT. Diagonal !TRUE. 
  END IF


  ! Set the harmonic system r.h.s
  b => Aharm % rhs
  
  IF( ASSOCIATED( Are % Rhs ) ) THEN
    b(1:2*n:2) = Are % RHS(1:n)
  ELSE
    b(1:2*n:2) = 0.0_dp
  END IF
  
  IF( ASSOCIATED( Are % Rhs_im ) ) THEN
    b(2:2*n:2) = Are % RHS_im(1:n)            
  ELSE
    b(2:2*n:2) = 0.0_dp
  END IF

  IF( ASSOCIATED(Are % MassValues) ) THEN
    CALL Info('ChangeToHarmonicSystem','We have mass matrix values',Level=12)
  ELSE
    CALL Warn('ChangeToHarmonicSystem','We do not have mass matrix values!')
  END IF

  IF( ASSOCIATED(Are % DampValues) ) THEN
    CALL Info('ChangeToHarmonicSystem','We have damp matrix values',Level=12)
    IF( Diagonal ) THEN
      CALL Fatal('ChangeToHarmonicSystem','Damping matrix cannot be block diagonal!')
    END IF
  ELSE
    CALL Info('ChangeToHarmonicSystem','We do not have damp matrix values',Level=12)
  END IF

  ! Set the harmonic system matrix
  IF( EigenMode ) THEN
    ALLOCATE(Aharm % MassValues(SIZE(Aharm % Values)))
    Aharm % MassValues = 0.0_dp

    DO k=1,n
      kr = Aharm % Rows(2*(k-1)+1)
      ki = Aharm % Rows(2*(k-1)+2)
      DO j=Are % Rows(k),Are % Rows(k+1)-1
        Aharm % Values(kr) = Are % Values(j)
        Aharm % Values(ki+1) = Are % Values(j)
        
        IF (ASSOCIATED(Are % DampValues)) THEN
          Aharm % Values(kr+1) = -Are % Dampvalues(j)
          Aharm % Values(ki)   =  Are % Dampvalues(j)
        END IF

        IF (ASSOCIATED(Are % MassValues)) THEN
          Aharm % MassValues(kr) = Are % MassValues(j)
          Aharm % MassValues(ki+1) = Are % MassValues(j)
        END IF

        kr = kr + 2
        ki = ki + 2
      END DO
    END DO
  ELSE IF( Diagonal ) THEN
    DO k=1,n
      kr = Aharm % Rows(2*(k-1)+1)
      ki = Aharm % Rows(2*(k-1)+2)
      DO j=Are % Rows(k),Are % Rows(k+1)-1
        val = Are % Values(j)
        IF (ASSOCIATED(Are % MassValues)) val = val - omega**2* Are % MassValues(j)
        
        Aharm % Values(kr) = val 
        Aharm % Values(ki) = val 
        kr = kr + 1
        ki = ki + 1
      END DO
    END DO
  ELSE
    DO k=1,n
      kr = Aharm % Rows(2*(k-1)+1)
      ki = Aharm % Rows(2*(k-1)+2)
      DO j=Are % Rows(k),Are % Rows(k+1)-1
        val = Are % Values(j)
        IF (ASSOCIATED(Are % MassValues)) val = val - omega**2* Are % MassValues(j)

        Aharm % Values(kr) = val
        Aharm % Values(ki+1) = val     
        
        IF (ASSOCIATED(Are % DampValues)) THEN
          Aharm % Values(kr+1) = -Are % Dampvalues(j) * omega
          Aharm % Values(ki)   =  Are % Dampvalues(j) * omega
        END IF

        kr = kr + 2
        ki = ki + 2
      END DO
    END DO
  END IF
    
  AnyDirichlet = .FALSE.
  
  ! Finally set the Dirichlet conditions for the solver    
  DO j=1,DOFs
    Name = ComponentName( Solver % Variable % Name, j ) 
    DO i=1,CurrentModel % NumberOFBCs
      BC => CurrentModel % BCs(i) % Values
      real_given = ListCheckPresent( BC, Name )
      imag_given = ListCheckPresent( BC, TRIM(Name) // ' im' )

      IF( real_given .OR. imag_given ) AnyDirichlet = .TRUE.

      IF ( real_given .AND. .NOT. imag_given ) THEN
        CALL Info('ChangeToHarmonicSystem','Setting zero >'//TRIM(Name)//' im< on BC '//I2S(i),Level=12)
        CALL ListAddConstReal( BC, TRIM(Name) // ' im', 0._dp)
      ELSE IF ( imag_given .AND. .NOT. real_given ) THEN
        CALL Info('ChangeToHarmonicSystem','Setting zero >'//TRIM(Name)//'< on BC '//I2S(i),Level=12)
        CALL ListAddConstReal( BC, Name, 0._dp )
      END IF
    END DO
  END DO


  IF( AnyDirichlet ) THEN
    DO j=1,DOFs
      Name = ComponentName( SaveVar % Name, j ) 
      
      CALL SetDirichletBoundaries( CurrentModel, Aharm, b, Name, &
          2*j-1, 2*DOFs, SaveVar % Perm )

      CALL SetDirichletBoundaries( CurrentModel, Aharm, b, TRIM(Name) // ' im', &
          2*j, 2*DOFs, SaveVar % Perm )
    END DO
  END IF

  
  ! Create the new fields, the total one and the imaginary one
  !-------------------------------------------------------------
  k = INDEX( SaveVar % name, '[' )
  ParseName = ( k > 0 ) 

  ! Name of the full complex variable not used for postprocessing
  IF( ParseName ) THEN
    Name = TRIM(SaveVar % Name(1:k-1))//' complex'
  ELSE
    Name = TRIM( SaveVar % Name )//' complex'
  END IF

  CALL Info('ChangeToHarmonicSystem','Harmonic system full name: '//TRIM(Name),Level=12)


  HarmVar => VariableGet( Solver % Mesh % Variables, Name )
  IF( ASSOCIATED( HarmVar ) ) THEN
    CALL Info('ChangeToHarmonicSystem','Reusing full system harmonic dofs',Level=12)
  ELSE
    CALL Info('ChangeToHarmonicSystem','Creating full system harmonic dofs',Level=12)
    CALL VariableAddVector( Solver % Mesh % Variables,Solver % Mesh,Solver, &
        Name,2*DOFs,Perm=SaveVar % Perm,Output=.FALSE.)
    HarmVar => VariableGet( Solver % Mesh % Variables, Name )
    IF(.NOT. ASSOCIATED( HarmVar ) ) CALL Fatal('ChangeToHarmonicSystem','New created variable should exist!')

    ! Repoint the values of the original solution vector
    HarmVar % Values(1:2*n:2) = SaveVar % Values(1:n)

    ! It beats me why this cannot be deallocated without some NaNs later
    !DEALLOCATE( SaveVar % Values )
    SaveVar % Values => HarmVar % Values(1:2*n:2)
    SaveVar % Secondary = .TRUE.

    ! Repoint the components of the original solution
    IF( Dofs > 1 ) THEN
      DO i=1,Dofs
        TmpVar => VariableGet( Solver % Mesh % Variables, ComponentName( SaveVar % Name, i ) )
        IF( ASSOCIATED( TmpVar ) ) THEN
          TmpVar % Values => HarmVar % Values(2*i-1::HarmVar % Dofs)
        ELSE
          CALL Fatal('ChangeToHarmonicSystem','Could not find re component '//I2S(i))
        END IF
      END DO
    END IF

    IF( ParseName ) THEN
      Name = ListGetString( Solver % Values,'Imaginary Variable',Found )
      IF(.NOT. Found ) THEN
        CALL Fatal('ChangeToHarmonicSystem','We need > Imaginary Variable < to create harmonic system!')
      END IF
    ELSE
      Name = TRIM( SaveVar % Name )//' im'
      CALL Info('ChangeToHarmonicSystem','Using derived name for imaginary component: '//TRIM(Name),Level=12)
    END IF

    TmpVals => HarmVar % Values(2:2*n:2)
    CALL VariableAdd( Solver % Mesh % Variables,Solver % Mesh,Solver, &
        Name, DOFs,TmpVals, Perm=SaveVar % Perm,Output=.TRUE.,Secondary=.TRUE.)        

    IF( Dofs > 1 ) THEN
      DO i=1,Dofs
        TmpVals => HarmVar % Values(2*i:2*n:2*Dofs)
        CALL VariableAdd( Solver % Mesh % Variables,Solver % Mesh,Solver, &
            ComponentName(Name,i),1,TmpVals,Perm=SaveVar % Perm,Output=.TRUE.,Secondary=.TRUE.)        
      END DO
    END IF
    
  END IF

  IF ( EigenMode ) THEN 
     IF ( ASSOCIATED( Solver % Variable % EigenValues ) ) THEN
       HarmVar % EigenValues  => Solver % Variable % Eigenvalues
       HarmVar % EigenVectors => Solver % Variable % EigenVectors
     END IF
  END IF


  ! Now change the pointers such that when we visit the linear solver
  ! the system will automatically be solved as complex
  Solver % Variable => HarmVar
  Solver % Matrix => Aharm

  IF(AnyDirichlet) THEN
    IF(ParEnv % PEs>1) CALL ParallelInitMatrix(Solver, Aharm)
    CALL EnforceDirichletConditions( Solver, Aharm, b )
  END IF

  ! Save the original matrix and variable in Ematrix and Evar
  Solver % Matrix % Ematrix => SaveMatrix
  Solver % Variable % Evar => SaveVar    

  ! Eliminate cyclic dependence that is a bummer when deallocating stuff
  ! We are toggling {Are,Aharm} in {Solver % Matrix, Solver % Matrix % Ematrix}
  NULLIFY( Solver % Matrix % EMatrix % Ematrix )

!------------------------------------------------------------------------------
END SUBROUTINE ChangeToHarmonicSystem
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Eliminate linear restriction only using the ListMatrix structure.
!> This is placed in a separate routine such that it can be called
!> when solving with and without block matrix being active.
!------------------------------------------------------------------------------
SUBROUTINE EliminateLinearRestriction( StiffMatrix, ForceVector, RestMatrix, &
    CollectionMatrix, Solver, CopyStiffMatrix, ExportUsePerm, ExportUseIPerm, ExportUseDiag )
  IMPLICIT NONE
  TYPE(Matrix_t) :: StiffMatrix
  REAL(KIND=dp) :: ForceVector(:) 
  TYPE(Matrix_t), POINTER :: RestMatrix
  TYPE(Matrix_t) :: CollectionMatrix  
  TYPE(Solver_t) :: Solver
  LOGICAL, OPTIONAL :: CopyStiffMatrix
  INTEGER, POINTER, OPTIONAL :: ExportUsePerm(:), ExportUseIPerm(:)
  REAL(KIND=dp), POINTER, OPTIONAL :: ExportUseDiag(:)
  
  INTEGER :: m,n,i,j,k,l,ix,p,q,Loop
  INTEGER, ALLOCATABLE, TARGET :: SlavePerm(:),MasterPerm(:),SlaveIPerm(:),MasterIPerm(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: SlaveDiag(:), MasterDiag(:), DiagDiag(:)
  INTEGER, POINTER :: UsePerm(:), UseIPerm(:)
  REAL(KIND=dp), POINTER :: UseDiag(:)
  REAL(KIND=dp) :: scl
  REAL(KIND=dp), POINTER :: TVals(:), Vals(:)
  REAL(KIND=dp), POINTER :: CollectionVector(:), RestVector(:)
  TYPE(ListMatrix_t), POINTER :: Lmat(:)
  TYPE(Matrix_t), POINTER :: Xmat, Tmat
  TYPE(ListMatrixEntry_t), POINTER :: cTmp
  LOGICAL :: Found, EliminateSlave, EliminateFromMaster, UseTranspose
  TYPE(ValueList_t), POINTER :: Params
  CHARACTER(*), PARAMETER :: Caller = 'EliminateLinearRestriction'
  

  CALL Info(Caller,'Eliminating Constraints from CollectionMatrix',Level=12)

  Params => Solver % Values


  EliminateSlave = ListGetLogical( Params, 'Eliminate Slave',Found )
  EliminateFromMaster = ListGetLogical( Params, 'Eliminate From Master',Found )

  UseTranspose = ListGetLogical(Params, 'Use Transpose values', Found)
  IF( UseTranspose ) THEN
    CALL Info(Caller,'Using transpose values in elimination',Level=15)            
  END IF

  
  n = StiffMatrix % NumberOfRows
  m = RestMatrix % NumberOfRows

  RestVector => NULL()
  IF(ASSOCIATED(RestMatrix)) RestVector => RestMatrix % RHS

  IF(.NOT. ASSOCIATED(CollectionMatrix % Rhs) ) THEN
    ALLOCATE(CollectionMatrix % Rhs(n) )
    CollectionMatrix % Rhs = 0.0_dp
  END IF    
  CollectionVector => CollectionMatrix % RHS  

  ! We may optionally ask that the stiffness matrix is copied to the base.
  IF( PRESENT(CopyStiffMatrix)) THEN
    IF(CopyStiffMatrix) THEN
      DO i=StiffMatrix % NumberOfRows,1,-1
        DO j=StiffMatrix % Rows(i+1)-1,StiffMatrix % Rows(i),-1
          CALL AddToMatrixElement( CollectionMatrix, &
              i, StiffMatrix % Cols(j), StiffMatrix % Values(j) )
        END DO
        CollectionVector(i) = CollectionVector(i) + ForceVector(i)
      END DO
    END IF
  END IF
    
  
  ALLOCATE(SlaveDiag(m),MasterDiag(m),SlavePerm(n),MasterPerm(n),&
      SlaveIPerm(m),MasterIPerm(m),DiagDiag(m))
  SlavePerm  = 0; SlaveIPerm  = 0; 
  MasterPerm = 0; MasterIPerm = 0
  SlaveDiag = 0.0_dp; MasterDiag = 0.0_dp
  DiagDiag = 0.0_dp
  
  Tvals => RestMatrix % TValues
  IF (.NOT.ASSOCIATED(Tvals)) Tvals => RestMatrix % Values 

  ! Extract diagonal entries for constraints:
  !------------------------------------------
  CALL Info(Caller,'Extracting diagonal entries for constraints',Level=15)

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

  IF(InfoActive(25)) THEN
    PRINT *,'SlaveSum:',SUM(SlaveDiag)
    PRINT *,'MasterSum:',SUM(MasterDiag) 
    PRINT *,'SlaveSum abs:',SUM(ABS(SlaveDiag))
    PRINT *,'MasterSum abs:',SUM(ABS(MasterDiag))
  END IF

  IF(EliminateFromMaster) THEN
    CALL Info(Caller,'Eliminating from master',Level=15)      
    UsePerm  => MasterPerm 
    UseDiag  => MasterDiag
    UseIPerm => MasterIPerm 
  ELSE
    CALL Info(Caller,'Eliminating from slave',Level=15)            
    UsePerm  => SlavePerm
    UseDiag  => SlaveDiag
    UseIPerm => SlaveIPerm
  END IF
      
  IF(UseTranspose) THEN
    Vals => Tvals
  ELSE
    Vals => RestMatrix % Values
  END IF

  ! The rest is done in List Matrix format so move to that in case not yet!
  IF( CollectionMatrix % FORMAT /= MATRIX_LIST ) THEN
    CALL List_ToListMatrix(CollectionMatrix)
  END IF

  ! Replace elimination equations by the constraints (could done be as a postprocessing
  ! step, if eq's totally eliminated from linsys.)
  ! ----------------------------------------------------------------------------------
  CALL Info(Caller,'Deleting rows from equation to be eliminated',Level=15)

  Lmat => CollectionMatrix % ListMatrix
  DO m=1,RestMatrix % NumberOfRows
    i = UseIPerm(m)
    CALL List_DeleteRow(Lmat, i, Keep=.TRUE.)
  END DO

  CALL Info(Caller,'Copying rows from constraint matrix to eliminate dofs',Level=15)
  DO m=1,RestMatrix % NumberOfRows
    i = UseIPerm(m)
    DO l=RestMatrix % Rows(m+1)-1, RestMatrix % Rows(m), -1
      j = RestMatrix % Cols(l)

      ! skip l-coefficient entries, handled separately afterwards:
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
    Found = (q>1)

    Tmat => Xmat
    IF(Found) THEN
      Loop = Loop + 1
      CALL Info(Caller,'Recursive elimination round: '//I2S(Loop),Level=15)

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
  CALL Info(Caller,'Eliminating Lagrange Coefficients',Level=15)

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
        scl = -Tvals(j) / UseDiag(m)
      ELSE
        k = UseIPerm(k-n)
        scl = -Tvals(j) / UseDiag(m)
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
    CALL Info(Caller,'Eliminate slave dofs using constraint equations',Level=15)

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

  ! Optimize bandwidth, if needed:
  ! ------------------------------
  IF(EliminateFromMaster) THEN
    CALL Info(Caller,'Optimizing bandwidth after elimination',Level=15)
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

  IF( PRESENT(ExportUsePerm) ) THEN
    CALL Info(Caller,'Export UsePerm outside elimination',Level=20) 
    ALLOCATE(ExportUsePerm(SIZE(UsePerm)))
    ExportUsePerm = UsePerm
  END IF
  IF( PRESENT(ExportUseIPerm) ) THEN
    CALL Info(Caller,'Export UseIPerm outside elimination',Level=20) 
    ALLOCATE(ExportUseIPerm(SIZE(UseIPerm)))
    ExportUseIPerm = UseIPerm
  END IF
  IF( PRESENT(ExportUseDiag) ) THEN
    CALL Info(Caller,'Export UseDiag outside elimination',Level=20) 
    ALLOCATE(ExportUseDiag(SIZE(UseDiag)))
    ExportUseDiag = UseDiag
  END IF

  IF(PRESENT(CopyStiffMatrix)) THEN
    IF(CopyStiffMatrix) CALL List_ToCRSMatrix(CollectionMatrix)
  END IF
      
  CALL Info(Caller,'Finished Eliminating Restrictions',Level=12)

END SUBROUTINE EliminateLinearRestriction

  

!------------------------------------------------------------------------------
!>  This subroutine will solve the system with some linear restriction.
!>  The restriction matrix is assumed to be in the ConstraintMatrix-field of 
!>  the StiffMatrix. The restriction vector is the RHS-field of the
!>  ConstraintMatrix.
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE SolveWithLinearRestriction( StiffMatrix, ForceVector, &
    Solution, Norm, DOFs, Solver )
!------------------------------------------------------------------------------  
  IMPLICIT NONE
  TYPE(Matrix_t), POINTER :: StiffMatrix !< Linear equation matrix information. 
  REAL(KIND=dp),TARGET :: ForceVector(:) !< The right hand side of the linear equation
  REAL(KIND=dp),TARGET :: Solution(:)    !< Previous solution as input, new solution as output.
  REAL(KIND=dp) :: Norm                  !< The L2 norm of the solution.
  INTEGER :: DOFs                        !< Number of degrees of freedom of the equation.
  TYPE(Solver_t), TARGET :: Solver       !< Linear equation solver options.
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: SolverPointer
  TYPE(Matrix_t), POINTER :: CollectionMatrix, RestMatrix, AddMatrix, RestMatrixTranspose 
  REAL(KIND=dp), POINTER CONTIG :: CollectionVector(:), RestVector(:), AddVector(:) 
  REAL(KIND=dp), POINTER  :: MultiplierValues(:), pSol(:),DiagScaling(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: CollectionSolution(:)
  INTEGER :: NumberOfRows, NumberOfValues, MultiplierDOFs, istat, NoEmptyRows 
  INTEGER :: i, j, k, l, m, n, p,q, ix, Loop, colj, nIter
  TYPE(Variable_t), POINTER :: MultVar, iterV
  REAL(KIND=dp) :: scl, rowsum, Relax, val
  LOGICAL :: Found, ExportMultiplier, NotExplicit, Refactorize, EnforceDirichlet, &
      NonEmptyRow, ComplexSystem, ConstraintScaling, UseTranspose, EliminateConstraints, &
      SkipConstraints, ResidualMode
  INTEGER, POINTER :: UseIPerm(:), UsePerm(:)
  REAL(KIND=dp), POINTER :: UseDiag(:) 
  LOGICAL  :: Parallel, UseTreeGauge, NeedMassDampValues, DoOwnScaling
  LOGICAL, ALLOCATABLE :: TrueDof(:)
  INTEGER, ALLOCATABLE :: Iperm(:)
  REAL(KIND=dp) :: t0,rt0,st,rst
  CHARACTER(:), ALLOCATABLE :: str,MultiplierName
  TYPE(ValueList_t), POINTER :: Params
  CHARACTER(*), PARAMETER :: Caller = 'SolveWithLinearRestriction'

  TYPE(ParEnv_t), POINTER :: ParEnvSave

  SAVE MultiplierValues, SolverPointer
  
!------------------------------------------------------------------------------
  CALL Info( Caller, ' ', Level=12 )
  ParEnvSave => ParEnv

  SolverPointer => Solver  
  Params => Solver % Values

  t0 = CPUTime()
  rt0 = RealTime()
    
  Parallel = Solver % Parallel

  ResidualMode = ListGetLogical( Params,'Restriction System Residual Mode',Found )    
  iterV => VariableGet(Solver % Mesh % Variables,'nonlin iter')
  nIter = NINT(iterV % Values(1))
    
  NotExplicit = ListGetLogical(Params,'No Explicit Constrained Matrix',Found)
  IF(.NOT. Found) NotExplicit=.FALSE.

  NeedMassDampValues = ListGetLogical(Params, 'Eigen Analysis', Found )

  RestMatrix => NULL()
  IF(.NOT.NotExplicit) RestMatrix => StiffMatrix % ConstraintMatrix
  RestVector => Null()
  IF(ASSOCIATED(RestMatrix)) RestVector => RestMatrix % RHS

  AddMatrix => StiffMatrix % AddMatrix
  AddVector => NULL()
  IF(ASSOCIATED(AddMatrix)) AddVector => AddMatrix % RHS

  EliminateConstraints = ListGetLogical( Params, 'Eliminate Linear Constraints', Found)
  
  NumberOfRows = StiffMatrix % NumberOfRows
  
  CollectionMatrix => StiffMatrix % CollectionMatrix
  Refactorize = ListGetLogical(Params,'Linear System Refactorize',Found)
  IF(.NOT.Found) THEN
    Refactorize = .NOT. ( ResidualMode .AND. nIter > 1) 
  END IF

  IF(ASSOCIATED(CollectionMatrix)) THEN
    IF(Refactorize .AND. .NOT.NotExplicit) THEN
      CALL Info( Caller,'Freeing previous collection matrix structures',Level=10)
      CALL FreeMatrix(CollectionMatrix)
      CollectionMatrix => NULL()
      CALL Info( Caller,'Refactoring requested, creating fully new matrix',Level=10)
    ELSE
      CALL Info( Caller,'Trying to keep previous collection matrix structures',Level=10)
    END IF
  END IF

  IF(.NOT.ASSOCIATED(CollectionMatrix)) THEN
    CollectionMatrix => AllocateMatrix()
    CollectionMatrix % FORMAT = MATRIX_LIST
  ELSE
    DEALLOCATE(CollectionMatrix % RHS)
    CollectionMatrix % Values = 0.0_dp

    IF(NeedMassDampValues) THEN
      IF(ASSOCIATED(CollectionMatrix % MassValues)) CollectionMatrix % MassValues = 0.0_dp
      IF(ASSOCIATED(CollectionMatrix % DampValues)) CollectionMatrix % DampValues = 0.0_dp
    END IF
  END IF
  IF(NotExplicit) CollectionMatrix % ConstraintMatrix => StiffMatrix % ConstraintMatrix  
  
  NumberOfRows = StiffMatrix % NumberOfRows
  IF(ASSOCIATED(AddMatrix)) NumberOfRows = MAX(NumberOfRows,AddMatrix % NumberOfRows)
  IF(ASSOCIATED(RestMatrix)) THEN
    IF(.NOT.EliminateConstraints) NumberOfRows = NumberOFRows + RestMatrix % NumberOfRows
  END IF

  ALLOCATE( CollectionMatrix % RHS( NumberOfRows ), &
       CollectionSolution( NumberOfRows ), STAT = istat )
  IF ( istat /= 0 ) CALL Fatal( Caller, 'Memory allocation error.' )

  CollectionVector => CollectionMatrix % RHS
  CollectionVector = 0.0_dp
  CollectionSolution = 0.0_dp

  ComplexSystem = StiffMatrix % COMPLEX .OR. &
      ListGetLogical(Params,'Linear System Complex', Found )
  
!------------------------------------------------------------------------------
! If multiplier should be exported,  allocate memory and export the variable.
!------------------------------------------------------------------------------

  ExportMultiplier = ListGetLogical(Params, 'Export Lagrange Multiplier', Found )
  IF ( ExportMultiplier ) THEN
    MultiplierName = LagrangeMultiplierName( Solver )
    MultVar => VariableGet(Solver % Mesh % Variables, MultiplierName)
    j = 0
    IF(ASSOCIATED(RestMatrix)) j = RestMatrix % NumberofRows
    IF(ASSOCIATED(AddMatrix))  j = j+MAX(0,AddMatrix % NumberofRows-StiffMatrix % NumberOfRows)

    IF ( .NOT. ASSOCIATED(MultVar) ) THEN
      CALL Info(Caller,'Creating variable for Lagrange multiplier',Level=8)
      ALLOCATE( MultiplierValues(j), STAT=istat )
      IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error.')

      MultiplierValues = 0.0_dp
      IF( ComplexSystem ) THEN
        CALL VariableAddVector(Solver % Mesh % Variables, Solver % Mesh, SolverPointer, &
            MultiplierName, 2, MultiplierValues)               
      ELSE
        CALL VariableAdd(Solver % Mesh % Variables, Solver % Mesh, SolverPointer, &
            MultiplierName, 1, MultiplierValues)      
      END IF
      MultVar => VariableGet(Solver % Mesh % Variables, MultiplierName)
    END IF

    IF( InfoActive( 20 ) ) THEN
      CALL VectorValuesRange(MultVar % Values,SIZE(MultVar % Values),TRIM(MultVar % Name))
    END IF

    MultiplierValues => MultVar % Values

    IF (j > SIZE(MultiplierValues)) THEN
      CALL Info(Caller,'Increasing Lagrange multiplier size to: '//I2S(j),Level=8)
      ALLOCATE(MultiplierValues(j)); MultiplierValues=0._dp       
      MultiplierValues(1:SIZE(MultVar % Values)) = MultVar % Values     

      ! If the Lagrange variable includes history also change its size.
      IF( ASSOCIATED( MultVar % PrevValues ) ) THEN
        MultVar % Values = MultVar % PrevValues(:,1)
        DEALLOCATE( MultVar % PrevValues )
        ALLOCATE( MultVar % PrevValues(j,1) )
        MultVar % PrevValues = 0.0_dp
        MultVar % PrevValues(:,1) = MultVar % Values
      END IF

      DEALLOCATE(MultVar % Values)
      MultVar % Values => MultiplierValues
    END IF

    IF( InfoActive(25) ) THEN
      CALL VectorValuesRange(MultVar % values,SIZE(MultVar % values),'MultVar')
    END IF
  ELSE
    MultiplierValues => NULL()
  END IF

  UseTreeGauge = ListGetlogical(Params, 'Use Tree Gauge', Found )

!------------------------------------------------------------------------------
! Put the RestMatrix to lower part of CollectionMatrix
!------------------------------------------------------------------------------

  EnforceDirichlet = ListGetLogical(Params, 'Enforce Exact Dirichlet BCs',Found)
  IF(.NOT.Found) EnforceDirichlet = .TRUE.
  EnforceDirichlet = EnforceDirichlet .AND. ALLOCATED(StiffMatrix % ConstrainedDOF)

  UseTranspose = ListGetLogical(Params, 'Use Transpose values', Found)
  IF( UseTranspose ) THEN
    CALL Info(Caller,'Using transpose values in elimination',Level=15)            
  END IF
    
  CALL Info(Caller,'Number of Rows / Nonzeros in original matrix: '&
      //I2S(StiffMatrix % NumberOfRows)//' / '//I2S(SIZE(StiffMatrix % Values)),Level=22) 
  
  IF(ASSOCIATED(RestMatrix) .AND. .NOT. EliminateConstraints) THEN

    CALL Info(Caller,'Adding ConstraintMatrix into CollectionMatrix',Level=10)

    CALL Info(Caller,'Number of Rows / Nonzeros in constraint matrix: '&
      //I2S(RestMatrix % NumberOfRows)//' / '//I2S(SIZE(RestMatrix % Values)),Level=12) 

    NoEmptyRows = 0
    ConstraintScaling = ListGetLogical(Params, 'Constraint Scaling',Found)
    IF(ConstraintScaling) THEN
      rowsum = ListGetConstReal(Params, 'Constraint Scale', Found)
      IF(Found) RestMatrix % Values = RestMatrix % Values * rowsum
    END IF

    ALLOCATE( iperm(SIZE(Solver % Variable % Perm)) )
    iperm = 0
    DO i=1,SIZE(Solver % Variable % Perm)
      IF ( Solver % Variable % Perm(i)>0) Iperm(Solver % Variable % Perm(i))=i
    END DO

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
      NonEmptyRow = .FALSE.

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
          ! Skip non-positive column indexes, why should there be any? 
          colj = RestMatrix % Cols(j)
          IF( colj <= 0 ) CYCLE

          ! Complex system requires that we have exact 2x2 block structure. Don't spoil that.
          IF ( .NOT. ComplexSystem ) THEN
            IF( ABS(RestMatrix % Values(j)) < EPSILON(1._dp)*rowsum ) CYCLE
          END IF

          ! If we have Dirichlet condition set for the matrix use that directly and do not add
          ! stuff to the row that would spoil the condition. 
          Found = .TRUE.
          IF (EnforceDirichlet .AND. colj <= StiffMatrix % NumberOfRows) THEN
            Found = .NOT. StiffMatrix % ConstrainedDOF(colj)
          END IF
            
          IF(Found) THEN
            IF (ASSOCIATED(RestMatrix % TValues)) THEN
              val = RestMatrix % TValues(j)
            ELSE
              val = RestMatrix % Values(j)
            END IF              
            CALL AddToMatrixElement( CollectionMatrix, colj, k, val ) 

            ! Add the Transpose part
            IF (UseTranspose .AND. ASSOCIATED(RestMatrix % TValues)) THEN
              val = RestMatrix % TValues(j)
            ELSE
              val = RestMatrix % Values(j)
            END IF

            ! Only add the transpose when it is associated to the unknowns of the initial matrix.
            ! Otherwise the entries related to largrange multipliers would be multiplied by factor 2!
            IF( colj <= StiffMatrix % NumberOfRows ) THEN            
              CALL AddToMatrixElement( CollectionMatrix, k, colj, val ) 
              NonEmptyRow = NonEmptyRow .OR. val /= 0
            END IF
          ELSE
            IF (UseTranspose .AND. ASSOCIATED(RestMatrix % TValues)) THEN
              val = RestMatrix % TValues(j)
            ELSE 
              val = RestMatrix % Values(j)
            END IF
            ! Use the value of the Dirichlet condition from "Dvalues"
            CollectionVector(k) = CollectionVector(k) - val * StiffMatrix % DValues(colj)  
          END IF
        END DO
      END IF
 
      Found = .TRUE.
      IF (EnforceDirichlet) THEN
        IF(ASSOCIATED(RestMatrix % InvPerm)) THEN
          l = RestMatrix % InvPerm(i)
          IF(l>0) THEN
            l = MOD(l-1,StiffMatrix % NumberOfRows)+1
            IF(StiffMatrix % ConstrainedDOF(l)) THEN
              l = iperm((l-1)/Solver % Variable % DOFs+1) 
              IF (l<=Solver % Mesh % NumberOfNodes) THEN
                Found = .FALSE.
                CALL ZeroRow(CollectionMatrix,k)
                CollectionVector(k) = 0
                CALL SetMatrixElement(CollectionMatrix,k,k,1._dp)
              END IF
            END IF
          END IF
        END IF
      END IF

      ! If there is no matrix entry, there can be no non-zero r.h.s.
      IF ( Found ) THEN
        IF( .NOT.NonEmptyRow ) THEN
          NoEmptyRows = NoEmptyRows + 1
          CollectionVector(k) = 0._dp
!          might not be the right thing to do in parallel!!
          IF(UseTreeGauge) THEN
            CALL SetMatrixElement( CollectionMatrix,k,k,1._dp )
          END IF 
        ELSE
          IF( ASSOCIATED( RestVector ) ) CollectionVector(k) = CollectionVector(k) + RestVector(i)
        END IF
      END IF
    END DO

    IF( NoEmptyRows > 0 ) THEN
      CALL Info(Caller,'Constraint Matrix has '//I2S(NoEmptyRows)// &
          ' empty rows out of '//I2S(RestMatrix % NumberOfRows),Level=6 )
    END IF

    CALL Info(Caller,'Finished Adding ConstraintMatrix',Level=12)
  END IF

!------------------------------------------------------------------------------
! Put the AddMatrix to upper part of CollectionMatrix
!------------------------------------------------------------------------------
  IF(ASSOCIATED(AddMatrix)) THEN
    CALL Info(Caller,'Adding AddMatrix into CollectionMatrix',Level=12)

    CALL Info(Caller,'Number of Rows / Nonzeros in additional matrix: '&
      //I2S(AddMatrix % NumberOfRows)//' / '//I2S(SIZE(AddMatrix % Values)),Level=12) 
    
    DO i=AddMatrix % NumberOfRows,1,-1

      Found = .TRUE.
      IF (EnforceDirichlet .AND. i<=StiffMatrix % NumberOFRows) THEN
        Found = .NOT.StiffMatrix % ConstrainedDOF(i)
      END IF
        
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
    CALL Info(Caller,'Finished Adding AddMatrix',Level=12)
  END IF

!------------------------------------------------------------------------------
! Put the StiffMatrix to upper part of CollectionMatrix
!------------------------------------------------------------------------------
  CALL Info(Caller,'Adding Stiffness Matrix into CollectionMatrix',Level=12)

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
! necessarily biorthogonal constraint equation test functions.
!------------------------------------------------------------------------------
  IF (ASSOCIATED(RestMatrix) .AND. EliminateConstraints) THEN
    IF ( ExportMultiplier ) THEN
      ! With the multiplier active we need to use it also for elimination in case the
      ! constraint shares some dofs with the multiplier. 
      CALL EliminateLinearRestriction( StiffMatrix, ForceVector, RestMatrix, &
          CollectionMatrix, Solver, ExportUseIPerm = UseIPerm, ExportUseDiag = UseDiag )
    ELSE
      CALL EliminateLinearRestriction( StiffMatrix, ForceVector, RestMatrix, &
        CollectionMatrix, Solver )
    END IF
  END IF
  
  IF(CollectionMatrix % FORMAT==MATRIX_LIST) THEN
    CALL Info(Caller,'Reverting CollectionMatrix back to CRS matrix',Level=10)
    CALL List_toCRSMatrix(CollectionMatrix)
  END IF
    
  ! CRS-format matrix needed here
  IF ( NeedMassDampValues ) THEN  ! Doesn't work with constraints, "AddMatrix" only !!
    CALL CopyMassDampValues(CollectionMatrix, StiffMatrix, AddMatrix)
  END IF
  
  CALL Info( Caller, 'CollectionMatrix done', Level=12 )

!------------------------------------------------------------------------------
! Assign values to CollectionVector
!------------------------------------------------------------------------------

  j = StiffMatrix % NumberOfRows  
  CollectionSolution(1:j) = Solution(1:j)
  
  i = StiffMatrix % NumberOfRows+1
  j = SIZE(CollectionSolution)
  IF( j >= i) THEN
    CollectionSolution(i:j) = 0._dp
    IF(ExportMultiplier) CollectionSolution(i:j) = MultiplierValues(1:j-i+1)
  END IF

  IF( InfoActive(30) ) THEN
    pSol => CollectionSolution
    CALL VectorValuesRange(pSol,j,'CollectionSolution')           
  END IF
  
  CollectionMatrix % ExtraDOFs = CollectionMatrix % NumberOfRows - &
                  StiffMatrix % NumberOfRows

  CollectionMatrix % ParallelDOFs = 0
  IF(ASSOCIATED(AddMatrix)) &
    CollectionMatrix % ParallelDOFs = MAX(AddMatrix % NumberOfRows - &
            StiffMatrix % NumberOfRows,0)
    
  CALL Info( Caller, 'CollectionVector done', Level=12 )

!------------------------------------------------------------------------------
! Solve the Collection-system 
!------------------------------------------------------------------------------

! Collectionmatrix % Complex = StiffMatrix % Complex


  CollectionMatrix % Comm = StiffMatrix % Comm

  
  st  = CPUTime() - t0;
  rst = RealTime() - rt0  
  WRITE(Message,'(a,f8.2,f8.2,a)') 'Collection matrix creation time (CPU,REAL): ',st,rst,' (s)'
  CALL Info(Caller,Message,Level=6)    
  
  i = CollectionMatrix % NumberOfRows - StiffMatrix % NumberOfRows
  j = SIZE(CollectionMatrix % Values) - SIZE(StiffMatrix % Values )
  CALL Info(Caller,'Collection matrix increased with '//I2S(i)//&
      ' rows and '//I2S(j)//' non-zeros',Level=8)
    
  IF( InfoActive( 30 ) ) THEN
    CALL VectorValuesRange(CollectionMatrix % Values,SIZE(CollectionMatrix % Values),'A')       
    CALL VectorValuesRange(CollectionMatrix % rhs,SIZE(CollectionMatrix % rhs),'b')       
  END IF
      
  IF( ResidualMode ) THEN
    BLOCK
      REAL(KIND=dp), POINTER :: Res(:)      
      ! If residual mode is requested make change of variables:
      ! Ax=b -> Adx = b-Ax0 = r
      IF( niter > 1 ) THEN
        CALL Info(Caller,'Changing the equation to residual based mode',Level=10)
        ALLOCATE( Res(SIZE(CollectionSolution)) ) 
        CALL LinearSystemResidual( CollectionMatrix, CollectionVector, CollectionSolution, res )
        CollectionVector = Res
        CollectionSolution = 0.0_dp
        DEALLOCATE(Res)
      END IF
    END BLOCK
  END IF
  
  ! We may want to skip ComputeChange including the constraints if we use certain other options
  SkipConstraints = ResidualMode .OR. &
      ListGetLogical( Params, 'Nonlinear System Convergence Without Constraints',Found ) .OR. &
      ListGetLogical( Params, 'NonLinear System Consistent Norm',Found )   
  str = ListGetString( Params, 'NonLinear System Convergence Measure',Found )
  IF( str == 'solution' ) THEN
    SkipConstraints = .TRUE.
    CALL Info(Caller,&
        'Nonlinear system convergence measure == "solution" must skip constraints',Level=10)
  END IF
  IF( SkipConstraints ) THEN
    CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)
    CALL ListAddLogical( Params,'Skip Advance Nonlinear iter',.TRUE.)
  END IF

  DoOwnScaling = ListGetLogical( Params,'Linear System Scaling',Found)
  IF(.NOT. Found) DoOwnScaling = .TRUE.
  IF(.NOT. ResidualMode) DoOwnScaling = .FALSE.
  IF(DoOwnScaling) THEN   
    CALL Info(Caller,'Performing special scaling with constraints',Level=10)
    DiagScaling => CollectionMatrix % DiagScaling
    IF(Niter == 1 ) THEN
      IF(.NOT. ASSOCIATED(DiagScaling) ) THEN
        ALLOCATE( DiagScaling(SIZE(CollectionVector)))
        CollectionMatrix % DiagScaling => DiagScaling
      END IF

      ! Should we scale only part or the full matrix? 
      IF(.FALSE.) THEN
        DiagScaling = 1.0_dp
        StiffMatrix % DiagScaling => DiagScaling
        ! Just build the scaling matrix using only the original stiffness matrix.
        CALL ScaleLinearSystem(Solver,StiffMatrix,ApplyScaling=.FALSE.)     
        CollectionMatrix % ScalingMethod = StiffMatrix % ScalingMethod    
        StiffMatrix % DiagScaling => NULL()      
      ELSE
        CALL ScaleLinearSystem(Solver,CollectionMatrix,ApplyScaling=.FALSE.)     
      END IF        
    END IF
    
    CALL ScaleLinearSystem(Solver,CollectionMatrix,CollectionVector,&
        CollectionSolution,DiagScaling=CollectionMatrix % DiagScaling)
    CALL ListAddLogical( Params,'Linear System Skip Scaling',.TRUE. ) 
  END IF
  
  
  !IF( ListGetLogical( Params,'Linear System Save',Found ) ) THEN        
  !  CALL SaveLinearSystem( Solver, CollectionMatrix,'RestrictedMat')
  !END IF
  

  CALL Info(Caller,'Now solving the linear system with constraints!',Level=10)
  Collectionmatrix % DGMatrix = StiffMatrix %  DGMatrix
  CALL SolveLinearSystem( CollectionMatrix, CollectionVector, &
     CollectionSolution, Norm, DOFs, Solver, StiffMatrix )


  IF(DoOwnScaling) THEN
    CALL BackScaleLinearSystem( Solver,CollectionMatrix,CollectionVector,&
        CollectionSolution,CollectionMatrix % DiagScaling)
    CALL ListAddLogical( Params,'Linear System Skip Scaling',.FALSE. ) 
  END IF

  !-------------------------------------------------------------------------------
  ! For restricted systems study the norm without some block components.
  ! For example, excluding gauge constraints may give valuable information
  ! of the real accuracy of the unconstrained system. Currently just for info.
  !-------------------------------------------------------------------------------
  IF( ListGetLogical( Params,'Restricted System Norm',Found ) ) THEN
    ALLOCATE( TrueDof( CollectionMatrix % NumberOfRows ) )
    TrueDof = .TRUE.
    
    Norm = LinearSystemMaskedResidualNorm( CollectionMatrix, CollectionVector, &
        CollectionSolution, TrueDof, TrueDof )
    
    WRITE( Message,'(A,ES13.6)') 'Residual norm of the original system:',Norm
    CALL Info(Caller,Message, Level = 5 )
    
    IF( ListGetLogical( Params,'Restricted System Norm Skip Nodes',Found ) ) THEN
      i = 1
      j = MAXVAL( Solver % Variable % Perm(1:Solver % Mesh % NumberOfNodes) )
      CALL Info(Caller,'Skipping nodal dof range: '//I2S(i)//'-'//I2S(j),Level=8)
      TrueDof(i:j) = .FALSE.
    END IF

    IF( ListGetLogical( Params,'Restricted System Norm Skip Constraints',Found ) ) THEN
      i = StiffMatrix % NumberOfRows + 1
      j = CollectionMatrix % NumberOfRows      
      CALL Info(Caller,'Skipping constraints dof range: '&
          //I2S(i)//'-'//I2S(j),Level=8)
      TrueDof(i:j) = .FALSE.
    END IF
    
    Norm = LinearSystemMaskedResidualNorm( CollectionMatrix, CollectionVector, &
        CollectionSolution, TrueDof, TrueDof )
    
    WRITE( Message,'(A,ES13.6)') 'Residual norm of the masked system:',Norm
    CALL Info(Caller,Message, Level = 5 )
    
    DEALLOCATE( TrueDof )
  END IF
    
  
!------------------------------------------------------------------------------
! Separate the solution from CollectionSolution
!------------------------------------------------------------------------------
    CALL Info(Caller,'Picking solution from collection solution',Level=10)

    j = StiffMatrix % NumberOfRows

    IF( ResidualMode .AND. nIter > 1) THEN
      Solution(1:j) = Solution(1:j) + CollectionSolution(1:j)
    ELSE
      Solution(1:j) = CollectionSolution(1:j)
    END IF
    
    IF ( ExportMultiplier ) THEN
      CALL Info(Caller,'Separating Lagrange multiplier from collection solution',Level=10)
      
      IF(ASSOCIATED(RestMatrix) .AND. EliminateConstraints) THEN        
        ! Compute eliminated l-coefficient values:
        ! ---------------------------------------
        MultiplierValues = 0.0_dp
        DO i=1,RestMatrix % NumberOfRows
          scl = 1._dp / UseDiag(i)
          m = UseIPerm(i)
          MultiplierValues(i) = scl * ForceVector(m)
          DO j=StiffMatrix % Rows(m), StiffMatrix % Rows(m+1)-1
            MultiplierValues(i) = MultiplierValues(i) - &
                scl * StiffMatrix % Values(j) * Solution(StiffMatrix % Cols(j))
          END DO
        END DO

        DEALLOCATE( UseIPerm, UseDiag )
      ELSE
        i = StiffMatrix % NumberOfRows
        j=0
        IF(ASSOCIATED(RestMatrix)) j = RestMatrix % NumberOfRows
        IF(ASSOCIATED(AddMatrix)) &
            j=j+MAX(0,AddMatrix % NumberOfRows - StiffMatrix % NumberOFRows)

        Relax = ListGetCReal( Params,'Lagrange Multiplier Relaxation Factor', Found )
        IF( ResidualMode .AND. nIter > 1 ) THEN
          IF( Found ) THEN          
            MultiplierValues(1:j) = MultiplierValues(1:j) + &
                Relax * CollectionSolution(i+1:i+j)
          ELSE
            MultiplierValues(1:j) = MultiplierValues(1:j) + CollectionSolution(i+1:i+j)
          END IF
        ELSE
          IF( Found ) THEN          
            MultiplierValues(1:j) = (1-Relax) * MultiplierValues(1:j) + &
                Relax * CollectionSolution(i+1:i+j)
          ELSE       
            MultiplierValues(1:j) = CollectionSolution(i+1:i+j)
          END IF
        END IF
                
      END IF      
    END IF

      
!------------------------------------------------------------------------------

    IF( SkipConstraints ) THEN
      CALL ListAddLogical( Params,'Skip Advance Nonlinear iter',.FALSE.)
      CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.FALSE.)
      CALL ComputeChange(Solver,.FALSE.,StiffMatrix % NumberOfRows,Matrix=StiffMatrix,Rhs=ForceVector)
    END IF
        
    DEALLOCATE(CollectionSolution)
    CollectionMatrix % ConstraintMatrix => NULL()
    StiffMatrix % CollectionMatrix => CollectionMatrix

    ParEnv => ParEnvSave
    
    CALL Info( Caller, 'All done', Level=10 )
CONTAINS


!------------------------------------------------------------------------------
   SUBROUTINE CopyMassDampValues(A,B,C)
!------------------------------------------------------------------------------
     TYPE(Matrix_t)  :: A,B,C
!------------------------------------------------------------------------------
     INTEGER :: i,j,n
!------------------------------------------------------------------------------
     REAL(KIND=dp), POINTER :: svals(:)
!------------------------------------------------------------------------------
     n = SIZE(A % Values)

     IF(ASSOCIATED(B % MassValues) .OR. ASSOCIATED(C % MassValues)) THEN
       IF(.NOT.ASSOCIATED(A % MassValues)) THEN
         ALLOCATE(A % MassValues(n)); A % MassValues = 0._dp
       END IF
     END IF

     IF(ASSOCIATED(B % DampValues) .OR. ASSOCIATED(C % DampValues)) THEN
       IF(.NOT.ASSOCIATED(A % DampValues)) THEN
         ALLOCATE(A % DampValues(n)); A % DampValues = 0._dp
       END IF
     END IF

     svals => A % Values

     IF(ASSOCIATED(B % MassValues)) THEN
       A % Values => A % MassValues
       DO i=B % NumberOfRows,1,-1
         DO j=B % Rows(i+1)-1,B % Rows(i),-1
           CALL AddToMatrixElement( A, i, B % Cols(j), B % Values(j) )
         END DO
       END DO
     END IF
 
     IF(ASSOCIATED(C % MassValues)) THEN
       A % Values => A % MassValues
       DO i=C % NumberOfRows,1,-1
         DO j=C % Rows(i+1)-1,C % Rows(i),-1
           CALL AddToMatrixElement( A, i, C % Cols(j), C % Values(j) )
         END DO
       END DO
     END IF
 
     IF(ASSOCIATED(B % DampValues)) THEN
       A % Values => A % DampValues
       DO i=B % NumberOfRows,1,-1
         DO j=B % Rows(i+1)-1,B % Rows(i),-1
           CALL AddToMatrixElement( A, i, B % Cols(j), B % DampValues(j) )
         END DO
       END DO
     END IF
 
     IF(ASSOCIATED(C % DampValues)) THEN
       A % Values => A % DampValues
       DO i=C % NumberOfRows,1,-1
         DO j=C % Rows(i+1)-1,C % Rows(i),-1
           CALL AddToMatrixElement( A, i, C % Cols(j), C % DampValues(j) )
         END DO
       END DO
     END IF

     A % Values => svals
!------------------------------------------------------------------------------
   END SUBROUTINE CopyMassDampValues
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE SolveWithLinearRestriction
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Get the node from on which the controlled value should be set. 
!------------------------------------------------------------------------------
  FUNCTION GetControlNode(Mesh,Perm,Params,iControl) RESULT ( ControlNode ) 
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Perm(:)
    TYPE(ValueList_t), POINTER :: Params    
    INTEGER :: iControl
    INTEGER :: ControlNode

    INTEGER :: i,m
    REAL(KIND=dp) :: Coord(3), MinDist
    REAL(KIND=dp), POINTER :: RealWork(:,:)
    LOGICAL :: Found
    CHARACTER(:), ALLOCATABLE :: str
    CHARACTER(*), PARAMETER :: Caller = 'GetControlNode'

    str = 'Control Node Index '//I2S(iControl)                
    ControlNode = ListGetInteger( Params,str,Found )
    IF(.NOT. Found .AND. iControl == 1 ) THEN
      str = 'Control Node Index'
      ControlNode = ListGetInteger( Params,str,Found )    
    END IF
   
    IF(.NOT. Found ) THEN        
      ControlNode = -1

      Coord = 0.0_dp
      str = 'Control Node Coordinates'
      RealWork => ListGetConstRealArray( Params,str,Found )           
      IF(Found) THEN
        i = iControl
      ELSE
        str = TRIM(str)//' '//I2S(iControl)        
        RealWork => ListGetConstRealArray( Params,str,Found )                         
        i = 1
      END IF
      
      IF( Found ) THEN
        IF(SIZE(RealWork,2)==1) THEN
          m = SIZE(RealWork,1)
          Coord(1:m) = RealWork(1:m,i)
        ELSE
          m = SIZE(RealWork,2)
          Coord(1:m) = RealWork(i,1:m)
        END IF
        
        CALL FindClosestNode(Mesh,Coord,MinDist,ControlNode,ParEnv % PEs>1,Perm=Perm)
        CALL Info(Caller,'Control Node located to index: '//I2S(ControlNode),Level=6)

        ! Add the index for future rounds since it takes time to make the search every time!
        str = 'Control Node Index '//I2S(iControl)  
        CALL ListAddInteger( Params, str, ControlNode )
      END IF
    END IF

  END FUNCTION GetControlNode


  FUNCTION GetControlValue(Mesh,Params,iControl,Var,dof) RESULT ( val )
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params  
    INTEGER :: iControl
    TYPE(Variable_t), POINTER, OPTIONAL :: Var
    INTEGER, OPTIONAL :: dof
    REAL(KIND=dp) :: val
    
    TYPE(Variable_t), POINTER :: pVar
    INTEGER :: i,j
    INTEGER :: dof0
    REAL(KIND=dp) :: val0
    LOGICAL :: Found
    CHARACTER(:), ALLOCATABLE :: str, varname
    CHARACTER(*), PARAMETER :: Caller = 'GetControlValue'

    IF(.NOT. ASSOCIATED(Mesh)) CALL Fatal(Caller,'Mesh not associated!')
    IF(.NOT. ASSOCIATED(Params)) CALL Fatal(Caller,'Params not associated!')
    
    str = 'Control Variable'
    varname = ListGetString( Params, str, Found )
    IF(.NOT. Found ) THEN
      str = 'Control Variable '//I2S(iControl)
      varname = ListGetString( Params, str, Found )
    END IF
    IF(.NOT. Found ) THEN
      IF( PRESENT(Var) ) THEN
        pVar => Var
      ELSE
        CALL Fatal(Caller,'Could not find keyword for: '//TRIM(str))
      END IF
    ELSE
      pVar => VariableGet( Mesh % Variables, varname )
      IF(.NOT. ASSOCIATED(pVar) ) THEN
        CALL Fatal(Caller,'Could not find control variable: '//TRIM(varname))
      END IF
    END IF

    i = GetControlNode(Mesh,pVar % Perm,Params,iControl)
    IF(i==-1) THEN
      CALL Fatal(Caller,'Could not find control node!')
    END IF
        
    dof0 = 1   
    IF(PRESENT(dof)) THEN
      dof0 = dof
    ELSE IF( pVar % Dofs > 1) THEN
      dof0 = ListGetInteger( Params,'Control Target Component',UnfoundFatal=.TRUE.)
    END IF

    IF( i == 0 ) THEN
      val = -HUGE(val)
    ELSE
      j = pVar % dofs*(pVar % Perm(i)-1)+dof0
      val = pVar % Values(j) 
    END IF

    val = ParallelReduction(val,2)
    
    str = 'Control Target Value'        
    val0 = ListGetCReal( Params, str, Found )
    IF(.NOT. Found ) THEN
      str = 'Control Target Value '//I2S(iControl)
      val0 = ListGetCReal( Params, str, Found )
    END IF    
    val = val - val0 
  
    !PRINT *,'Control value:',val,val0,i,j

  END FUNCTION GetControlValue
    

  SUBROUTINE ApplyExplicitControl(Solver)
    TYPE(Solver_t), POINTER :: Solver

    INTEGER :: i, n
    REAL(KIND=dp), POINTER :: Fvec(:)
    LOGICAL :: Found
    TYPE(Variable_t), POINTER :: FVar
    
    n = ListGetInteger(Solver % Values,'Number Of Controls',Found )
    IF( n == 0 ) THEN
      CALL Warn('ApplyExplicitControl','Explicit control points requested but no controls available!')
      RETURN
    END IF
    
    FVar => VariableGet( Solver % Mesh % Variables,'cpar' )
    IF(.NOT. ASSOCIATED( FVar ) ) THEN
      CALL VariableAddVector( Solver % Mesh % Variables,Solver % Mesh,Solver,'cpar',n,Global=.TRUE.)
      FVar => VariableGet( Solver % Mesh % Variables,'cpar' )
    END IF

    Fvec => FVar % Values    
    DO i=1,n
      Fvec(i) = GetControlValue(Solver % Mesh,Solver % Values,i,Solver % Variable) 
    END DO

    !PRINT *,'Control values:',Fvec
       
  END SUBROUTINE ApplyExplicitControl

  
!------------------------------------------------------------------------------
!> Given the operation point and an additional r.h.s. source vector find the
!> amplitude for the latter one such that the control problem is resolved.
!> We can request a field value at given point, for example. This tries to
!> mimic some ideas of the "Smart Heater Control" of "HeatSolver" available
!> long ago. This would hopefully be applicable to wider set of modules. 
!------------------------------------------------------------------------------
  SUBROUTINE ControlLinearSystem(Solver,PreSolve)
    TYPE(Solver_t), POINTER :: Solver
    LOGICAL :: PreSolve

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: A    
    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), POINTER :: x0(:),b(:),BulkRhsSave(:),dr(:),r0(:),dy(:),y0(:)
    REAL(KIND=dp), ALLOCATABLE, TARGET :: dx(:),f(:,:)
    INTEGER, POINTER :: Perm(:)
    INTEGER :: dofs, i, j, nsize, ControlNode, dof0, nControl, iControl,jControl
    REAL(KIND=dp) :: Nrm, val, cand, mincand, Relax
    LOGICAL :: GotF, Found, UseLoads, ExtremumMode, DiagControl    
    REAL(KIND=dp), ALLOCATABLE :: cAmp(:), cTarget(:), cVal(:), dc(:), cSens(:,:)
    INTEGER, ALLOCATABLE :: cDof(:)
    
    CHARACTER(:), ALLOCATABLE :: str
    CHARACTER(*), PARAMETER :: Caller = 'ControlLinearSystem'

    
    SAVE f, cAmp, cTarget, cSens, cVal, dc, cDof

    IF( ParEnv % PEs > 1 ) THEN
      CALL Fatal(Caller,'Controlling of source terms implemented only in serial!')
    END IF
    
    Params => Solver % Values
    Mesh => Solver % Mesh 
    A => Solver % Matrix
    Var => Solver % Variable    
    b => A % RHS
    x0 => Var % Values
    dofs = Var % Dofs
    Perm => Var % Perm
    nsize = SIZE(x0)


    nControl = ListGetInteger(Params,'Number of Controls',Found ) 
    IF(.NOT. Found ) nControl = 1


    IF( PreSolve ) THEN
      CALL Info(Caller,'Applying controlled sources',Level=7)     
      ALLOCATE(f(nsize,nControl),cAmp(nControl),cTarget(nControl),cVal(nControl),&
          dc(nControl),cSens(nControl,nControl),cDof(nControl))      
      cAmp = 0.0_dp; cTarget = 0.0_dp; cVal = 0.0_dp; dc = 0.0_dp; cSens = 0.0_dp; cDof = 0
      
      f = 0.0_dp
      DO iControl = 1, nControl
       ! This is inhereted from previous control iterations.
       str = 'Control Amplitude'
       IF(nControl > 1) str = TRIM(str)//' '//I2S(iControl)
       cAmp(iControl) = ListGetCReal( Params, str, Found )

       IF(.NOT. Found ) THEN
         str = 'Initial Control Amplitude'
         IF(nControl > 1) str = TRIM(str)//' '//I2S(iControl)
         cAmp(iControl) = ListGetCReal( Params, str, Found )          
        END IF
      END DO
      
      DO iControl = 1, Ncontrol            
        ! Default name for controlled source term
        str = TRIM(Var % Name)//' Control'
        IF(Ncontrol>1) str = TRIM(str)//' '//I2S(iControl)
              
        ! We need to add the control source here in order to be able to use
        ! standard means for convergence monitoring. 
        CALL Info(Caller,'Computing source term for: '//TRIM(str),Level=7)
        CALL SetNodalSources( CurrentModel,Mesh,str,dofs, Perm, GotF, f(:,iControl) )

        ! The additional source needs to be nullified for Dirichlet conditions
        IF( ALLOCATED( A % ConstrainedDOF ) ) THEN
          WHERE( A % ConstrainedDOF ) f(:,iControl) = 0.0_dp
        END IF

        IF(InfoActive(10)) THEN
          DO i=1,dofs
            PRINT *,'ranges b:',i,MINVAL(b(i::dofs)),MAXVAL(b(i::dofs)),SUM(b(i::dofs))
            PRINT *,'ranges f:',i,MINVAL(f(i::dofs,iControl)),&
                MAXVAL(f(i::dofs,iControl)),SUM(f(i::dofs,iControl))
          END DO
        END  IF
       
        IF( ABS(cAmp(iControl)) > 1.0e-20 ) THEN
          b(1:nsize) = b(1:nsize) + cAmp(iControl) * f(1:nsize,iControl)
        END IF
      END DO
    END IF

      
    IF(.NOT. PreSolve ) THEN
      CALL Info(Caller,'Dertermining source term amplitude',Level=7)     
      
      CALL ListPushNameSpace('control:')
      CALL ListAddLogical( Params,'control: Skip Compute Nonlinear Change',.TRUE.)
      CALL ListAddLogical( Params,'control: Skip Advance Nonlinear iter',.TRUE.)

      ALLOCATE(dx(nsize))      
      UseLoads = ListGetLogical( Params,'Control Use Loads', Found )
      IF(UseLoads) THEN        
        ALLOCATE(r0(nsize),dr(nsize))      
      END IF

      DiagControl = ListGetLogical( Params,'Control Diagonal', Found )
      
      dof0 = 1
      IF( dofs > 1) THEN
        dof0 = ListGetInteger( Params,'Control Target Component',UnfoundFatal=.TRUE.)
      END IF
            
      ! Get the target values for control
      DO iControl = 1, Ncontrol            
        str = 'Control Target Value'
        IF(nControl > 1) str = TRIM(str)//' '//I2S(iControl)        
        val = ListGetCReal( Params,str,UnfoundFatal=.TRUE.)
        cTarget(iControl) = val

        i = GetControlNode(Mesh,Perm,Params,iControl) 

        IF(i>0) THEN
          i = dofs*(Perm(i)-1)+dof0
          cDof(iControl) = i
        END IF
      END DO
      
      ! The possibility to use control for extremum temperature is here included.
      ExtremumMode = .FALSE.
      IF( ANY(cDof==0) ) THEN
        IF( Ncontrol == 1 ) THEN
          ExtremumMode = .TRUE.
        ELSE
          CALL Fatal(Caller,'Extremum control cannot be used with '//I2S(Ncontrol)//' controls!')
        END IF
      END IF
      
      
      DO iControl = 1, Ncontrol            
        ! We already know the sources, now compute their affect
        dx = 0.0_dp
        CALL SolveSystem(A,ParMatrix,f(:,iControl),dx,Nrm,dofs,Solver)
        
        ! We use either solution or reaction force for (y0,dy) so that we can
        ! generalize the control procedures for both. 
        IF( UseLoads ) THEN
          ! Nodal loads with the base case
          IF(iControl==1) THEN
            CALL CalculateLoads( Solver, A, x0, dofs, .TRUE., NodalValues = r0 ) 
          END IF

          ! We we use loads then compute the effect of the controlled source to the
          ! reaction force. Hence some hassle with the temporal pointers.          
          BulkRhsSave => A % BulkRhs
          A % BulkRhs => f(:,iControl)
          CALL CalculateLoads( Solver, A, dx, dofs, .TRUE., NodalValues = dr ) 
          A % BulkRhs => BulkRhsSave
          y0 => r0
          dy => dr
        ELSE
          y0 => x0
          dy => dx
        END IF
       
        val = cTarget(iControl) 
        ControlNode = cDof(iControl)

        IF( ExtremumMode ) THEN
          ! We basically do tuning here already but for the sake of uniformity lets just
          ! register the sensitivity and current value. 
          mincand = HUGE(mincand)
          DO i=1,nsize
            j = dofs*(i-1)+dof0          
            IF(ABS(dy(j)) < TINY(dy(j))) CYCLE
            cand = (val-y0(j))/dy(j)
            IF( ABS(cand) < ABS(mincand) ) THEN
              mincand = cand
              cSens(1,1) = dy(j)
              cVal(iControl) = y0(j)
            END IF
          END DO
          CALL Info(Caller,'Extremum value is easiest found in dof: '//I2S(ControlNode),Level=7)    
        ELSE                       
          DO jControl=1,nControl           
            IF(DiagControl .AND. jControl /= iControl) CYCLE
            cSens(jControl,iControl) = dy(cDof(jControl))
          END DO
          cVal(iControl) = y0(cDof(iControl))
        END IF
      END DO

      IF( InfoActive(20) ) THEN
        PRINT *,'cVal:',cVal
        PRINT *,'cTarget:',cTarget
        
        DO i=1,NControl
          PRINT *,'Sens',i,':',cSens(i,:)
        END DO
      END IF
                  
      ! Here we solve the control equation without any assumption of diagonal dominance etc. 
      dc = cTarget - cVal      
      CALL LuSolve(nControl,cSens,dc)

      Relax = ListGetCReal( Params,'Control Relaxation Factor', Found ) 
      IF( Found ) dc = Relax * dc
      

      DO iControl = 1, Ncontrol                    
        str = 'Control Amplitude'
        IF(nControl > 1) str = TRIM(str)//' '//I2S(iControl)        
        cAmp(iControl) = ListGetCReal( Params,str,Found)
        
        cAmp(iControl) = cAmp(iControl) + dc(iControl)
        CALL ListAddConstReal( Params, str, cAmp(iControl) )
        
        ! Apply control, this always to the solution - not to load
        x0(1:nsize) = x0(1:nsize) + dc(iControl) * dx(1:nsize)

        WRITE(Message,'(A,ES15.6)') 'Applied '//TRIM(str)//': ',cAmp(iControl)      
        CALL Info(Caller,Message,Level=5)
      END DO
        
      CALL ListPopNamespace()
      
      DEALLOCATE(f,dx,cAmp,cTarget,cVal,dc,cSens,cDof)
      IF(UseLoads) DEALLOCATE(dr,r0)      
    END IF

    CALL Info(Caller,'All done for now',Level=15)
    
  END SUBROUTINE ControlLinearSystem



!------------------------------------------------------------------------------
!> Given the operation point and an additional r.h.s. source vector find the
!> amplitude for the latter one such that the control problem is resolved.
!> We can request a field value at given point, for example. This tries to
!> mimic some ideas of the "Smart Heater Control" of "HeatSolver" available
!> long ago. This would hopefully be applicable to wider set of modules. 
!------------------------------------------------------------------------------
  SUBROUTINE ControlNonlinearSystem(Solver,PreSolve)
    TYPE(Solver_t), POINTER :: Solver
    LOGICAL :: PreSolve

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: A    
    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), POINTER :: x0(:),b(:),dr(:),r0(:),dy(:),y0(:),prevvalues(:),x(:),dx(:,:)
    INTEGER, POINTER :: Perm(:)
    INTEGER :: dofs, i, j, nsize, ControlNode, dof0, nControl, iControl=0,jControl
    REAL(KIND=dp) :: Nrm, val, cand, mincand, Relax, Eps
    LOGICAL :: GotF, Found, UseLoads, ExtremumMode, DiagControl, Multiply    
    REAL(KIND=dp), ALLOCATABLE :: cAmp(:), cTarget(:), cVal(:), dc(:), cSens(:,:)
    INTEGER, ALLOCATABLE :: cDof(:)
    TYPE(Model_t), POINTER :: Model    
    CHARACTER(:), ALLOCATABLE :: str
    CHARACTER(*), PARAMETER :: Caller = 'ControlNonlinearSystem'
    
    SAVE cAmp, cTarget, cSens, cVal, dc, cDof, iControl, &
        UseLoads, DiagControl, ExtremumMode, Eps, &
        dy, dx, x0, r0, dr, prevvalues 

    IF( ParEnv % PEs > 1 ) THEN
      CALL Fatal(Caller,'Controlling of source terms implemented only in serial!')
    END IF

    Model => CurrentModel
    Params => Solver % Values
    Mesh => Solver % Mesh 
    A => Solver % Matrix
    Var => Solver % Variable    
    b => A % RHS
    x => Var % Values
    dofs = Var % Dofs
    Perm => Var % Perm
    nsize = SIZE(x)


    nControl = ListGetInteger(Params,'Number of Controls',Found ) 
    IF(.NOT. Found ) nControl = 1

    Multiply = .TRUE.
    
    IF( PreSolve ) THEN
      IF( iControl == 0 ) THEN
        CALL Info(Caller,'Applying controlled sources',Level=7)     
        nsize = SIZE(x)
        ALLOCATE(x0(nsize),dx(nsize,nControl),prevvalues(nsize),&
            cAmp(nControl),cTarget(nControl),cVal(nControl),&
            dc(nControl),cSens(nControl,nControl),cDof(nControl))      
        cAmp = 1.0_dp; cTarget = 0.0_dp; cVal = 0.0_dp; dc = 0.0_dp; cSens = 0.0_dp; cDof = 0

        ! Save previous values
        prevvalues = x

        UseLoads = ListGetLogical( Params,'Control Use Loads', Found )
        IF( UseLoads ) THEN
          ALLOCATE(r0(nsize),dr(nsize))
        END IF
                  
        dof0 = 1
        IF( dofs > 1) THEN
          dof0 = ListGetInteger( Params,'Control Target Component',UnfoundFatal=.TRUE.)
        END IF

        ! Get the target values for control
        DO jControl = 1, Ncontrol            
          str = 'Control Target Value'
          IF(nControl > 1) str = TRIM(str)//' '//I2S(jControl)        
          val = ListGetCReal( Params,str,UnfoundFatal=.TRUE.)
          cTarget(jControl) = val
          !i = GetControlNode(jControl)

          i = GetControlNode(Mesh,Perm,Params,jControl) 

          IF(i>0) i = dofs*(Perm(i)-1)+dof0
          cDof(jControl) = i 
        END DO

        ! The possibility to use control for extremum temperature is here included.
        ExtremumMode = .FALSE.
        IF( ANY(cDof==0) ) THEN
          IF( Ncontrol == 1 ) THEN
            ExtremumMode = .TRUE.
          ELSE
            CALL Fatal(Caller,'Extremum control cannot be used with multiple controls!')
          END IF
        END IF

        DiagControl = ListGetLogical( Params,'Control Diagonal', Found ) 

        Eps = ListGetCReal( Params,'Control Epsilon',Found )
        IF(.NOT. Found ) Eps = 0.01_dp
      ELSE IF( iControl == 1 ) THEN        
        CALL ListPushNameSpace('control:')
        CALL ListAddLogical( Params,'control: Skip Compute Nonlinear Change',.TRUE.)
        CALL ListAddLogical( Params,'control: Skip Advance Nonlinear iter',.TRUE.)        
      END IF
    END IF
    

    IF(.NOT. PreSolve) THEN
      IF(iControl == 0 ) THEN
        x0 = Var % Values

        IF(UseLoads) THEN        
          ! Reaction force for the base case
          CALL CalculateLoads( Solver, A, x, dofs, .TRUE., NodalValues = r0 )
        END IF
      ELSE
        ! Remove variation of the parameters
        val = 1.0/(1.0_dp + eps)
        CALL ListSetParameters( Model, iControl, val, multiply, Found )            

        dx(:,iControl) = x - x0
        
        IF(UseLoads) THEN
          ! Reaction force for the variation
          y0 => r0
          CALL CalculateLoads( Solver, A, x, dofs, .TRUE., NodalValues = dr )
          dr = dr - r0
          dy => dr
        ELSE
          y0 => x0
          dy => dx(:,iControl)
        END IF

        val = cTarget(iControl) 
        ControlNode = cDof(iControl)

        IF( ExtremumMode ) THEN
          ! We basically do tuning here already but for the sake of uniformity lets just
          ! register the sensitivity and current value. 
          mincand = HUGE(mincand)
          DO i=1,nsize
            j = dofs*(i-1)+dof0          
            IF(ABS(dy(j)) < TINY(dy(j))) CYCLE
            cand = (val-y0(j))/dy(j)
            IF( ABS(cand) < ABS(mincand) ) THEN
              mincand = cand
              cSens(1,1) = dy(j) / eps
              cVal(iControl) = y0(j)
            END IF
          END DO
          CALL Info(Caller,'Extremum value is easiest found in dof: '//I2S(ControlNode),Level=7)    
        ELSE                       
          DO jControl=1,nControl           
            IF(DiagControl .AND. jControl /= iControl) CYCLE
            cSens(jControl,iControl) = dy(cDof(jControl)) / eps
          END DO
          cVal(iControl) = y0(cDof(iControl))
        END IF
      END IF
                
      IF(iControl == nControl ) THEN
        CALL Info(Caller,'Dertermining source term amplitude',Level=7)     
        
        ! Here we solve the control equation without any assumption of diagonal dominance etc. 
        dc = cTarget - cVal      
        CALL LuSolve(nControl,cSens,dc)

        IF( InfoActive(20) ) THEN
          PRINT *,'cVal:',cVal
          PRINT *,'cTarget:',cTarget          
          DO i=1,NControl
            PRINT *,'Sens',i,':',cSens(i,:)
          END DO
        END IF
        
        Relax = ListGetCReal( Params,'Control Relaxation Factor', Found ) 
        IF( Found ) dc = Relax * dc        
        
        x = x0
        
        DO jControl = 1, Ncontrol                    
          str = 'Control Amplitude'
          IF(nControl > 1) str = TRIM(str)//' '//I2S(iControl)        
          cAmp(jControl) = ListGetCReal( Params,str,Found)
          IF(.NOT. Found) cAmp(jControl) = 1.0_dp

          val = 1.0_dp + dc(jControl)
          cAmp(jControl) = val * cAmp(jControl) 

          IF( .NOT. multiply ) val = cAmp(jControl)             
          CALL ListSetParameters( Model, jControl, val, multiply, Found )            
            
          CALL ListAddConstReal( Params, str, cAmp(jControl) )
        
          ! Apply control, this always to the solution - not to load
          x = x + dc(jControl) * dx(:,jControl)

          WRITE(Message,'(A,ES15.6)') 'Applied '//TRIM(str)//': ',cAmp(iControl)      
          CALL Info(Caller,Message,Level=5)
        END DO
        
        DEALLOCATE(prevvalues,x0,dx)
        IF(UseLoads) DEALLOCATE(r0,dr)
        DEALLOCATE(cAmp,cTarget,cVal,dc,cSens,cDof)
        
        CALL ListPopNamespace()
        iControl = 0
      ELSE
        iControl = iControl + 1
        ! Add variation from the parameters
        val = (1.0_dp + eps)
        CALL ListSetParameters( Model, iControl, val, multiply, Found )                    

        ! Start from the same base case with the matrix assembly
        x = prevvalues
      END IF

    END IF
      
    CALL Info(Caller,'All done for now',Level=15)
    
  END SUBROUTINE ControlNonlinearSystem

  

!------------------------------------------------------------------------------
  SUBROUTINE SaveLinearSystem( Solver, Ain, LinSysName, OffsetInd )
!------------------------------------------------------------------------------
    TYPE( Solver_t ) :: Solver
    TYPE(Matrix_t), POINTER, OPTIONAL :: Ain
    CHARACTER(LEN=*), OPTIONAL :: LinSysName
    INTEGER, OPTIONAL :: OffsetInd
!------------------------------------------------------------------------------    
    TYPE(Matrix_t), POINTER :: A
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(:), ALLOCATABLE :: dumpfile, dumpprefix
    INTEGER, POINTER :: Perm(:)
    REAL(KIND=dp), POINTER :: Sol(:)
    INTEGER :: i
    LOGICAL :: SaveMass, SaveDamp, SavePerm, SaveSol, Found , Parallel, CNumbering, SkipZeros, SaveSum
    CHARACTER(*), PARAMETER :: Caller = 'SaveLinearSystem'
!------------------------------------------------------------------------------

    CALL Info(Caller,'Saving linear system',Level=6)

    Parallel = ParEnv % PEs > 1

    Params => Solver % Values
    IF(.NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal(Caller,'Parameter list not associated!')
    END IF

    CNumbering = ListGetLogical(Params, 'Linear System Save Continuous Numbering',Found)

    IF( PRESENT(Ain)) THEN
      A => Ain
    ELSE
      A => Solver % Matrix
    END IF

    IF(.NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal(Caller,'Matrix not associated!')
    END IF

    SaveSum = ListGetLogical( Params,'Linear System Save Sum',Found)

    SaveMass = ListGetLogical( Params,'Linear System Save Mass',Found)

    SaveDamp = ListGetLogical( Params,'Linear System Save Damp',Found)   

    SkipZeros = ListGetLogical( Params,'Linear System Save Skip Zeros', Found ) 
    
    IF( PRESENT( LinSysName ) ) THEN
      dumpprefix = TRIM(LinSysName) 
    ELSE
      dumpprefix = ListGetString( Params, 'Linear System Save Prefix', Found)
      IF(.NOT. Found ) dumpprefix = 'linsys'
    END IF
          
    dumpfile = TRIM(dumpprefix)//'_a.dat'
    IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)
    CALL Info(Caller,'Saving matrix to: '//TRIM(dumpfile),Level=5)
    OPEN(1,FILE=dumpfile, STATUS='Unknown')
    CALL PrintMatrix(A,Parallel,Cnumbering,SaveMass=SaveMass,SaveDamp=SaveDamp,SkipZeros=SkipZeros)
    CLOSE(1)

    IF( ASSOCIATED(A % rhs) .OR. SaveSum ) THEN
      dumpfile = TRIM(dumpprefix)//'_b.dat'
      IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)
      CALL Info(Caller,'Saving matrix rhs to: '//TRIM(dumpfile),Level=5)
      OPEN(1,FILE=dumpfile, STATUS='Unknown')
      CALL PrintRHS(A, Parallel, CNumbering, SaveSum )
      CLOSE(1)
    END IF
      
    SavePerm = ListGetLogical( Params,'Linear System Save Perm',Found)
    IF( SavePerm ) THEN
      Perm => Solver % Variable % Perm
      IF( .NOT. ASSOCIATED( Perm ) ) THEN
        CALL Warn(Caller,'Permuation not associated!')
        SavePerm = .FALSE.
      ELSE
        dumpfile = TRIM(dumpprefix)//'_perm.dat'
        IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)
        CALL Info(Caller,'Saving permutation to: '//TRIM(dumpfile),Level=5)
        OPEN(1,FILE=dumpfile, STATUS='Unknown')
        DO i=1,SIZE(Perm)
          WRITE(1,'(I0,A,I0)') i,' ',Perm(i)
        END DO
        CLOSE( 1 ) 
      END IF
    END IF


    SaveSol = ListGetLogical( Params,'Linear System Save Solution',Found)
    IF( SaveSol ) THEN
      Sol => Solver % Variable % Values
      IF( .NOT. ASSOCIATED( Sol ) ) THEN
        CALL Warn(Caller,'Solution not associated!')
        SaveSol = .FALSE.
      ELSE
        dumpfile = TRIM(dumpprefix)//'_sol.dat'
        IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)
        CALL Info(Caller,'Saving solution to: '//TRIM(dumpfile),Level=5)
        OPEN(1,FILE=dumpfile, STATUS='Unknown')
        DO i=1,SIZE(Sol)
          WRITE(1,'(I0,ES15.6)') i,Sol(i)
        END DO
        CLOSE( 1 ) 
      END IF
    END IF
    
    
    dumpfile = TRIM(dumpprefix)//'_sizes.dat'
    IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)
    CALL Info(Caller,'Saving matrix sizes to: '//TRIM(dumpfile),Level=5)
    OPEN(1,FILE=dumpfile, STATUS='Unknown')
    WRITE(1,*) A % NumberOfRows
    WRITE(1,*) SIZE(A % Values)
    i = 0
    IF( SavePerm ) i = SIZE( Perm ) 
    WRITE(1,*) i        
    WRITE(1,*) MINVAL(A % Cols)
    WRITE(1,*) MAXVAL(A % Cols)
    IF(PRESENT(OffsetInd)) WRITE(1,*) OffsetInd
    CLOSE(1)

    IF(Parallel) THEN
      dumpfile = TRIM(dumpprefix)//'_sizes.dat'
      CALL Info(Caller,'Saving matrix sizes to: '//TRIM(dumpfile),Level=6)
      OPEN(1,FILE=dumpfile, STATUS='Unknown')
      WRITE(1,*) ParallelReduction(A % ParMatrix % &
                           SplittedMatrix % InsideMatrix % NumberOfRows)
      WRITE(1,*) ParallelReduction(SIZE(A % Values))
      IF( SavePerm ) WRITE(1,*) ParallelReduction(SIZE( Perm ))
      CLOSE(1)
    END IF
    
    IF( ListGetLogical( Params,'Linear System Save and Stop',Found ) ) THEN
      CALL Info(Caller,'Just saved matrix and stopped!',Level=4)
      STOP EXIT_OK
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
    CHARACTER(:), ALLOCATABLE :: BoundaryName
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:), rhs(:)
    REAL(KIND=dp) :: detJ, val
    LOGICAL :: Stat
    
    
    Mesh => Solver % Mesh
        
    N = Mesh % MaxElementNodes 
    ALLOCATE( Basis(n), dBasisdx(n, 3), FORCE(N), STIFF(N,N), &
        Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        STAT=istat)

    ALLOCATE( rhs(A % NumberOfRows) )
    rhs = 0.0_dp
    
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
      
      CALL UpdateGlobalEquations( A,STIFF,rhs,FORCE,n,1,Perm(Indexes(1:n)) )
    END DO
    
    DEALLOCATE( Basis, dBasisdx, FORCE, STIFF, & 
        Nodes % x, Nodes % y, Nodes % z)
    DEALLOCATE( rhs ) 

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

    ! We don't currently need to do this as we copy the ParallelInfo from primary matrix.
    IF ( Parenv  % PEs > 1 ) THEN
      ! CALL ParallelInitMatrix( Solver, A )
    END IF
    
  END SUBROUTINE MassMatrixAssembly



!------------------------------------------------------------------------------
!> Create diagonal matrix from P (not square) by summing the entries together
!> and multiplying with a constant.
!------------------------------------------------------------------------------
  SUBROUTINE DiagonalMatrixSumming( Solver, P, A )
    
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), POINTER :: P, A
    !------------------------------------------------------------------------------
    INTEGER :: i,j,k,n
    REAL(KIND=dp) :: val, rowsum, minsum, maxsum, sumsum

    CALL Info('DiagonalMatrixSumming','Creating diagonal matrix from absolute rowsums',Level=8)

    IF(.NOT. ASSOCIATED(P) ) CALL Fatal('DiagonalMatrixSumming','Matrix P not associated!')
    IF(.NOT. ASSOCIATED(A) ) CALL Fatal('DiagonalMatrixSumming','Matrix A not associated!')

    
    n = P % NumberOfRows 
    CALL Info('DiagonalMatrixSumming','Number of rows in matrix: '//I2S(n),Level=10)

    A % FORMAT = MATRIX_CRS
    
    A % NumberOfRows = n
    ALLOCATE( A % Cols(n), A % Rows(n+1), A % Values(n) )

    A % Cols = 0
    A % Rows = 0
    A % Values = 0.0_dp
    
    minsum = HUGE(minsum)
    maxsum = 0.0_dp
    sumsum = 0.0_dp
    
    DO i = 1, n
      rowsum = 0.0_dp
      DO j=P % Rows(i),P % Rows(i+1)-1
        k = P % Cols(j)
        val = P % Values(j) 
        rowsum = rowsum + ABS( val )
      END DO

      A % Values(i) = rowsum
      A % Cols(i) = i
      A % Rows(i) = i

      minsum = MIN(minsum, rowsum)
      maxsum = MAX(maxsum, rowsum)
      sumsum = sumsum + rowsum
    END DO
    A % Rows(n+1) = n+1

    PRINT *,'diagonal sums:',minsum,maxsum,sumsum/n
    
  END SUBROUTINE DiagonalMatrixSumming



!------------------------------------------------------------------------------
!> Assemble coupling matrix related to fluid-structure interaction
!------------------------------------------------------------------------------
  SUBROUTINE FsiCouplingAssembly( Solver, FVar, SVar, A_f, A_s, A_fs, A_sf, &
      IsPlate, IsShell, IsNS )
    
    TYPE(Solver_t) :: Solver          ! leading solver
    TYPE(Variable_t), POINTER :: FVar ! fluid variable
    TYPE(Variable_t), POINTER :: SVar ! structure variable
    TYPE(Matrix_t), POINTER :: A_fs, A_sf, A_f, A_s
    LOGICAL :: IsPlate, IsShell, IsNS
   !------------------------------------------------------------------------------
    LOGICAL, POINTER :: ConstrainedF(:), ConstrainedS(:)
    INTEGER, POINTER :: FPerm(:), SPerm(:)
    INTEGER :: FDofs, SDofs
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Indexes(:), pIndexes(:)
    INTEGER :: i,j,ii,jj,k,n,t,istat,pn,ifluid,jstruct,pcomp
    TYPE(Element_t), POINTER :: Element, Parent
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:)
    REAL(KIND=dp), POINTER :: Basis(:)
    REAL(KIND=dp) :: detJ, val, c(3), pc(3), Normal(3), coeff, Omega, Rho, area, fdiag
    LOGICAL :: Stat, IsHarmonic, IsTransient
    INTEGER :: dim,mat_id,tcount
    LOGICAL :: FreeF, FreeS, FreeFim, FreeSim, UseDensity, Found
    LOGICAL :: DoMass, DoDamp
    LOGICAL, ALLOCATABLE :: NodeDone(:)
    REAL(KIND=dp) :: MultSF, MultFS, dt
    REAL(KIND=dp), POINTER :: A_fs_values(:), A_sf_values(:)
    CHARACTER(*), PARAMETER :: Caller = 'FsiCouplingAssembly'
    TYPE(Variable_t), POINTER :: dtVar
    
    
    CALL Info(Caller,'Creating coupling matrix for FSI',Level=6)

    
    IF( A_fs % FORMAT /= MATRIX_LIST ) THEN
      A_fs % Values = 0.0_dp
      A_sf % Values = 0.0_dp      
    END IF
      
    
    Mesh => Solver % Mesh
    FPerm => FVar % Perm
    SPerm => SVar % Perm
    
    fdofs = FVar % Dofs
    sdofs = SVar % Dofs


    n = COUNT( FPerm>0 .AND. SPerm>0 ) 
    IF( n == 0 ) THEN
      CALL List_toCRSMatrix(A_fs)
      CALL List_toCRSMatrix(A_sf)
      CALL Info(Caller,'No shared nodes between fluid and structure! Nothing to do!',Level=6)
      RETURN
    END IF

    
    IF( IsPlate ) THEN
      CALL Info(Caller,'Assuming structure to be plate',Level=8)
    ELSE IF( IsShell ) THEN
      CALL Info(Caller,'Assuming structure to be shell',Level=8)
    ELSE
      CALL Info(Caller,'Assuming structure to be solid',Level=8)      
    END IF
      
    IF( IsNS ) THEN
      CALL Info(Caller,'Assuming fluid to have velocities',Level=8)
    ELSE
      CALL Info(Caller,'Assuming fluid to have pressure',Level=8)
    END IF
      

    UseDensity = .FALSE.
    DO i=1,CurrentModel % NumberOfSolvers
      PSolver => CurrentModel % Solvers(i)      
      IF( ASSOCIATED( PSolver % Variable, FVar ) ) THEN
        UseDensity = ListGetLogical( PSolver % Values,'Use Density',Found ) 
        EXIT
      END IF
    END DO
    IF( UseDensity ) THEN
      CALL Info(Caller,'The Helmholtz equation is multiplied by density',Level=10)
    END IF
    
    DoMass = .FALSE.
    DoDamp = .FALSE.
    
    ConstrainedF => A_f % ConstrainedDof
    ConstrainedS => A_s % ConstrainedDof
    
    ! Here we assume harmonic coupling if there are more then 3 structure dofs
    dim = 3
    IsHarmonic = .FALSE.
    IF( IsPlate ) THEN
      IF( sdofs == 6 ) THEN
        IsHarmonic = .TRUE.
      ELSE IF( sdofs /= 3 ) THEN
        CALL Fatal(Caller,'Invalid number of dofs in plate solver: '//I2S(sdofs))
      END IF
    ELSE IF( IsShell ) THEN
      IF( sdofs == 12 ) THEN
        IsHarmonic = .TRUE.
      ELSE IF( sdofs /= 6 ) THEN
        CALL Fatal(Caller,'Invalid number of dofs in shell solver: '//I2S(sdofs))
      END IF
    ELSE
      IF( sdofs == 4 .OR. sdofs == 6 ) THEN
        IsHarmonic = .TRUE.
      ELSE IF( sdofs /= 2 .AND. sdofs /= 3 ) THEN
        CALL Fatal(Caller,'Invalid number of dofs in elasticity solver: '//I2S(sdofs))
      END IF
      IF( sdofs == 4 .OR. sdofs == 2 ) dim = 2
    END IF

    ! The elasticity solver defines whether the system is real or harmonic
    IF( IsHarmonic ) THEN
      CALL Info(Caller,'Assuming harmonic coupling matrix',Level=10)
      IsTransient = .FALSE.
    ELSE
      CALL Info(Caller,'Assuming real valued coupling matrix',Level=10)
      IsTransient = ( ListGetString( CurrentModel % Simulation,&
          'Simulation Type' ) == 'transient' ) 
      IF( IsTransient ) THEN
        CALL Info(Caller,'Assuming transient coupling matrix',Level=10)
      ELSE
        CALL Info(Caller,'Assuming steady-state coupling matrix',Level=10)       
      END IF      
    END IF

    
    ! The fluid system must be consistent with elasticity system
    IF( IsNS ) THEN
      IF( IsHarmonic ) THEN
        IF( fdofs /= 2*(dim+2) .AND. fdofs /= 2*(dim+1) ) THEN
          CALL Fatal(Caller,&
              'Inconsistent number of harmonic dofs in NS solver: '//I2S(fdofs))
        END IF
        ! pressure component
        pcomp = fdofs / 2
      ELSE
        IF( fdofs /= (dim+2) .AND. fdofs /= (dim+1) ) THEN
          CALL Fatal(Caller,&
              'Inconsistent number of real dofs in NS solver: '//I2S(fdofs))
        END IF
        pcomp = fdofs
      END IF
      ALLOCATE( NodeDone(MAXVAL(FPerm)) )
      NodeDone = .FALSE.
    ELSE
      IF( IsHarmonic ) THEN
        IF( fdofs /= 2 ) CALL Fatal(Caller,&
            'Inconsistent number of harmonic dofs in pressure solver: '//I2S(fdofs))
      ELSE
        IF( fdofs /= 1 ) CALL Fatal(Caller,&
            'Inconsistent number of real dofs in pressure solver: '//I2S(fdofs))
      END IF
      pcomp = 1
    END IF

    dt = 0.0_dp
    Omega = 0.0_dp
    
    IF( IsHarmonic ) THEN
      Omega = 2 * PI * ListGetCReal( CurrentModel % Simulation,'Frequency',Stat ) 
      IF( .NOT. Stat) THEN
        CALL Fatal(Caller,'Frequency in Simulation list not found!')
      END IF
    ELSE IF( IsTransient ) THEN
      dtVar => VariableGet( Solver % Mesh % Variables,'timestep size')
      dt = dtVar % Values(1)
    END IF
    
    i = SIZE( FVar % Values ) 
    j = SIZE( SVar % Values ) 
    
    CALL Info(Caller,'Fluid dofs '//I2S(i)//&
        ' with '//I2S(fdofs)//' components',Level=10)
    CALL Info(Caller,'Structure dofs '//I2S(j)//&
        ' with '//I2S(sdofs)//' components',Level=10)   
    CALL Info(Caller,'Assuming '//I2S(dim)//&
        ' active dimensions',Level=10)   

    ! Add the lasrgest entry that allocates the whole list matrix structure
    CALL AddToMatrixElement(A_fs,i,j,0.0_dp)
    CALL AddToMatrixElement(A_sf,j,i,0.0_dp)
    
    N = Mesh % MaxElementNodes 
    ALLOCATE( Basis(n), MASS(N,N), Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        STAT=istat)

    tcount = 0
    area = 0.0_dp
    
    MultFS = ListGetCReal( Solver % Values,'FS multiplier',Found)
    IF( .NOT. Found ) MultFS = 1.0_dp

    MultSF = ListGetCReal( Solver % Values,'SF multiplier',Found)
    IF( .NOT. Found ) MultSF = 1.0_dp
    
    FreeS = .TRUE.; FreeSim = .TRUE.
    FreeF = .TRUE.; FreeFim = .TRUE.    
    
    
100 DO t=Mesh % NumberOfBulkElements+1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      
      Element => Mesh % Elements(t)
      n = Element % TYPE % NumberOfNodes
      Indexes => Element % NodeIndexes
      
      IF( ANY( FPerm(Indexes) == 0 ) ) CYCLE
      IF( ANY( SPerm(Indexes) == 0 ) ) CYCLE
      IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE
      
      Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
      Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
      Nodes % z(1:n) = Mesh % Nodes % z(Indexes)
      
      Normal = NormalVector( Element, Nodes )

    
      ! The following is done in order to check that the normal points to the fluid      
      Parent => Element % BoundaryInfo % Left
      IF( ASSOCIATED( Parent ) ) THEN
        IF( ANY( FPerm(Parent % NodeIndexes) == 0 ) ) Parent => NULL()
      END IF

      IF(.NOT. ASSOCIATED( Parent ) ) THEN
        Parent => Element % BoundaryInfo % Right
        IF( ASSOCIATED( Parent ) ) THEN
          IF( ANY( FPerm(Parent % NodeIndexes) == 0 ) ) Parent => NULL()
        END IF
      END IF
                
      ! Could not find a proper fluid element to define the normal 
      IF(.NOT. ASSOCIATED( Parent ) ) CYCLE

      tcount = tcount + 1

      
      pn = Parent % TYPE % NumberOfNodes
      pIndexes => Parent % NodeIndexes
      
      c(1) =  SUM( Nodes % x(1:n) ) / n
      c(2) =  SUM( Nodes % y(1:n) ) / n
      c(3) =  SUM( Nodes % z(1:n) ) / n
      
      pc(1) =  SUM( Mesh % Nodes % x(pIndexes) ) / pn
      pc(2) =  SUM( Mesh % Nodes % y(pIndexes) ) / pn
      pc(3) =  SUM( Mesh % Nodes % z(pIndexes) ) / pn
      
      ! The normal vector has negative projection to the vector drawn from center of
      ! boundary element to the center of bulk element. 
      IF( SUM( (pc-c) * Normal ) < 0.0_dp ) THEN
        Normal = -Normal
      END IF
      
      MASS(1:n,1:n) = 0.0_dp
      
      mat_id = ListGetInteger( CurrentModel % Bodies(Parent % BodyId) % Values,'Material' )
      rho = ListGetConstReal( CurrentModel % Materials(mat_id) % Values,'Density',Stat)
      IF(.NOT. Stat) rho = ListGetConstReal( CurrentModel % Materials(mat_id) % Values, &
          'Equilibrium Density',Stat)

      IF( .NOT. Stat) THEN
        CALL Fatal(Caller,'Fluid density not found in material :'//I2S(mat_id))
      END IF
      
      ! The sign depends on the convection of the normal direction
      ! If density is divided out already in the Helmholtz equation the multiplier will
      ! be different. 
      IF( UseDensity ) THEN
        coeff = omega**2
      ELSE
        coeff = rho * omega**2
      END IF
      
      ! Numerical integration:
      !----------------------
      IP = GaussPoints( Element )
      
      DO k=1,IP % n        
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
            IP % W(k),  detJ, Basis )
        
        ! The mass matrix of the boundary element
        !----------------------------------------
        val = IP % s(k) * detJ
        DO i=1,n
          DO j=1,n
            MASS(i,j) = MASS(i,j) + val * Basis(i) * Basis(j)
          END DO
        END DO
        area = area + val
      END DO

      ! A: fs
      ! Effect of structure on fluid           
      IF( IsNs ) THEN
        ! For the N-S equation the condition applies directly on the velocity components
        
        DO i=1,n
          ii = Indexes(i)
          j = i
          jj = Indexes(j) ! one-to-one mapping


          IF( NodeDone( Fperm(ii) ) ) CYCLE
          NodeDone( FPerm(ii) ) = .TRUE.
          
          
          DO k=1,dim
            
            ! The velocity component of the fluid
            IF( IsHarmonic ) THEN
              ifluid = fdofs*(FPerm(ii)-1)+2*(k-1)+1
              !IF( ASSOCIATED( ConstrainedF ) ) THEN
              !  FreeF = .NOT. ConstrainedF(ifluid)
              !  FreeFim = .NOT. ConstrainedF(ifluid+1)            
              !END IF
            ELSE
              ifluid = fdofs*(FPerm(ii)-1)+k
              !IF( ASSOCIATED( ConstrainedF ) ) THEN
              !  FreeF = .NOT. ConstrainedF(ifluid)
              !END IF
            END IF
                          
            ! Shell and 3D elasticity are both treated with the same routine
            IF( .NOT. IsPlate ) THEN

              val = 1.0_dp
              IF( IsHarmonic ) THEN
                val = omega
                jstruct = sdofs*(SPerm(jj)-1)+2*(k-1)+1  
              ELSE IF( IsTransient ) THEN
                val = 1.0/dt
                jstruct = sdofs*(SPerm(jj)-1)+k 
                CALL Fatal(Caller,'NS coupling only done for harmonic elasticity!')               
              END IF

                
            ELSE ! If IsPlate
              IF( IsHarmonic ) THEN              
                val = omega * Normal(k)
                
                ! By default the plate should be oriented so that normal points to z
                ! If there is a plate then fluid is always 3D
                IF( Normal(3) < 0 ) val = -val
                jstruct = sdofs*(SPerm(jj)-1)+1
              ELSE
                CALL Fatal(Caller,'NS coupling only done for harmonic plates!')               
              END IF
            END IF

            IF( IsHarmonic ) THEN                                             
              ! Structure load on the fluid: v = i*omega*u
              fdiag = A_f % Values( A_f % diag(ifluid) )
              IF( FreeF ) THEN
                CALL AddToMatrixElement(A_fs,ifluid,jstruct+1,MultFS*val*fdiag)     ! Re 
              ELSE
                CALL AddToMatrixElement(A_fs,ifluid,jstruct+1,0.0_dp)
              END IF
              
              fdiag = A_f % Values( A_f % diag(ifluid+1) )
              IF( FreeFim ) THEN
                CALL AddToMatrixElement(A_fs,ifluid+1,jstruct,-MultFS*val*fdiag)      ! Im
              ELSE                
                CALL AddToMatrixElement(A_fs,ifluid+1,jstruct,0.0_dp )
              END IF

              ! These must be created for completeness because the matrix topology of complex
              ! matrices must be the same for all components.
              CALL AddToMatrixElement(A_fs,ifluid,jstruct,0.0_dp)
              CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,0.0_dp)
            ELSE IF( IsTransient ) THEN
              ! Structure load on the fluid: v = (u-u_prev)/dt              
              fdiag = A_f % Values( A_f % diag(ifluid) )
              IF( FreeF ) THEN
                CALL AddToMatrixElement(A_fs,ifluid,jstruct+1,MultFS*val*fdiag) 
                A_f % rhs(ifluid) = A_f % rhs(ifluid) + MultFS*val*fdiag * &
                    SVar % PrevValues(jstruct,1)
              ELSE
                CALL AddToMatrixElement(A_fs,ifluid,jstruct,0.0_dp)
              END IF                            
              CALL Fatal(Caller,'NS coupling only done for harmonic system!')
            END IF
          END DO
        END DO

      ELSE ! .NOT. IsNS
        ! For pressure equations (Helmholtz) the structure applies a Neumann condition
        
        DO i=1,n
          ii = Indexes(i)

          ! The pressure component of the fluid
          IF( IsHarmonic ) THEN
            ifluid = fdofs*(FPerm(ii)-1)+2*(pcomp-1)+1
            IF( ASSOCIATED( ConstrainedF ) ) THEN
              FreeF = .NOT. ConstrainedF(ifluid)
              FreeFim = .NOT. ConstrainedF(ifluid+1)            
            END IF
          ELSE
            ifluid = fdofs*(FPerm(ii)-1)+pcomp
            IF( ASSOCIATED( ConstrainedF ) ) THEN
              FreeF = .NOT. ConstrainedF(ifluid)
            END IF
          END IF


          DO j=1,n
            jj = Indexes(j)

            ! Shell and 3D elasticity are both treated with the same routine
            IF( .NOT. IsPlate ) THEN

              DO k=1,dim

                val = MASS(i,j) * Normal(k)

                IF( IsHarmonic ) THEN
                  jstruct = sdofs*(SPerm(jj)-1)+2*(k-1)+1  

                  ! Structure load on the fluid: This assembles
                  !
                  !    -1/rho <dp/dn,v> = -omega^2 <u.n,v> = omega^2 <u.m,v> 
                  !
                  ! with the normal vectors satisfying m = -n. Note that the density (rho) 
                  ! must be defined for Helmholtz solver to make it assemble a system
                  ! consistent with the boundary integral -1/rho <dp/dn,v>.
                  IF( FreeF ) THEN
                    CALL AddToMatrixElement(A_fs,ifluid,jstruct,MultFS*val*coeff)     ! Re 
                  ELSE
                    CALL AddToMatrixElement(A_fs,ifluid,jstruct,0.0_dp)
                  END IF

                  IF( FreeFim ) THEN
                    CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,MultFS*val*coeff) ! Im
                  ELSE                
                    CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,0.0_dp )
                  END IF

                  ! These must be created for completeness because the matrix topology of complex
                  ! matrices must be the same for all components.
                  CALL AddToMatrixElement(A_fs,ifluid,jstruct+1,0.0_dp)     
                  CALL AddToMatrixElement(A_fs,ifluid+1,jstruct,0.0_dp)
                ELSE
                  jstruct = sdofs*(SPerm(jj)-1)+k

                  ! Structure load on the fluid: dp/dn = -u. (This seems strange???)                  
                  IF( FreeF ) THEN
                    IF( IsTransient ) THEN                    
                      CALL AddToMatrixElement(A_fs,ifluid,jstruct,-MultFS*val)
                    ELSE
                      IF( DoMass ) THEN
                        CALL AddToMatrixElement(A_fs,ifluid,jstruct,-MultFS*val)
                        !PRINT *,'fs:',ifluid,jstruct,MultFS,val
                      ELSE
                        ! Ensure correct matrix structure
                        CALL AddToMatrixElement(A_fs,ifluid,jstruct,-MultFS*val*0.001)
                      END IF
                    END IF
                  ELSE
                    CALL AddToMatrixElement(A_fs,ifluid,jstruct,0.0_dp)                  
                  END IF
                END IF
              END DO

            ELSE ! If IsPlate

              val = MASS(i,j) 

              ! By default the plate should be oriented so that normal points to z
              ! If there is a plate then fluid is always 3D
              IF( Normal(3) < 0 ) val = -val

              IF( IsHarmonic ) THEN
                jstruct = sdofs*(SPerm(jj)-1)+1

                ! Structure load on the fluid: -1/rho dp/dn = -omega^2 u.n = omega^2 u.m
                IF( FreeF ) THEN
                  CALL AddToMatrixElement(A_fs,ifluid,jstruct,MultFS*val*coeff)     ! Re 
                ELSE
                  CALL AddToMatrixElement(A_fs,ifluid,jstruct,0.0_dp)
                END IF

                IF( FreeFim ) THEN
                  CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,MultFS*val*coeff) ! Im
                ELSE                
                  CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,0.0_dp )
                END IF

                ! These must be created for completeness because the matrix topology of complex
                ! matrices must be the same for all components.
                CALL AddToMatrixElement(A_fs,ifluid,jstruct+1,0.0_dp)
                CALL AddToMatrixElement(A_fs,ifluid+1,jstruct,0.0_dp)
              ELSE
                jstruct = sdofs*(SPerm(jj)-1)+1

                ! Structure load on the fluid: dp/dn = -u. (This seems strange???)
                IF( FreeF ) THEN
                  CALL AddToMatrixElement(A_fs,ifluid,jstruct,-MultFS*val)           
                END IF
              END IF

            END IF

          END DO
        END DO
      END IF

      
      IF( DoMass ) CYCLE

      ! A_sf:
      ! Effect of fluid (pressure) on structure.
      ! Each component get the normal component of the pressure as a r.h.s. term.
      ! The plate equation just gets the full load and is treated separately. 
      !----------------------------------------------------------------------------
      DO i=1,n
        ii = Indexes(i)

        ! The pressure component of the fluid
        IF( IsHarmonic ) THEN
          ifluid = fdofs*(FPerm(ii)-1)+2*(pcomp-1)+1
        ELSE
          ifluid = fdofs*(FPerm(ii)-1)+pcomp
        END IF

        DO j=1,n
          jj = Indexes(j)
          
          ! Shell and 3D elasticity are both treated with the same routine
          IF( .NOT. IsPlate ) THEN

            DO k=1,dim
              
              val = MASS(i,j) * Normal(k)
              
              IF( IsHarmonic ) THEN
                jstruct = sdofs*(SPerm(jj)-1)+2*(k-1)+1  
                
                IF( ASSOCIATED( ConstrainedS ) ) THEN
                  FreeS = .NOT. ConstrainedS(jstruct)
                  FreeSim = .NOT. ConstrainedS(jstruct+1)
                END IF

                ! Fluid load on the structure: tau \cdot n = p * n
                IF( FreeS ) THEN
                  CALL AddToMatrixElement(A_sf,jstruct,ifluid,MultSF*val)           ! Re terms coupling
                ELSE
                  CALL AddToMatrixElement(A_sf,jstruct,ifluid,0.0_dp)
                END IF
                
                IF( FreeSim ) THEN
                  CALL AddToMatrixElement(A_sf,jstruct+1,ifluid+1,MultSF*val)       ! Im
                ELSE
                  CALL AddToMatrixElement(A_sf,jstruct+1,ifluid+1,0.0_dp)
                END IF

                ! These must be created for completeness because the matrix topology of complex
                ! matrices must be the same for all components.
                CALL AddToMatrixElement(A_sf,jstruct,ifluid+1,0.0_dp)
                CALL AddToMatrixElement(A_sf,jstruct+1,ifluid,0.0_dp)
              ELSE
                jstruct = sdofs*(SPerm(jj)-1)+k
                
                IF( ASSOCIATED( ConstrainedS ) ) THEN
                  FreeS = .NOT. ConstrainedS(jstruct)
                END IF

                ! Fluid load on the structure: tau \cdot n = p * n
                IF( FreeS ) THEN
                  CALL AddToMatrixElement(A_sf,jstruct,ifluid,MultSF*val)           
                END IF
                
              END IF
            END DO
            
          ELSE ! If IsPlate
            
            val = MASS(i,j) 
            
            ! By default the plate should be oriented so that normal points to z
            ! If there is a plate then fluid is always 3D
            IF( Normal(3) < 0 ) val = -val
            
            IF( IsHarmonic ) THEN
              jstruct = sdofs*(SPerm(jj)-1)+1
              
              IF( ASSOCIATED( ConstrainedS ) ) THEN
                FreeS = .NOT. ConstrainedS(jstruct)
                FreeSim = .NOT. ConstrainedS(jstruct+1)
              END IF
              
              ! Fluid load on the structure: tau \cdot n = p * n
              IF( FreeS ) THEN
                CALL AddToMatrixElement(A_sf,jstruct,ifluid,MultSF*val)           ! Re terms coupling
              ELSE
                CALL AddToMatrixElement(A_sf,jstruct,ifluid,0.0_dp)
              END IF

              IF( FreeSim ) THEN
                CALL AddToMatrixElement(A_sf,jstruct+1,ifluid+1,MultSF*val)       ! Im
              ELSE
                CALL AddToMatrixElement(A_sf,jstruct+1,ifluid+1,0.0_dp)
              END IF

              ! These must be created for completeness because the matrix topology of complex
              ! matrices must be the same for all components.
              CALL AddToMatrixElement(A_sf,jstruct,ifluid+1,0.0_dp)
              CALL AddToMatrixElement(A_sf,jstruct+1,ifluid,0.0_dp)
            ELSE
              jstruct = sdofs*(SPerm(jj)-1)+1

              IF( ASSOCIATED( ConstrainedS ) ) THEN
                FreeS = .NOT. ConstrainedS(jstruct)
              END IF
              
              ! Fluid load on the structure: tau \cdot n = p * n
              IF( FreeS ) THEN
                CALL AddToMatrixElement(A_sf,jstruct,ifluid,MultSF*val)           
              END IF
            END IF
            
          END IF

        END DO
      END DO

    END DO ! Loop over boundary elements
          
    IF( A_fs % FORMAT == MATRIX_LIST ) THEN
      CALL List_toCRSMatrix(A_fs)
      CALL List_toCRSMatrix(A_sf)
    END IF

    ! The list matrix type only includes one field. Hence when we want to also assembly
    ! mass values we create a CRS, allocate the mass values and do another round.
    ! Alternatively we could have two list matrices being worked at the same time. 
    !------------------------------------------------------------------------------------
    IF(.NOT. DoMass ) THEN
      DoMass = .NOT. ( IsHarmonic .OR. IsTransient) .AND. ASSOCIATED( A_s % MassValues ) 
                 
      IF( DoMass ) THEN
        !PRINT *,'Damp:',ASSOCIATED(A_s % DampValues), ASSOCIATED(A_f % DampValues ) 
        !PRINT *,'Mass:',ASSOCIATED(A_s % MassValues), ASSOCIATED(A_f % MassValues ) 
        
        IF( .NOT. ASSOCIATED( A_f % MassValues ) ) THEN
          CALL Warn(Caller,'Both models should have MassValues if one has!')
        END IF

        IF( InfoActive(20) ) THEN
          PRINT *,'f sum:',SUM(A_f % Values), SUM( ABS( A_f % Values ) )
          !PRINT *,'f damp sum:',SUM(A_f % DampValues), SUM( ABS( A_f % DampValues ) )
          PRINT *,'f mass sum:',SUM(A_f % MassValues), SUM( ABS( A_f % MassValues ) )
          PRINT *,'s sum:',SUM(A_s % Values), SUM( ABS( A_s % Values ) )
          !PRINT *,'s damp sum:',SUM(A_s % DampValues), SUM( ABS( A_s % DampValues ) )
          PRINT *,'s mass sum:',SUM(A_s % MassValues), SUM( ABS( A_s % MassValues ) )
        END IF
          
      ALLOCATE( A_fs % MassValues( SIZE( A_fs % Values ) ) )
        A_fs % MassValues = 0.0_dp
        A_fs_values => A_fs % Values
        A_fs % Values => A_fs % MassValues
        ALLOCATE( A_sf % MassValues( SIZE( A_sf % Values ) ) )
        A_sf % MassValues = 0.0_dp           
        A_sf_values => A_sf % Values
        A_sf % Values => A_sf % MassValues
        CALL Info(Caller,'Looping again over the FSI boundary, now for mass values!',Level=10)      
        GOTO 100 
      END IF
    ELSE 
      ! We repointed the values in order to use the AddToMatrixElement subroutine.
      ! Now let's revert to the original vectors. 
      CALL Info(Caller,'Mass values assembled for the coupling matrices!',Level=10)
      A_fs % Values => A_fs_values
      A_sf % Values => A_sf_values     
    END IF

    DEALLOCATE( Basis, MASS, Nodes % x, Nodes % y, Nodes % z)
    
    IF( InfoActive(20) ) THEN        
      PRINT *,'interface area:',area
      PRINT *,'interface fs sum:',SUM(A_fs % Values), SUM( ABS( A_fs % Values ) )
      PRINT *,'interface sf sum:',SUM(A_sf % Values), SUM( ABS( A_sf % Values ) )
      IF( ASSOCIATED( A_fs % MassValues ) ) THEN
        PRINT *,'interface fs mass sum:',SUM(A_fs % MassValues), SUM( ABS( A_fs % MassValues ) )
      END IF
      IF( ASSOCIATED( A_sf % MassValues ) ) THEN
        PRINT *,'interface fs mass sum:',SUM(A_sf % MassValues), SUM( ABS( A_sf % MassValues ) )
      END IF
    END IF
      
    CALL Info(Caller,'Number of elements on interface: '&
        //I2S(tcount),Level=10)    
    CALL Info(Caller,'Number of entries in fluid-structure matrix: '&
        //I2S(SIZE(A_fs % Values)),Level=10)
    CALL Info(Caller,'Number of entries in structure-fluid matrix: '&
        //I2S(SIZE(A_sf % Values)),Level=10)
    
    CALL Info(Caller,'All done',Level=20)

    
  END SUBROUTINE FsiCouplingAssembly





  
  ! The following function is a copy from ShellSolver.F90.
  ! The suffix Int is added for unique naming.
  !-------------------------------------------------------------------------------
  FUNCTION GetElementalDirectorInt(Mesh, Element, &
      ElementNodes, node) RESULT(DirectorValues) 
  !-------------------------------------------------------------------------------    
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    TYPE(Nodes_t), OPTIONAL, INTENT(IN) :: ElementNodes
    INTEGER, OPTIONAL :: node
    REAL(KIND=dp), POINTER :: DirectorValues(:)
    !-------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Visited = .FALSE., UseElementProperty = .FALSE., UseNormalSolver = .FALSE.
    REAL(KIND=dp) :: Normal(3)
    REAL(KIND=dp), POINTER :: NodalNormals(:)
    REAL(KIND=dp), POINTER :: DirectorAtNode(:)
    REAL(KIND=dp), POINTER :: PropertyValues(:)
    INTEGER :: n    

    SAVE Visited, UseElementProperty, NodalNormals, DirectorAtNode, Nodes
    !-------------------------------------------------------------------------------
        
    IF (.NOT. Visited) THEN
      DirectorValues => GetElementPropertyInt('director', Element)
      UseElementProperty = ASSOCIATED( DirectorValues ) 

      IF (.NOT. UseElementProperty) THEN
        n = CurrentModel % MaxElementNodes
        ALLOCATE( NodalNormals(3*n), Nodes % x(n), Nodes % y(n), Nodes % z(n) ) 
      END IF
      ALLOCATE( DirectorAtNode(3) )
      Visited = .TRUE.
    END IF

    IF ( UseElementProperty ) THEN    
      PropertyValues => GetElementPropertyInt('director', Element)
      IF( PRESENT( node ) ) THEN        
        DirectorAtNode(1:3) = PropertyValues(3*(node-1)+1:3*node)
        DirectorValues => DirectorAtNode
      ELSE
        DirectorValues => PropertyValues
      END IF
      
    ELSE
      IF( PRESENT( ElementNodes ) ) THEN
        Normal = NormalVector( Element, ElementNodes, Check = .TRUE. ) 
      ELSE
        n = Element % Type % NumberOfNodes
        Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
        Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
        Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)
        Normal = NormalVector( Element, Nodes, Check = .TRUE. ) 
      END IF

      IF( PRESENT( node ) ) THEN
        !PRINT *,'Normal:',Normal
        DirectorAtNode(1:3) = Normal(1:3)
        DirectorValues => DirectorAtNode
      ELSE              
        n = Element % TYPE % NumberOfNodes
        NodalNormals(1:3*n:3) = Normal(1)
        NodalNormals(2:3*n:3) = Normal(2)
        NodalNormals(3:3*n:3) = Normal(3)      
        DirectorValues => NodalNormals
      END IF
    END IF

  CONTAINS
        
    FUNCTION GetElementPropertyInt( Name, Element ) RESULT(Values)
      CHARACTER(LEN=*) :: Name
      TYPE(Element_t), POINTER :: Element
      REAL(KIND=dp), POINTER :: Values(:)

      TYPE(ElementData_t), POINTER :: p

      Values => NULL()
      p=> Element % PropertyData

      DO WHILE( ASSOCIATED(p) )
        IF ( Name==p % Name ) THEN
          Values => p % Values
          RETURN
        END IF
        p => p % Next
      END DO
    END FUNCTION GetElementPropertyInt
    
  !-------------------------------------------------------------------------------    
  END FUNCTION GetElementalDirectorInt
  !-------------------------------------------------------------------------------


  
!------------------------------------------------------------------------------
!> Assemble coupling matrices related to structure-structure interaction.
!> A possible scenario is that the diagonal blocks are the matrices of the 
!> solvers listed using the keyword "Block Solvers". The (1,1)-block is then
!> tied up with the value of the first entry in the "Block Solvers" array. 
!------------------------------------------------------------------------------
  SUBROUTINE StructureCouplingAssembly(Solver, FVar, SVar, A_f, A_s, A_fs, A_sf, &
      IsSolid, IsPlate, IsShell, IsBeam, DrillingDOFs)
!------------------------------------------------------------------------------    
    TYPE(Solver_t) :: Solver          !< The leading solver defining block structure 
    TYPE(Variable_t), POINTER :: FVar !< "Slave" structure variable
    TYPE(Variable_t), POINTER :: SVar !< "Master" structure variable
    TYPE(Matrix_t), POINTER :: A_f    !< (2,2)-block for the "slave" variable
    TYPE(Matrix_t), POINTER :: A_s    !< (1,1)-block for the "master" variable
    TYPE(Matrix_t), POINTER :: A_fs   !< (2,1)-block for interaction
    TYPE(Matrix_t), POINTER :: A_sf   !< (1,2)-block for interaction
    LOGICAL :: IsSolid, IsPlate, IsShell, IsBeam !< The type of the slave variable
    LOGICAL :: DrillingDOFs           !< Use drilling rotation formulation for shells
   !------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, POINTER :: ConstrainedF(:), ConstrainedS(:)
    LOGICAL :: Found, DoDamp, DoMass, DoLumping 
    INTEGER, POINTER :: FPerm(:), SPerm(:)
    INTEGER :: FDofs, SDofs
    INTEGER :: i,j,k,n,jf,js,kf,ks,nf,ns,dim,ncount
    REAL(KIND=dp) :: vdiag
    LOGICAL :: EnsureTrans
    CHARACTER(*), PARAMETER :: Caller = 'StructureCouplingAssembly'
   !------------------------------------------------------------------------------

    CALL Info(Caller,'Creating coupling matrix for structures',Level=6)
    
    Mesh => Solver % Mesh
    dim = Mesh % MeshDim

    ! S refers to the first and F to the second block (was fluid):
    FPerm => FVar % Perm
    SPerm => SVar % Perm
    
    fdofs = FVar % Dofs
    sdofs = SVar % Dofs

    IF( IsSolid ) CALL Info(Caller,'Assuming coupling with solid solver',Level=8)
    IF( IsBeam )  CALL Info(Caller,'Assuming coupling with beam solver',Level=8)
    IF( IsPlate ) CALL Info(Caller,'Assuming coupling with plate solver',Level=8)
    IF( IsShell ) CALL Info(Caller,'Assuming coupling with shell solver',Level=8)
    
    ConstrainedF => A_f % ConstrainedDof
    ConstrainedS => A_s % ConstrainedDof
                  
    nf = SIZE( FVar % Values ) 
    ns = SIZE( SVar % Values ) 
    
    CALL Info(Caller,'Slave structure dofs '//I2S(nf)//&
        ' with '//I2S(fdofs)//' components',Level=10)
    CALL Info(Caller,'Master structure dofs '//I2S(ns)//&
        ' with '//I2S(sdofs)//' components',Level=10)   
    CALL Info(Caller,'Assuming '//I2S(dim)//&
        ' active mesh dimensions',Level=10)   

    n = COUNT( FPerm>0 .AND. SPerm>0 ) 
    IF( n == 0 ) THEN
      CALL List_toCRSMatrix(A_fs)
      CALL List_toCRSMatrix(A_sf)
      CALL Info(Caller,'No shared nodes between two structures! Nothing to do!',Level=6)
      RETURN
    END IF
    
    IF( A_fs % FORMAT == MATRIX_LIST ) THEN
      ! Add the outmost entry that allocates the whole list matrix structure
      CALL AddToMatrixElement(A_fs,nf,ns,0.0_dp)
      CALL AddToMatrixElement(A_sf,ns,nf,0.0_dp)
    ELSE
      ! If we are revisiting then initialize the CRS matrices to zero
      A_fs % Values = 0.0_dp
      A_sf % Values = 0.0_dp      
    END IF

    DoMass = .FALSE.
    IF( ASSOCIATED( A_f % MassValues ) ) THEN
      IF( ASSOCIATED( A_s % MassValues ) ) THEN
        DoMass = .TRUE.        
      ELSE
        CALL Warn(Caller,'Both models should have MassValues if one has!')
      END IF
    END IF

    DoDamp = ASSOCIATED( A_f % DampValues )
    IF( DoDamp ) THEN
      CALL Warn(Caller,'Damping matrix values at a coupling interface will be dropped!')
    END IF

    DoLumping = ListGetLogical( CurrentModel % Solver % Values,'Block System Mass Lumping',Found )
    IF(.NOT. Found) DoLumping = .TRUE.

    EnsureTrans = ListGetLogical( CurrentModel % Solver % Values,'Block System Topo Symmetric',Found )
    
    ! This is still under development and not used for anything
    ! Probably this will not be needed at all but rather we need the director.
    ! IF( IsShell ) CALL DetermineCouplingNormals()

    !
    ! Create coupling conditions between a solid and a shell:
    !
    IF (IsShell) THEN
      !
      ! In this block, the directional derivative DOFs of the shell model are set to follow
      ! the displacement of the solid. The interaction conditions for the corresponding forces are 
      ! also created.
      !
      BLOCK
        INTEGER, POINTER :: Indexes(:)
        INTEGER, ALLOCATABLE :: NodeHits(:), InterfacePerm(:)
        INTEGER, ALLOCATABLE :: EdgePerm(:),EdgeShellCount(:),EdgeShellTable(:,:),&
            EdgeSolidCount(:),EdgeSolidTable(:,:)
        INTEGER :: MaxEdgeSolidCount, MaxEdgeShellCount, NoFound
        INTEGER :: InterfaceN, EdgeCount, Phase, Interface_Edge_Id
        INTEGER :: p,lf,ls,ii,jj,m,t,l,k1,k2
        INTEGER :: NormalDir
        REAL(KIND=dp), POINTER :: Director(:)
        REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
        REAL(KIND=dp), ALLOCATABLE :: A_f0(:), rhs0(:), Mass0(:)
        REAL(KIND=dp) :: u,v,w,weight,weight0,detJ,val
        REAL(KIND=dp) :: x, y, z 
        REAL(KIND=dp) :: d(3)
        REAL(KIND=dp), ALLOCATABLE :: TestArray(:)
        TYPE(Element_t), POINTER :: Element, ShellElement, Edge
        TYPE(Nodes_t) :: Nodes
        LOGICAL :: Stat
        LOGICAL :: Reject_Imperfect, CheckWeights

        ! 
        ! Tag the nodes which are shared by a solid and a shell
        !
        n = Mesh % NumberOfNodes
        ALLOCATE( InterfacePerm( n ) )

        InterfacePerm = 0
        InterfaceN = 0
        DO i=1,Mesh % NumberOfNodes
          jf = FPerm(i)      
          js = SPerm(i)
          IF( jf == 0 .OR. js == 0 ) CYCLE

          ! also number the interface
          InterfaceN = InterfaceN + 1
          InterfacePerm(i) = InterfaceN
        END DO
        CALL Info(Caller,'Number of nodes shared by the shell and solid: ' &
            //I2S(InterfaceN),Level=10)

        EdgeCount = 0
        IF (InterfaceN > 0) THEN
          !
          ! We need to create mesh edges to simplify many things
          !
          CALL FindMeshEdges( Mesh, FindFaces=.FALSE. ) 
          !
          ! We now assume that an edgewise construction of the coupling interface is made, 
          ! so that the array EdgePerm will be the primary way to check whether we are on
          ! the coupling interface. Therefore, sharing an isolated node is not enough for
          ! activating the coupling.
          !
          ALLOCATE( EdgePerm( Mesh % NumberOfEdges ) )
          EdgePerm = 0
          DO t=1,Mesh % NumberOfEdges
            Edge => Mesh % Edges(t)
            Indexes => Edge % NodeIndexes

            k = COUNT( InterfacePerm(Indexes(1:2)) > 0 )
            IF( k == 2 ) THEN
              EdgeCount = EdgeCount + 1
              EdgePerm(t) = EdgeCount
            END IF
          END DO

          CALL Info(Caller,'Number of edges at interface: '//I2S(EdgeCount),Level=10)
        END IF

        IF( EdgeCount == 0 ) THEN
          CALL List_toCRSMatrix(A_fs)
          CALL List_toCRSMatrix(A_sf)
          CALL Info(Caller,'Coupling matrices are empty! Nothing to do!',Level=6)
          DEALLOCATE(InterfacePerm)
          IF (InterfaceN > 0) DEALLOCATE(EdgePerm)
          RETURN
        END IF
        
        ! Find the shell elements that have an edge on the solid-shell interface.
        ! Count how many times an interface edge is associated with such shell element.

        ALLOCATE( EdgeShellCount(EdgeCount), EdgeSolidCount(EdgeCount) )
        ALLOCATE( NodeHits( Mesh % NumberOfNodes ) )
        NodeHits = 0

        !
        ! If the user does not prevent, we shall pre-examine interface edges to find whether a limitation of 
        ! implementation is met. Creating imperfect constraints can then be avoided by neglecting 
        ! problematic edges
        !
        Reject_Imperfect = ListGetLogical(CurrentModel % Solver % Values, &
            'Reject Imperfect Couplings', Found)
        IF (.NOT. Found) Reject_Imperfect = .TRUE. 
        Reject_Imperfect = Reject_Imperfect .AND. DrillingDOFs

        DO Phase = 0,1         
          EdgeShellCount = 0                  
          
          DO t=1,Mesh % NumberOfBulkElements
            
            Element => Mesh % Elements(t)
            Indexes => Element % NodeIndexes 
            
            n = Element % TYPE % ElementCode
            ! Shell element must be a triangle or quadrilateral
            IF( n > 500 .OR. n < 300 ) CYCLE
            
            ! We must have at least two interface nodes
            k = COUNT( InterfacePerm(Indexes) > 0 )
            IF( k < 2 ) CYCLE
            
            ! We must have shell equation present everywhere 
            IF(ANY( FPerm(Indexes) == 0 ) ) CYCLE 
            
            ! We should not have the shell immersed in solid
            IF(ALL( SPerm(Indexes) /= 0 ) ) CYCLE 
          
            ! Ok, now associate the shell element with the interface edge

            DO i = 1, Element % TYPE % NumberOfEdges
              j = Element % EdgeIndexes(i)

              k = EdgePerm(j)
              IF( k == 0 ) CYCLE
              Edge => Mesh % Edges(j)

              IF (Reject_Imperfect) THEN
                DO ii = 1, Edge % TYPE % NumberOfNodes
                  DO jj = 1, Element % TYPE % NumberOfNodes
                    IF (Indexes(jj) == Edge % NodeIndexes(ii)) EXIT
                  END DO
                  
                  ! To pre-examine, get the director of the shell element
                  Director => GetElementalDirectorInt(Mesh, Element, Node=jj)

                  !
                  ! With the drilling DOFs, the implementation is limited to cases where the director is
                  ! aligned with one of the global coordinate axes. Now, deactivate interface edges for
                  ! which this condition is not satisfied.
                  !
                  NormalDir = 1              
                  DO p=2,3
                    IF (ABS(Director(p)) > ABS(Director(NormalDir))) NormalDir = p
                  END DO
                  IF (1.0_dp - ABS(Director(NormalDir)) > 1.0d-5) THEN
                    !WRITE(Message,'(A,I0,A,F7.4)') 'Director not properly aligned with axis ',&
                    !    NormalDir,': ',Director(NormalDir)
                    !CALL Warn(Caller, Message)
                    CALL Info(Caller, 'No constraint created for the edge '//I2S(j), Level=10)
                    EdgePerm(j) = 0
                    EXIT
                  END IF
                END DO

                IF (EdgePerm(j) == 0) CYCLE
              END IF

              ! Ok, we have an edge on the interface

              EdgeShellCount(k) = EdgeShellCount(k) + 1
              IF( Phase == 1 ) THEN
                EdgeShellTable(k,EdgeShellCount(k)) = t
              END IF
            END DO

          END DO

          IF( Phase == 0 ) THEN
            MaxEdgeShellCount = MAXVAL( EdgeShellCount )
            ALLOCATE( EdgeShellTable(EdgeCount,MaxEdgeShellCount) )
            EdgeShellTable = 0
            Reject_Imperfect = .FALSE.            

            CALL Info(Caller,'Maximum number of owner-shell elements for one edge: '&
                //I2S(MaxEdgeShellCount),Level=10)
            CALL Info(Caller,'Total number of owner-shell elements for all edges: '&
                //I2S(SUM(EdgeShellCount)),Level=10)
          END IF
        END DO

        NoFound = COUNT( EdgeShellCount == 0 )
        IF (NoFound > 0)  CALL Warn(Caller,'The number of deactivated interface edges: '//I2S(NoFound))
        m = COUNT( EdgeShellCount > 1 ) 
        CALL Info(Caller,'Number of edges with several owner-shell elements: '//I2S(m),Level=10)
 
       
        ! Find the solid elements that have edges on the solid-shell interface.
        ! Count how many times an interface edge is associated with such solid element.
        ! Also, count how many times a node is called to construct an interface edge
        ! when looping over the solid elements.
        !
        DO Phase = 0,1
          EdgeSolidCount = 0
          
          DO t=1,Mesh % NumberOfBulkElements 
            Element => Mesh % Elements(t)
            Indexes => Element % NodeIndexes 
            
            n = Element % TYPE % ElementCode
            ! Solid element must be 3D element
            IF( n < 500 ) CYCLE

            ! We must have at least two interface nodes
            k = COUNT( InterfacePerm(Indexes) > 0 )
            IF( k < 2 ) CYCLE

            ! We must have solid equation present everywhere 
            IF(ANY( SPerm(Indexes) == 0 ) ) CYCLE 
            
            ! Ok, now associate the solid element with the mesh edge
            DO i = 1, Element % Type % NumberOfEdges
              j = Element % EdgeIndexes(i)

              ! Is this an active edge?               
              k = EdgePerm(j)
              IF( k == 0 ) CYCLE

              Edge => Mesh % Edges(j)
              EdgeSolidCount(k) = EdgeSolidCount(k) + 1
              
              IF( Phase == 1 ) THEN
                NodeHits(Edge % NodeIndexes) = NodeHits(Edge % NodeIndexes) + 1
                EdgeSolidTable(k,EdgeSolidCount(k)) = t
              END IF                                           
            END DO
          END DO

          IF( Phase == 0 ) THEN
            MaxEdgeSolidCount = MAXVAL( EdgeSolidCount )
            CALL Info(Caller,'Maximum number of owner-solid elements for one edge: '&
                //I2S(MaxEdgeSolidCount),Level=10)
            CALL Info(Caller,'Total number of owner-solid elements for all edges: '&
                //I2S(SUM(EdgeSolidCount)),Level=10)
            ALLOCATE( EdgeSolidTable(EdgeCount,MaxEdgeSolidCount) )
            EdgeSolidTable = 0
          END IF
        END DO

        !
        ! Finally, start to create entries related to the coupling.
        ! "s" refers to solid and "f" to shell.
        !
        n = Mesh % MaxElementNodes 
        ALLOCATE( Basis(n), dBasisdx(n,3), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

        ! Memorize the original values
        ALLOCATE( A_f0( SIZE( A_f % Values ) ) )
        A_f0 = A_f % Values
                
        IF (DrillingDOFs) THEN
          ALLOCATE(rhs0(SIZE(A_f % rhs)))
          rhs0 = A_f % rhs
          IF (DoMass) THEN
            ALLOCATE(Mass0(SIZE(A_f % MassValues)))
            Mass0 = A_f % MassValues
          END IF
        END IF

        ! First, zero the rows related to directional derivative dofs, 
        ! i.e. the components 4,5,6. 
        !
        DO t=1,Mesh % NumberOfEdges
          IF( EdgePerm(t) == 0 ) CYCLE

          Edge => Mesh % Edges(t)
          Indexes => Edge % NodeIndexes
          DO i=1,Edge % Type % NumberOfNodes
            jf = FPerm(Indexes(i))
            DO lf = 4, 6
              kf = fdofs*(jf-1)+lf

              IF( ConstrainedF(kf) ) CYCLE

              DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
                A_f % Values(k) = 0.0_dp
                IF (DoMass) THEN
                  A_f % MassValues(k) = 0.0_dp
                END IF
                IF( DoDamp) THEN
                  A_f % DampValues(k) = 0.0_dp
                END IF
              END DO
              A_f % rhs(kf) = 0.0_dp
            END DO
          END DO
        END DO

        Reject_Imperfect = ListGetLogical(CurrentModel % Solver % Values, &
            'Reject Imperfect Couplings', Found)
        IF (.NOT. Found) Reject_Imperfect = .TRUE.

        CheckWeights = .TRUE.
        IF (CheckWeights) THEN
          ! We shall perform an additional check that the logic works as planned.
          ALLOCATE( TestArray(6*Mesh % NumberOfNodes) )
          TestArray = 0.0_dp          
        END IF

        ! Create weights for the averaging of constraints when the drilling DOFs
        ! formulation is used. For the default shell formulation the counting of
        ! constraints is more straightforward and the array NodeHits will serve
        ! for the calculation of weights.
        !
        IF (DrillingDOFs) THEN
          DO t=1,Mesh % NumberOfEdges
            
            Interface_Edge_Id = EdgePerm(t)
            IF( Interface_Edge_Id == 0 ) CYCLE

            Edge => Mesh % Edges(t)

            DO k1 = 1, EdgeShellCount(Interface_Edge_Id)

              ShellElement => Mesh % Elements( EdgeShellTable(Interface_Edge_Id,k1) )
              
              ! Get the director at the center of the interface edge
              d(:) = 0.0_dp
              DO ii = 1, Edge % Type % NumberOfNodes
                DO jj = 1, ShellElement % TYPE % NumberOfNodes
                  IF (ShellElement % NodeIndexes(jj) == Edge % NodeIndexes(ii)) EXIT
                END DO
                Director => GetElementalDirectorInt(Mesh, ShellElement, Node=jj)
                d(1:3) = d(1:3) + Director(1:3)
              END DO

              d = 1.0_dp/SQRT(SUM(d**2)) * d
              Director(1:3) = d(1:3)

              NormalDir = 1              
              DO i=2,3
                IF (ABS(Director(i)) > ABS(Director(NormalDir))) NormalDir = i
              END DO

              DO k2 = 1, EdgeSolidCount(Interface_Edge_Id)        
                DO ii = 1, 2
                  i = Edge % NodeIndexes(ii)           
                  jf = FPerm(i)      
                  DO lf = 4, 6
                    kf = fdofs*(jf-1)+lf

                    IF( ConstrainedF(kf) ) CYCLE

                    IF ((lf-3) /= NormalDir) THEN
                      IF (CheckWeights) TestArray(kf) = TestArray(kf) + 1.0_dp
                    END IF
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF

        PICK_INTERFACE_EDGES: DO t=1,Mesh % NumberOfEdges

          Interface_Edge_Id = EdgePerm(t)
          IF( Interface_Edge_Id == 0 ) CYCLE         

          Edge => Mesh % Edges(t)
          Indexes => Edge % NodeIndexes          
          
          ! We shall retrieve tying information at the centre of the interface edge.
          ! By construction an edge is identified by two nodes, but let's be generic.
          n = Edge % Type % NumberOfNodes
          x = SUM( Mesh % Nodes % x(Indexes) ) / n
          y = SUM( Mesh % Nodes % y(Indexes) ) / n 
          z = SUM( Mesh % Nodes % z(Indexes) ) / n 

          !
          ! If the interface edge has several owner-shell elements, distinct shell 
          ! elements may have different director definitions and then different constraints 
          ! may arise. The following loop is used to construct an average constraint.
          !
          ! The count of the different approximations of a rotation is 
          ! EdgeShellCount(Interface_Edge_Id) per an edge. To get a single
          ! approximation per an edge, we should assemble a fraction 
          ! 1/EdgeShellCount(Interface_Edge_Id) per an edge to average. We shall
          ! average the rotations at an edge as constants, so that the same
          ! weighting works for the weighting of nodal rotations.
          !
          weight0 = 1.0_dp / EdgeShellCount(Interface_Edge_Id)
          SHELLS_SHARING_INTERFACE_EDGE: DO k1 = 1, EdgeShellCount(Interface_Edge_Id)

            ShellElement => Mesh % Elements( EdgeShellTable(Interface_Edge_Id,k1) )

            ! Get the director at the center of the interface edge
            d(:) = 0.0_dp
            DO ii = 1, n
              DO jj = 1, ShellElement % TYPE % NumberOfNodes
                IF (ShellElement % NodeIndexes(jj) == Edge % NodeIndexes(ii)) EXIT
              END DO
              Director => GetElementalDirectorInt(Mesh, ShellElement, Node=jj)
              d(1:3) = d(1:3) + Director(1:3)
            END DO

            d = 1.0_dp/SQRT(SUM(d**2)) * d
            Director(1:3) = d(1:3)

            IF (DrillingDOFs) THEN
              !
              ! In the case of drilling rotation formulation, the tangential components
              ! trace of the global rotations ROT is related to the directional derivative
              ! of the displacement field u by -Du[d] x d = d x ROT x d. This implementation 
              ! is limited to cases where the director is aligned with one of the global
              ! coordinate axes. Find the closest one and use that. 
              !
              NormalDir = 1              
              DO i=2,3
                IF (ABS(Director(i)) > ABS(Director(NormalDir))) NormalDir = i
              END DO
              IF (.NOT. Reject_Imperfect) THEN
                ! This is not good, but maybe not bad enough to throw the whole analysis away...
                IF (1.0_dp - ABS(Director(NormalDir)) > 1.0d-5) THEN
                  WRITE(Message,'(A,I0,A,F7.4)') 'Director not properly aligned with axis ',&
                      NormalDir,': ',Director(NormalDir)
                  CALL Warn(Caller,Message)
                END IF
              END IF
            END IF
        
            ! Then go through each solid element associated with the interface edge and
            ! create matrix entries defining the interaction conditions for the
            ! directional derivatives and corresponding forces. 
            !
            SOLIDS_SHARING_INTERFACE_EDGE: DO k2 = 1, EdgeSolidCount(Interface_Edge_Id)

              Element => Mesh % Elements( EdgeSolidTable(Interface_Edge_Id,k2) )
              Indexes => Element % NodeIndexes 

              n = Element % TYPE % NumberOfNodes
              Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
              Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
              Nodes % z(1:n) = Mesh % Nodes % z(Indexes)

              ! TO DO: The following call may not work for p-elements!
              CALL GlobalToLocal( u, v, w, x, y, z, Element, Nodes )

              ! Evaluate the basis functions for the solid at the center of the interface edge:
              stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )          

              NODES_ON_EDGE: DO ii = 1, 2
                i = Edge % NodeIndexes(ii)           
                jf = FPerm(i)      

                ! It is not self-evident how we should sum up conditions when several
                ! shell and solid elements are related to a single edge. Ideally, only
                ! one visit at a node would be sufficient to write a constraint. Now
                ! we average by the count of alternative constraints which can be created
                ! at the node, so that we end up to a direct expression for the nodal
                ! rotation.

                IF (.NOT. DrillingDOFs) THEN
                  Weight = 1.0_dp / NodeHits(i) 
                  weight = weight0 * weight
                END IF

                ROTATION_DOFS: DO lf = 4, 6
                  kf = fdofs*(jf-1)+lf

                  IF( ConstrainedF(kf) ) CYCLE

                  DRILLING_OR_NOT: IF (DrillingDOFs) THEN
                    Weight = 1.0_dp/TestArray(kf)
                    IF ((lf-3) /= NormalDir) THEN

                      DO p = 1,n
                        js = SPerm(Indexes(p))

                        IF (NormalDir == 1) THEN
                          SELECT CASE(lf)
                          CASE(5)
                            ks = sdofs*(js-1)+3
                            val = dBasisdx(p,1)
                          CASE(6)
                            ks = sdofs*(js-1)+2
                            val = -dBasisdx(p,1)
                          END SELECT
                        ELSE IF (NormalDir == 2) THEN
                          SELECT CASE(lf)
                          CASE(4)
                            ks = sdofs*(js-1)+3
                            val = -dBasisdx(p,2)
                          CASE(6)
                            ks = sdofs*(js-1)+1
                            val = dBasisdx(p,2)
                          END SELECT
                        ELSE IF (NormalDir == 3) THEN
                          SELECT CASE(lf)
                          CASE(4)
                            ks = sdofs*(js-1)+2
                            val = dBasisdx(p,3)
                          CASE(5)
                            ks = sdofs*(js-1)+1
                            val = -dBasisdx(p,3)
                          END SELECT
                        END IF
                        
                        CALL AddToMatrixElement(A_fs,kf,ks,weight*val)
                        DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
                          CALL AddToMatrixElement(A_sf,ks,A_f % Cols(k),-weight*val*A_f0(k)) 
                        END DO

                        IF( EnsureTrans ) THEN
                          CALL AddToMatrixElement(A_sf,ks,kf,0.0_dp)
                          DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
                            CALL AddToMatrixElement(A_fs,A_f % Cols(k),ks,0.0_dp)
                          END DO
                        END IF                          
                        
                      END DO

                    ELSE
                      !
                      ! Return one row of deleted values to the shell matrix
                      !
                      ! TO DO: Check what happens when several shell elements
                      ! with different directors share a node
                      !
                      DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
                        A_f % Values(k) = A_f0(k)
                        IF (DoMass) A_f % MassValues(k) = Mass0(k)
                      END DO

                      A_f % rhs(kf) = rhs0(kf)

                      ! TO DO: Return also damp values if used

                    END IF

                  ELSE
                    !
                    ! Directional derivative dofs D_{i+3} of the shell equations: 
                    ! We try to enforce the condition D_{i+3}=-<(grad u)d,e_i> 
                    ! where i=1,2,3; i+3=lf, d is director, e_i is unit vector, and 
                    ! u is the displacement field of the solid.
                    !
                    DO p = 1, n
                      js = SPerm(Indexes(p))
                      ks = sdofs*(js-1)+lf-3
                      DO ls = 1, dim
                        val = Director(ls) * dBasisdx(p,ls)

                        !PRINT *,'elem:',Element % ElementIndex, ShellElement % ElementIndex, &
                        !    ls,js,ks,kf,weight,val,Director(ls)

                        CALL AddToMatrixElement(A_fs,kf,ks,weight*val)
                        IF( EnsureTrans ) THEN
                          CALL AddToMatrixElement(A_sf,ks,kf,0.0_dp)
                        END IF
                        
                        ! Here the idea is to distribute the implicit moments of the shell solver
                        ! to forces for the solid solver. So even though the stiffness matrix related to the
                        ! directional derivatives is nullified, the forces are not forgotten.
                        ! This part may be thought of as being based on two (Råback's) conjectures: 
                        ! in the first place the Lagrange variable formulation should bring us to a symmetric 
                        ! coefficient matrix and the values of Lagrange variables can be estimated as nodal 
                        ! reactions obtained by performing a matrix-vector product.
                        !
                        ! Note that no attempt is currently made to transfer external moment
                        ! loads of the shell model to loads of the coupled model. Likewise
                        ! rotational inertia terms of the shell model are not transformed
                        ! to inertia terms of the coupled model. Neglecting the rotational
                        ! inertia might be acceptable in many cases.
                        !
                        ! Note that the minus sign of the entries is correct here:
                        DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
                          CALL AddToMatrixElement(A_sf,ks,A_f % Cols(k),-weight*val*A_f0(k)) 
                        END DO
                        
                        IF( EnsureTrans ) THEN
                          DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
                            CALL AddToMatrixElement(A_fs,A_f % Cols(k),ks,0.0_dp)
                          END DO
                        END IF
                      END DO
                    END DO
                  END IF DRILLING_OR_NOT

                  ! The diagonal entry should sum up to unity!
                  IF (DrillingDOFs) THEN
                    IF ((lf-3) /= NormalDir) THEN
                      CALL AddToMatrixElement(A_f,kf,kf,weight)
                    END IF
                  ELSE
                    CALL AddToMatrixElement(A_f,kf,kf,weight)
                    IF (CheckWeights) TestArray(kf) = TestArray(kf) + weight
                  END IF

                END DO ROTATION_DOFS
              END DO NODES_ON_EDGE
            END DO SOLIDS_SHARING_INTERFACE_EDGE
          END DO SHELLS_SHARING_INTERFACE_EDGE
        END DO PICK_INTERFACE_EDGES

        IF (.NOT. DrillingDOFs .AND. CheckWeights) THEN
          Found = .FALSE.
          DO i=1,SIZE(TestArray)
            IF (ABS(TestArray(i)) > 1.0d-8) THEN
              IF (ABS(TestArray(i)-1.0_dp) > 1.0d-8) THEN
                PRINT *, 'weight (index,value)',i,testarray(i)
                Found = .TRUE.
              END IF
            END IF
          END DO
          IF (Found) THEN
            CALL Fatal(Caller, 'The coupling creates a weight inconsistency')
          ELSE
            CALL Info(Caller, 'Weight check done', Level=10)
          END IF
        END IF
        
        DEALLOCATE( Basis, dBasisdx, Nodes % x, Nodes % y, Nodes % z )
        DEALLOCATE(A_f0, NodeHits, InterfacePerm, EdgePerm, EdgeShellCount, EdgeSolidCount, &
            EdgeShellTable, EdgeSolidTable)
        IF (DrillingDOFs) THEN
          DEALLOCATE(rhs0)
          IF (DoMass) DEALLOCATE(Mass0)
        END IF
        IF (CheckWeights) DEALLOCATE(TestArray)

      END BLOCK
    END IF

    
    ! Note: we may have to recheck this coupling if visiting for 2nd time!
    !
    ! Three DOFs for both shells and solids are the real Cartesian components of
    ! the displacement. Hence we can deal with the common parts of solid-solid and 
    ! solid-shell coupling in same subroutine.
    !
    IF( IsSolid .OR. IsShell ) THEN  
      ncount = 0
      DO i=1,Mesh % NumberOfNodes
        jf = FPerm(i)      
        js = SPerm(i)

        ! We set coupling at nodes that belong to both equations.
        IF( jf == 0 .OR. js == 0 ) CYCLE
        ncount = ncount + 1

        ! For the given node go through all displacement components. 
        DO j = 1, dim
          ! Indices for matrix rows
          kf = fdofs*(jf-1)+j
          ks = sdofs*(js-1)+j

          ! This is the original diagonal entry of the stiffness matrix.
          ! Let's keep it so that Dirichlet conditions are ideally set. 
          vdiag = A_f % Values( A_f % Diag(kf) ) 

          ! Copy the force from rhs from "F" to "S" and zero it
          A_s % rhs(ks) = A_s % rhs(ks) + A_f % rhs(kf)
          A_f % rhs(kf) = 0.0_dp

          ! Copy the force in implicit form from "F" to "SF" coupling matrix, and zero it.
          ! Now the solid equation includes forces of both equations. 
          DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
            IF( .NOT. ConstrainedS(ks) ) THEN        
              CALL AddToMatrixElement(A_sf,ks,A_f % Cols(k), A_f % Values(k) )
              IF( EnsureTrans ) THEN
                CALL AddToMatrixElement(A_fs,A_f % Cols(k),ks,0.0_dp)
              END IF
            END IF
            A_f % Values(k) = 0.0_dp

            ! We zero the mass associated to the Dirichlet conditions since
            ! otherwise the inertia will affect the condition.
            ! We use mass lumping since not all dofs of shell are present in the solid. 
            IF( DoLumping ) THEN
              IF( DoMass ) THEN
                IF( .NOT. ConstrainedS(ks) ) THEN                     
                  A_s % MassValues(A_s % Diag(ks)) = A_s % MassValues(A_s % Diag(ks)) + A_f % MassValues(k)
                END IF
                A_f % MassValues(k) = 0.0_dp
              END IF
              IF( DoDamp) THEN
                A_f % DampValues(k) = 0.0_dp
              END IF
            END IF
          END DO
          
          ! Set Dirichlet Condition to "F" such that it is equal to "S".
          ! Basically we could eliminate displacement condition and do this afterwards
          ! but this is much more economical. 
          A_f % Values( A_f % Diag(kf)) = vdiag          
          CALL AddToMatrixElement(A_fs,kf,ks, -vdiag )

          IF( EnsureTrans ) THEN
            CALL AddToMatrixElement(A_sf,ks,kf,0.0_dp )
          END IF
        END DO
      END DO
    ELSE
      CALL Fatal(Caller,'Coupling type not implemented yet!')
    END IF

    IF( A_fs % FORMAT == MATRIX_LIST ) THEN
      CALL List_toCRSMatrix(A_fs)
      CALL List_toCRSMatrix(A_sf)
    END IF

    CALL Info(Caller,'Number of nodes on interface: '&
        //I2S(ncount),Level=10)    
    CALL Info(Caller,'Number of entries in slave-master coupling matrix: '&
        //I2S(SIZE(A_fs % Values)),Level=10)
    CALL Info(Caller,'Number of entries in master-slave coupling matrix: '&
        //I2S(SIZE(A_sf % Values)),Level=10)

    
    IF(.NOT. DoLumping ) THEN
      ! This is just summary version of the previous for mass & damp values
      ! Now we know that the CRS matrix has been created. 
      BLOCK
        REAL(KIND=dp), POINTER :: TmpValues(:),DerValues(:)
        INTEGER :: der
        
        TmpValues => A_sf % Values

        DO der=1,2
          IF( der == 1 .AND. .NOT. DoDamp ) CYCLE
          IF( der == 2 .AND. .NOT. DoMass ) CYCLE
          
          ncount = 0
          
          IF( der == 1 ) THEN
            CALL Info(Caller,'Creating cross-terms for damping matrix',Level=10)
            IF(.NOT. ASSOCIATED( A_sf % DampValues ) ) THEN
              ALLOCATE( A_sf % DampValues(SIZE(TmpValues) ) )
            END IF
            A_sf % DampValues = 0.0_dp
            ! We set pointer to DampValues so we can use AddToMatrixElement routine
            A_sf % Values => A_sf % DampValues
            DerValues => A_f % DampValues 
          ELSE
            CALL Info(Caller,'Creating cross-terms for mass matrix',Level=10)
            IF(.NOT. ASSOCIATED( A_sf % MassValues ) ) THEN
              ALLOCATE( A_sf % MassValues(SIZE(TmpValues) ) )
            END IF
            A_sf % MassValues = 0.0_dp
            A_sf % Values => A_sf % MassValues
            DerValues => A_f % MassValues 
          END IF

          DO i=1,Mesh % NumberOfNodes
            jf = FPerm(i)
            js = SPerm(i)          
            IF( jf == 0 .OR. js == 0 ) CYCLE
            
            DO j = 1, dim
              kf = fdofs*(jf-1)+j
              ks = sdofs*(js-1)+j                                    
              DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
                IF( .NOT. ConstrainedS(ks) ) THEN        
                  CALL AddToMatrixElement(A_sf,ks,A_f % Cols(k), DerValues(k) )
                  ncount = ncount + 1
                END IF
                DerValues(k) = 0.0_dp              
              END DO
            END DO
          END DO
          CALL Info(Caller,'Number of entries at interface: '//I2S(ncount),Level=10)    
        END DO

        A_sf % Values => TmpValues
        
      END BLOCK
    END IF     
    
    CALL Info(Caller,'Structural coupling matrices created!',Level=20)
    
  CONTAINS


    ! This routine determines normals of the solid at the common nodes with shell solver.
    ! The normals are determined by summing up potential outer normals and thereafter
    ! subtracting projections to the shell normals.
    !------------------------------------------------------------------------------------
    SUBROUTINE DetermineCouplingNormals()
      INTEGER, ALLOCATABLE :: CouplingPerm(:)
      REAL(KIND=dp), ALLOCATABLE, TARGET :: CouplingNormals(:,:)
      REAL(KIND=dp), POINTER :: WallNormal(:)
      REAL(KIND=dp) :: Normal(3), sNormal
      INTEGER :: CouplingNodes, n, t, nbulk, nbound
      TYPE(Element_t), POINTER :: Element, Parent1, Parent2
      TYPE(Nodes_t), SAVE :: Nodes
      LOGICAL :: Solid1,Solid2
      

      ! allocate elemental stuff
      n = Mesh % MaxElementNodes
      IF ( .NOT. ASSOCIATED( Nodes % x ) ) THEN
        ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )
      END IF
     
      ! Generate the permutation for the common nodes
      n = Mesh % NumberOfNodes
      ALLOCATE(CouplingPerm(n))
      WHERE( FVar % Perm(1:n) > 0 .AND. SVar % Perm(1:n) > 0 )
        CouplingPerm = 1
      ELSE WHERE
        CouplingPerm = 0
      END WHERE
      j = 0 
      DO i=1,n
        IF( CouplingPerm(i) > 0 ) THEN
          j = j + 1
          CouplingPerm(i) = j
        END IF
      END DO
      CouplingNodes = j      
      PRINT *,'number of common nodes:',j

      ALLOCATE( CouplingNormals(j,3) )
      CouplingNormals = 0.0_dp
      
      nbulk = Mesh % NumberOfBulkElements
      nbound = Mesh % NumberOfBoundaryElements
      
      ! Sum up all the wall normals associated to coupling nodes together
      DO t=nbulk+1, nbulk+nbound
        Element => Mesh % Elements(t)

        ! If there a node for which we need normal? 
        IF( COUNT( CouplingPerm( Element % NodeIndexes ) > 0 ) < 2 ) CYCLE

        IF( ANY( SVar % Perm( Element % NodeIndexes ) == 0 ) ) CYCLE
        
        ! This needs to be an element where normal can be defined
        !IF( GetElementDim(Element) /= 2 ) CYCLE
        IF( Element % TYPE % ElementCode > 500 ) CYCLE
        IF( Element % TYPE % ElementCode < 300 ) CYCLE
   
        IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE
        
        n = Element % TYPE % NumberOfNodes

        !CALL GetElementNodes(Nodes,Element)
        Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
        Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
        Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)
        
        ! Determine whether parents also are active on the solid
        Solid1 = .FALSE.
        Parent1 => Element % BoundaryInfo % Left
        IF( ASSOCIATED( Parent1 ) ) THEN
          Solid1 = ALL(  SVar % Perm( Parent1 % NodeIndexes ) > 0 )
        END IF
        
        Solid2 = .FALSE.
        Parent2 => Element % BoundaryInfo % Right
        IF( ASSOCIATED( Parent2 ) ) THEN
          Solid2 = ALL(  SVar % Perm( Parent2 % NodeIndexes ) > 0 )
        END IF        

        ! Only consider external walls with just either parent in solid
        IF( .NOT. ( Solid1 .NEQV. Solid2 ) ) CYCLE
        
        ! Check that the normal points outward of the solid
        IF( Solid1 ) THEN
          Normal = NormalVector(Element,Nodes,Parent=Parent1)
        ELSE
          Normal = NormalVector(Element,Nodes,Parent=Parent2)
        END IF
        
        n = Element % TYPE % NumberOfNodes
        DO i=1,n          
          j = CouplingPerm( Element % NodeIndexes(i) )
          IF( j == 0 ) CYCLE
          
          ! Note that we assume that normals are consistent in a way that they can be summed up
          ! and do not cancel each other
          WallNormal => CouplingNormals(j,:)
          WallNormal = WallNormal + Normal 
        END DO                  
      END DO

      ! Remove the shell normal from the wall normal
      DO t=1, nbulk+nbound
        Element => Mesh % Elements(t)

        ! If there a node for which we need normal? 
        IF( COUNT( CouplingPerm( Element % NodeIndexes ) > 0 ) < 2 ) CYCLE

        ! Shell must be active for all nodes
        IF( ANY( FVar % Perm( Element % NodeIndexes ) == 0 ) ) CYCLE

        ! This needs to be an element where shell can be solved
        !IF( GetElementDim(Element) /= 2 ) CYCLE
        IF( Element % TYPE % ElementCode > 500 ) CYCLE
        IF( Element % TYPE % ElementCode < 300 ) CYCLE
        
        n = Element % TYPE % NumberOfNodes

        !CALL GetElementNodes(Nodes,Element)
        Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
        Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
        Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)

        ! Normal vector for shell, no need check the sign
        Normal = NormalVector(Element,Nodes)
 
        DO i=1,n
          j = CouplingPerm( Element % NodeIndexes(i) )
          IF( j == 0 ) CYCLE
          WallNormal => CouplingNormals(j,:)
          WallNormal = WallNormal - SUM( WallNormal * Normal ) * Normal
        END DO
      END DO

      ! Finally normalize the normals such that their length is one
      j = 0
      DO i=1,CouplingNodes
        WallNormal => CouplingNormals(i,:)
        sNormal = SQRT( SUM( WallNormal**2) )
        IF( sNormal > 1.0d-3 ) THEN
          WallNormal = WallNormal / sNormal 
          PRINT *,'WallNormal:',WallNormal
        ELSE
          j = j + 1
        END IF
      END DO
      
      IF( j > 0 ) THEN
        CALL Fatal('DetermineCouplingNormals','Could not define normals count: '//I2S(j))
      END IF
       
      
    END SUBROUTINE DetermineCouplingNormals


  END SUBROUTINE StructureCouplingAssembly

  
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
    CHARACTER(:), ALLOCATABLE :: str, VarName

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
      
      str = 'Linear System Rhs Factor'
      Coeff = ListGetCReal( Params, str, Found )
      IF( Found ) THEN
        Mode = 1
        WRITE( Message,'(A,ES12.3)') 'Multiplying the rhs with ',Coeff
        CALL Info('LinearSystemMultiply',Message, Level=6 )
      ELSE
        Coeff = ListGetCReal( Params, str//' '//I2S(k), Found )
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
          VarName = ListGetString( Params,str//' '//I2S(k),Found )
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
          Coeff = ListGetCReal( Params, str//' '//I2S(k)//I2S(l), Found )
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
            VarName = ListGetString( Params,str//' '//I2S(k)//I2S(l),Found )
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
        Coeff = ListGetCReal( Params, str//' '//I2S(k), Found )
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
          VarName = ListGetString( Params,str//' '//I2S(k),Found )
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
!> Set the diagonal entry to given minimum.
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
    CHARACTER(:), ALLOCATABLE :: str
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

    ! Set the minimum value for each component, only nodel dofs considered
    !---------------------------------------------------------------------
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
        Coeff = ListGetCReal( Params, str//' '//I2S(k), Found )
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

    CALL Info('LinearSystemMinDiagonal',&
        'Number of diagonal values set to minimum: '//I2S(NoSet),Level=5)
    WRITE( Message,'(A,ES12.3)') 'Average abs(diagonal) entry: ',DiagSum / n
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
    REAL(KIND=dp), POINTER :: BulkValues(:),u(:),M_L(:),udot(:), &
        pp(:),pm(:),qp(:),qm(:),corr(:),ku(:),ulow(:)
    REAL(KIND=dp), ALLOCATABLE :: mc_udot(:), fct_u(:)
    REAL(KIND=dp), POINTER CONTIG :: M_C(:), SaveValues(:)
    REAL(KIND=dp) :: rsum, Norm,m_ij,k_ij,du,d_ij,f_ij,c_ij,Ceps,CorrCoeff,&
        rmi,rmj,rpi,rpj,dt
    TYPE(Variable_t), POINTER :: Var, Variables
    LOGICAL :: Found, Symmetric, SaveFields, SkipCorrection
    CHARACTER(:), ALLOCATABLE :: VarName, TmpVarName

    REAL(KIND=dp), POINTER :: mmc(:), mmc_h(:), fct_d(:)
    TYPE(Element_t), POINTER :: Element
    LOGICAL, ALLOCATABLE :: ActiveNodes(:)

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
    Perm => Solver % Variable % Perm

    M_L => A % MassValuesLumped 
    IF (ParEnv % PEs>1) CALL ParallelSumVector(A,M_L)
    
    Var => VariableGet( Variables,'timestep size')
    dt = Var % Values(1) 

    ! low order solution at the start, high order in the end
    u => Solver % Variable % Values
    VarName = GetVarName(Solver % Variable)

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

    TmpVarName = TRIM( VarName )//' fctmm'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    Var % Values = M_L

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
#if 0
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
    CALL ListAddLogical( Params,'fct: Skip Compute Nonlinear Change',.TRUE.)
    CALL ListAddLogical( Params,'fct: Skip Advance Nonlinear iter',.TRUE.)
    SaveValues => A % Values
    A % Values => M_C
    CALL SolveLinearSystem( A, ku, udot, Norm, 1, Solver )
    A % Values => SaveValues
    CALL ListPopNamespace()
#else

  BLOCK
    REAL(KIND=dp), ALLOCATABLE :: TmpRhsVec(:), TmpXVec(:)

    SaveValues => A % Values

    IF (Parenv % PEs>1) THEN
      ALLOCATE(TmpRHSVec(n), TmpXVec(n))
      TmpxVec = 0._dp; tmpRHSVec = 0._dp

      A % Values => BulkValues
      CALL ParallelInitSolve(A,TmpXVec,TmpRhsVec,u)
      CALL ParallelVector(A,TmpRhsvec,u)

      CALL ParallelMatrixVector(A,TmpRhsvec,TmpXVec)

      CALL PartitionVector(A,Ku,TmpXVec)
      DEALLOCATE(TmpRhsVec, TmpXVec)
    ELSE
      DO i=1,n
        rsum = 0._dp
        DO k=Rows(i),Rows(i+1)-1
          j = Cols(k)
          K_ij = BulkValues(k)
          rsum = rsum + K_ij * u(j) 
        END DO
        ku(i) = rsum
      END DO
    END IF

    CALL ListPushNameSpace('fct:')
    CALL ListAddLogical( Params,'fct: Skip Compute Nonlinear Change',.TRUE.)
    CALL ListAddLogical( Params,'fct: Skip Advance Nonlinear iter',.TRUE.)
  
    A % Values => M_C
    udot = 0._dp
    CALL SolveLinearSystem(A,Ku,Udot,Norm,1,Solver)
    CALL ListPopNamespace()

    A % Values => SaveValues
  END BLOCK
#endif

    ! Computation of correction factors (Zalesak's limiter)
    ! Code derived initially from Kuzmin's subroutine   
    !---------------------------------------------------------
    CALL Info('FCT_Correction','Compute correction factors',Level=10)
    pp = 0 
    pm = 0
    qp = 0 
    qm = 0

    IF(ParEnv % PEs>1) THEN
      fct_d => A % FCT_D
      mmc    => A % MassValues
      mmc_h  => A % HaloMassValues

      ALLOCATE(ActiveNodes(n)); activeNodes=.FALSE.
      DO i=1,Solver % NumberOfActiveElements
        Element => Solver % Mesh % Elements(Solver % ActiveElements(i))
        IF ( Element % PartIndex /= ParEnv % MyPE ) CYCLE
        Activenodes(Solver % Variable % Perm(Element % NodeIndexes)) = .TRUE.
      END DO
    ELSE
      fct_d => A % FCT_D
      mmc => A % MassValues
    END IF
    DO i=1,n
      IF (ParEnv % PEs > 1 ) THEN
        IF ( .NOT. ActiveNodes(i) ) CYCLE
      end if

      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)

        IF (ParEnv % PEs>1) THEN
          IF ( .NOT.ActiveNodes(j) ) CYCLE
        END IF

        ! Compute the raw antidiffusive fluxes
        ! f_ij=m_ij*[udot(i)-udot(j)]+d_ij*[ulow(i)-ulow(j)]
        !-----------------------------------------------------
        ! d_ij and m_ij are both symmetric
        ! Hence F_ji = -F_ij
           
        f_ij = mmc(k)*(udot(i)-udot(j)) + fct_d(k)*(u(i)-u(j))
        IF ( ParEnv % PEs>1 ) f_ij=f_ij+mmc_h(k)*(udot(i)-udot(j))
        ! Compared to Kuzmin's paper F_ij=-F_ij since d_ij and 
        ! udot have different signs. 
        f_ij = -f_ij

        ! Antidiffusive fluxes to be limited
        du = u(j)-u(i)

        ! Prelimiting of antidiffusive fluxes i.e. du and the flux have different signs
        IF (f_ij*du >= TINY(du)) THEN
          f_ij = 0._dp
        ELSE        
          ! Positive/negative edge contributions
          pp(i) = pp(i) + MAX(0._dp,f_ij)
          pm(i) = pm(i) + MIN(0._dp,f_ij)
        END IF

        ! Maximum/minimum solution increments
        qp(i) = MAX(qp(i),du)
        qm(i) = MIN(qm(i),du)
      END DO
    END DO

    ! Computation of nodal correction factors
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
    ! Symmetric flux limiting
    ! Correction of the right-hand side


    !   IF (ParEnv % PEs>1) THEN
    !     CALL ParallelSumVector(A,pm)
    !     CALL ParallelSumVector(A,pp)
    !     CALL ParallelSumVector(A,qm)
    !     CALL ParallelSumVector(A,qp)
    !   END IF
    
    CorrCoeff = ListGetCReal( Params,'FCT Correction Coefficient',Found )
    IF( .NOT. Found ) CorrCoeff = 1._dp

    Ceps = ListGetCReal( Params,'FCT Correction Epsilon',Found )
    IF(.NOT. Found ) Ceps = EPSILON( Ceps )
    corr = 0._dp
    DO i=1,n
      IF (ParEnv % PEs>1) THEN
        IF( .NOT. ActiveNodes(i)) CYCLE
      END IF

      IF( pp(i) > Ceps ) THEN
        rpi = MIN( 1._dp, M_L(i)*qp(i)/pp(i) )
      ELSE
        rpi = 0._dp
      END IF

      IF( pm(i) < -Ceps ) THEN
        rmi = MIN( 1._dp, M_L(i)*qm(i)/pm(i) )
      ELSE
        rmi = 0._dp
      END IF

      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)
        IF(ParEnv % PEs>1) THEN
          IF(.NOT.ActiveNodes(j)) CYCLE
        END IF

        f_ij = mmc(k)*(udot(i)-udot(j))  + fct_d(k)*(u(i)-u(j))
        IF (ParEnv % PEs>1) f_ij = f_ij + mmc_h(k)*(udot(i)-udot(j))
        f_ij = -f_ij

        IF (f_ij > 0) THEN 
          IF( pm(j) < -Ceps ) THEN
            rmj = MIN( 1.0_dp, M_L(j)*qm(j)/pm(j) )
          ELSE
            rmj = 0._dp
          END IF
          c_ij = MIN(rpi,rmj)
        ELSE 
          IF( pp(j) > Ceps ) THEN
            rpj = MIN( 1._dp, M_L(j)*qp(j)/pp(j) )
          ELSE
            rpj = 0._dp
          END IF
          c_ij = MIN(rmi,rpj)
        END IF
        corr(i) = corr(i) + c_ij * f_ij
      END DO
      corr(i) = CorrCoeff * corr(i) / M_L(i)
    END DO

    IF (ParEnv % PEs>1) THEN
!     CALL ParallelSumVector(A,corr)
      DEALLOCATE(A % HaloValues, A % HaloMassValues)
      A % HaloValues => Null(); A % HaloMassValues => Null()
    END IF

    ! Optionally skip applying the correction, just for debugging purposes
    IF( SkipCorrection ) THEN
      CALL Info('FCT_Correction','Skipping Applying corrector',Level=6)
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
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL, OPTIONAL :: Nonlinear, SteadyState

     LOGICAL :: IsNonlinear,IsSteadyState,Timing, RequireNonlinear, ContactBC, &
         MortarBC, IntegralBC
     LOGICAL :: ApplyMortar, ApplyContact, ApplyIntegral, StoreCyclic, Found, &
         StaticProj, IsBodyForce
     INTEGER :: i,j,k,l,n,dsize,size0,col,row,dim
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Matrix_t), POINTER :: CM, CMP, CM0, CM1
     TYPE(Variable_t), POINTER :: DispVar
     REAL(KIND=dp) :: t0,rt0,rst,st,ct
     CHARACTER(*), PARAMETER :: Caller = 'GenerateProjectors'
     TYPE(Solver_t), POINTER :: PSolver
     TYPE(Matrix_t), POINTER :: Proj
     
     ApplyIntegral = ListGetLogical( Solver % Values,'Apply Integral BCs',Found)
     ApplyMortar = ListGetLogical(Solver % Values,'Apply Mortar BCs',Found) 
     ApplyContact = ListGetLogical(Solver % Values,'Apply Contact BCs',Found)     
     
     IF( .NOT. ( ApplyMortar .OR. ApplyContact .OR. ApplyIntegral ) ) RETURN
     
     ! Here we give the option to block out cyclic projector if not wanted. 
     StoreCyclic = ListGetLogical( Solver % Values,'Store Cyclic Projector', Found)
     IF(.NOT. Found ) StoreCyclic = ListGetLogical( Solver % Values,'Store Cyclic System', Found)
     PSolver => Solver
     
     i = ListGetInteger( Solver % Values,'Mortar BC Master Solver',Found ) 
     IF( Found ) THEN
       Solver % MortarBCs => CurrentModel % Solvers(i) % MortarBCs
       IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) THEN
         CALL Fatal(Caller,'Could not reuse projectors from solver: '//I2S(i))
       END IF
       CALL Info(Caller,'Reusing projectors from solver: '//I2S(i),Level=8)
       RETURN
     END IF

     CALL Info(Caller,'Generating various boundary projectors',Level=8)

     Timing = ListCheckPrefix(Solver % Values,'Projector Timing')
     IF( Timing ) THEN
       t0 = CPUTime(); rt0 = RealTime()      
     END IF

     IsNonlinear = .FALSE.
     IF( PRESENT( Nonlinear ) ) IsNonlinear = Nonlinear
     IsSteadyState = .NOT. IsNonlinear

     IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) THEN
       ALLOCATE( Solver % MortarBCs( Model % NumberOfBCs + Model % NumberOfBodyForces ) )
       DO i=1, Model % NumberOfBCs + Model % NumberOfBodyForces
         Solver % MortarBCs(i) % Projector => NULL()
       END DO
     END IF
     
     dim = CoordinateSystemDimension()

     DO i=1,Model % NumberOFBCs + Model % NumberOfBodyForces
       IF(i > Model % NumberOfBCs ) THEN
         BC => Model % BodyForces(i-Model % NumberOfBCs) % Values
         IsBodyForce = .TRUE.
       ELSE         
         BC => Model % BCs(i) % Values
         IsBodyForce = .FALSE.
       END IF
         
       k = 0
       j = ListGetInteger( BC,'Mortar BC',MortarBC)       
       j = j + ListGetInteger( BC,'Contact BC',ContactBC)       
       IntegralBC = ListGetLogical( BC,'Integral BC',Found ) 

       IF( MortarBC ) k = k+1
       IF( ContactBC ) k = k+1
       IF( IntegralBC ) k = k+1
       IF(k==0) CYCLE
         
       IF( k > 1 ) THEN
         CALL Fatal(Caller,'Boundary '//I2S(i)//' can only be one of mortar, contact and integral!')
       END IF     

       IF(IsBodyForce .AND. .NOT. IntegralBC ) THEN
         CALL Fatal(Caller,'Body Force '//I2S(i)//' can only have integral bc!')
       END IF     

       IF( InfoActive(10) ) THEN
         IF( MortarBC ) CALL Info(Caller,'Generating mortar conditions for BC: '//I2S(i))
         IF( ContactBC ) CALL Info(Caller,'Generating contact conditions for BC: '//I2S(i))
         IF( IntegralBC ) CALL Info(Caller,'Generating integral conditions for BC: '//I2S(i))
       END IF
       
       RequireNonlinear = ListGetLogical( BC,'Mortar BC Nonlinear',Found)
       IF( .NOT. Found ) THEN
         RequireNonlinear = ContactBC .AND. .NOT. ListGetLogical( BC,'Tie Contact',Found )
       END IF

       IF( IntegralBC ) RequireNonlinear = .FALSE.
       
       IF( IsNonlinear ) THEN
         IF( .NOT. RequireNonlinear ) CYCLE
       ELSE
         IF( RequireNonlinear ) CYCLE
       END IF             

       StaticProj = ListGetLogical( BC,'Mortar BC Static',Found)
       
       Proj => Solver % MortarBCs(i) % Projector
       IF( ASSOCIATED( Proj ) ) THEN
         IF( StaticProj ) CYCLE         

         IF( StoreCyclic ) THEN
           ! Don't release projectors in case they are cyclic 
           ! Instead reassign the pointer.
           CALL StoreCyclicProjector(PSolver,Proj,Found)
           IF(Found) THEN
             Solver % MortarBCs(i) % Projector => Proj
             Solver % MortarBCsChanged = .TRUE.
             CYCLE
           END IF
         ELSE  
           IF( ASSOCIATED( Proj % Ematrix ) ) THEN
             CALL FreeMatrix( Proj % Ematrix )
           END IF
           CALL FreeMatrix( Proj ) 
         END IF
       END IF

       ! Compute new projector
       IF( IntegralBC ) THEN         
         Proj => IntegralProjector(Model,Solver % Mesh, i, IsBodyForce ) 
       ELSE
         ! This is the same for mortar and contact!
         Proj => PeriodicProjector(Model,Solver % Mesh,i,j,dim,.TRUE.)
       END IF
         
       Solver % MortarBCs(i) % Projector => Proj       
       IF( ASSOCIATED( Proj ) ) THEN
         Solver % MortarBCsChanged = .TRUE.
       END IF

       ! Store new projector to the cyclic set
       IF( StoreCyclic ) THEN
         IF(.NOT. StaticProj ) CALL StoreCyclicProjector(PSolver,Proj)
       END IF
       
     END DO


     IF( Timing ) THEN
       st  = CPUTime() - t0;
       rst = RealTime() - rt0
       
       WRITE(Message,'(a,f8.2,f8.2,a)') 'Projector creation time (CPU,REAL) for '&
           //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
       CALL Info(Caller,Message,Level=6)    
       
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



  ! Create special type of projector after the linear system has been created. 
  ! -------------------------------------------------------------------------
   SUBROUTINE GenerateRobinProjectors(Model,Solver)

     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver

     TYPE(Solver_t), POINTER :: PSolver
     INTEGER :: i
     LOGICAL :: Found, GotSome, DoIt, IsBodyForce
     CHARACTER(*), PARAMETER :: Caller="GenerateRobinProjector"

     
     IF(.NOT. ListGetLogical( Solver % Values,'Apply Integral BCs',Found) ) RETURN

     PSolver => Solver
     GotSome = .FALSE.
     DO i=1,Model % NumberOFBCs + Model % NumberOfBodyForces
       IF(i<=Model % NumberOfBCs) THEN
         IF(.NOT. ListGetLogical( Model % BCs(i) % Values,'Flux Integral BC',Found ) ) CYCLE
         IsBodyForce = .FALSE.
       ELSE
         IF(.NOT. ListGetLogical( Model % BodyForces(i-Model % NumberOfBCs) % Values,'Flux Integral BC',Found ) ) CYCLE
         IsBodyForce = .TRUE.
       END IF
       
       CALL Info(Caller,'Generating flux integral conditions for BC: '//I2S(i))
       CALL RobinProjector(Model,PSolver, i, IsBodyForce) 
       GotSome = .TRUE.
     END DO

     IF( GotSome ) THEN
       Solver % MortarBCsChanged = .TRUE.
       ! We may want to export the lagrange multiplier as it has a physical meaning. 
       CALL ListAddNewLogical(Solver % Values,'Export Lagrange Multiplier',.TRUE.) 
     END IF
            
   END SUBROUTINE GenerateRobinProjectors



   ! This creates a projector that integrates over the BCs on the boundary such that
   ! an integral constraint for Robin type of BCs may be applied.
   !--------------------------------------------------------------------------------------
   SUBROUTINE RobinProjector(Model, Solver, BCInd, IsBodyForce ) 

     TYPE(Model_t) :: Model
     TYPE(Solver_t), POINTER :: Solver    
     INTEGER :: BCInd
     LOGICAL :: IsBodyForce

     TYPE(Matrix_t), POINTER :: Proj        
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(ValueList_t), POINTER :: BC
     LOGICAL :: Found
     INTEGER :: dofs
     TYPE(MortarBC_t), POINTER :: MortarBC
     CHARACTER(*), PARAMETER :: Caller="RobinProjector"

     IF(IsBodyForce) THEN
       BC => Model % BodyForces(BCind-Model % NumberOfBCs) % Values
     ELSE
       BC => Model % BCs(BCInd) % Values
     END IF
     IF( .NOT. ListGetLogical( BC,'Flux Integral BC', Found ) ) RETURN

     Proj => Solver % MortarBCs(BCInd) % Projector
     IF(.NOT. ASSOCIATED(Proj)) THEN
       Proj => AllocateMatrix()
       Proj % FORMAT = MATRIX_LIST
       Proj % ProjectorType = PROJECTOR_TYPE_ROBIN
       CALL Info(Caller,'Creating Robin flux constraint matrix for boundary: '//I2S(BCind),Level=6)
     ELSE
       CALL Info(Caller,'Updating Robin flux constraint matrix for boundary: '//I2S(BCind),Level=6)
       Proj % Values = 0.0_dp
     END IF

     Mesh => Solver % Mesh

     CALL CreateRobinProjector()

     IF(Proj % FORMAT == MATRIX_LIST ) THEN
       CALL List_toCRSMatrix(Proj)
     END IF

     IF( InfoActive(20) ) THEN
       WRITE(Message,'(A,ES12.3)') 'Sum of constraint matrix entries: ',SUM(Proj%Values)
       CALL Info(Caller,Message)
       CALL Info(Caller,'Constraint matrix cols min: '//I2S(MINVAL(Proj%Cols)))
       CALL Info(Caller,'Constraint matrix cols max: '//I2S(MAXVAL(Proj%Cols)))
       CALL Info(Caller,'Constraint matrix rows min: '//I2S(MINVAL(Proj%Rows)))
       CALL Info(Caller,'Constraint matrix rows max: '//I2S(MAXVAL(Proj%Rows)))
     END IF

     MortarBC => Solver % MortarBCs(BCind) 
     MortarBC % Projector => Proj       

     IF(.NOT. ASSOCIATED(MortarBC % Diag ) ) THEN      
       dofs = Solver % Variable % Dofs
       ALLOCATE(MortarBC % Diag(dofs))
       MortarBC % Diag = 1.0_dp
     END IF

     
   CONTAINS

     SUBROUTINE CreateRobinProjector()

       INTEGER :: i,j,k,n,t,dofs,idof,t1,t2
       REAL(KIND=dp) :: dval
       TYPE(Matrix_t), POINTER :: A
       TYPE(Variable_t), POINTER :: Var
       TYPE(Element_t), POINTER :: Element
       LOGICAL, ALLOCATABLE :: ActiveDof(:)
       INTEGER, POINTER :: Indexes(:)
       
       A => Solver % Matrix
       Var => Solver % Variable       
       n = A % NumberOfRows
       dofs = Var % dofs

       IF(.NOT. ASSOCIATED(A % BulkValues)) THEN
         CALL Fatal(Caller,'BulkValues are needed to aveluate Robin terms!')
       END IF

       IF(IsBodyForce) THEN
         t1 = 1
         t2 = Mesh % NumberOfBulkElements 
       ELSE
         t1 = Mesh % NumberOfBulkElements + 1
         t2 = (t1-1) + Mesh % NumberOfBoundaryElements
       END IF
       
       ! Mark the dofs of the matrix that are on the boundary.
       ! ActiveDof table will directly refer to the indexes of the matrix.
       ALLOCATE(ActiveDof(n))
       ActiveDof = .FALSE.       
       DO t = t1, t2
         Element => Mesh % Elements(t)
         IF(IsBodyForce) THEN
           i = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force',Found)
           IF(i /= BCind - Model % NumberOfBCs) CYCLE
         ELSE
           IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BCInd) % Tag ) CYCLE
         END IF
         Indexes => Element % NodeIndexes      

         IF(ANY(Var % Perm(Indexes) == 0 ) ) CYCLE
         DO i=1,dofs
           ActiveDof(dofs*(Var % Perm(Indexes)-1)+i) = .TRUE.
         END DO
       END DO
              
       DO i=1, n / dofs 
         DO idof = 1, dofs
           j = dofs*(i-1)+idof
           IF(.NOT. ActiveDof(j)) CYCLE
           
           DO k=A % Rows(j),A % Rows(j+1)-1
             IF(.NOT. ActiveDof(A % Cols(k))) CYCLE
             dval = A % Values(k) - A % BulkValues(k)
             ! We use the negative sign so that we get desired sign for the lagrange multiplier!
             CALL AddToMatrixElement(Proj, idof, A % Cols(k), -dval )
           END DO
         END DO
       END DO       
       
     END SUBROUTINE CreateRobinProjector

   END SUBROUTINE RobinProjector
   

   ! Generate constraint matrix from mortar projectors. 
   ! This routine takes each boundary projector and applies it 
   ! to the current field variable (scalar or vector) merging 
   ! all into one single projector. 
   !---------------------------------------------------------
   SUBROUTINE GenerateConstraintMatrix( Model, Solver )

     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver

     INTEGER, POINTER :: Perm(:)
     INTEGER :: i,j,j2,k,k2,l,l2,dofs,maxperm,permsize,bc_ind,constraint_ind,row,col,col2,&
         mcount,bcount,kk,cdofs,dim
     TYPE(Matrix_t), POINTER :: Atmp,Btmp, Ctmp
     LOGICAL :: AllocationsDone, CreateSelf, ComplexMatrix, TransposePresent, Found, &
         SetDof, SomeSet, SomeSkip, SumProjectors, NewRow, SumThis
     INTEGER, ALLOCATABLE :: SumPerm(:),SumCount(:)
     LOGICAL, ALLOCATABLE :: ActiveComponents(:), SetDefined(:)
     TYPE(ValueList_t), POINTER :: BC
     TYPE(MortarBC_t), POINTER :: MortarBC
     REAL(KIND=dp) :: wsum, Scale
     INTEGER :: rowoffset, arows, sumrow, EliminatedRows, NeglectedRows, sumrow0, k20
     LOGICAL :: ThisIsMortar, ThisIsRobin, Reorder
     LOGICAL :: AnyPriority
     INTEGER :: Priority, PrevPriority
     INTEGER, ALLOCATABLE :: BCOrdering(:), BCPriority(:)
     LOGICAL :: NeedToGenerate, ComplexSumRow 

     LOGICAL :: HaveMortarDiag, LumpedDiag, PerFlipActive
     REAL(KIND=dp) :: MortarDiag, val, valsum, EpsVal
     LOGICAL, POINTER :: PerFlip(:)
     CHARACTER(*), PARAMETER :: Caller = 'GenerateConstraintMatrix'

     LOGICAL :: IntegralBC
     REAL(KIND=dp) :: SetVal(6)
     REAL(KIND=dp), ALLOCATABLE :: PrevValues(:)
     INTEGER, ALLOCATABLE :: PrevInvPerm(:)
     TYPE(Variable_t), POINTER :: Var
     CHARACTER(:), ALLOCATABLE :: Str,MultName
     REAL(KIND=dp), ALLOCATABLE :: rsum(:)
     LOGICAL :: IsDg, IsBodyForce
     INTEGER, ALLOCATABLE :: DgSome(:)
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Element
     
     ! Should we genarete the matrix
     NeedToGenerate = Solver % MortarBCsChanged
     
     Mesh => Solver % Mesh 
     IsDg = Solver % DG     

     PerFlipActive = Solver % PeriodicFlipActive
     IF( PerFlipActive ) THEN
       CALL Info(Caller,'Periodic flip is active',Level=8)
       PerFlip => Mesh % PeriodicFlip           
     END IF
     
     ! Set pointers to save the initial constraint matrix
     ! We assume that the given ConstraintMatrix is constant but we have consider it the 1st time
     IF(.NOT. Solver % ConstraintMatrixVisited ) THEN       
       IF( ASSOCIATED( Solver % Matrix % ConstraintMatrix ) ) THEN
         CALL Info(Caller,'Saving initial constraint matrix to Solver',Level=12)
         Solver % ConstraintMatrix => Solver % Matrix % ConstraintMatrix
         Solver % Matrix % ConstraintMatrix => NULL()
         NeedToGenerate = .TRUE. 
       END IF
       Solver % ConstraintMatrixVisited = .TRUE.
     END IF
     
     IF( NeedToGenerate ) THEN
       CALL Info(Caller,'Building constraint matrix',Level=12)
     ELSE     
       CALL Info(Caller,'Nothing to do for now',Level=12)
       RETURN
     END IF

     
     ! Compute the number and size of initial constraint matrices
     !-----------------------------------------------------------
     row    = 0
     mcount = 0
     bcount = 0
     Ctmp => Solver % ConstraintMatrix
     IF( ASSOCIATED( Ctmp ) ) THEN
       DO WHILE(ASSOCIATED(Ctmp))
         mcount = mcount + 1
         row = row + Ctmp % NumberOfRows
         Ctmp => Ctmp % ConstraintMatrix

         IF( InfoActive(32) ) THEN           
           CALL VectorValuesRange(Ctmp % Values,SIZE(Ctmp % Values),'Ctmp'//I2S(mcount))
           IF( ASSOCIATED( Ctmp % InvPerm ) ) THEN
             PRINT *,'InvPerm range:',MINVAL(Ctmp % InvPerm), MAXVAL(Ctmp % InvPerm), SUM(Ctmp % InvPerm)
           END IF
         END IF
    
       END DO
       CALL Info(Caller,'Number of initial constraint matrices: '//I2S(mcount),Level=12)       
       CALL Info(Caller,'Number of rows in constraint matrices: '//I2S(row),Level=20)       
     END IF
       
     
     ! Compute the number and size of mortar matrices
     !-----------------------------------------------
     IF( ASSOCIATED( Solver % MortarBCs ) ) THEN
       DO bc_ind=1,Model % NumberOFBCs + Model % NumberOfBodyForces
         Atmp => Solver % MortarBCs(bc_ind) % Projector         
         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE
         bcount = bcount + 1
         row = row + Atmp % NumberOfRows
       END DO
       CALL Info(Caller,'Number of mortar matrices: '//I2S(bcount),Level=12)       
     END IF
     
     IF( row==0 ) THEN
       CALL Info(Caller,'Nothing to do since there are no constrained dofs!',Level=12)       
       RETURN
     END IF

     MortarDiag = ListGetCReal( Solver % Values,'Mortar Diag',HaveMortarDiag )
     LumpedDiag = ListGetLogical( Solver % Values,'Lumped Diag',Found )

     IF( HaveMortarDiag ) THEN
       CALL Info(Caller,'Adding diagonal entry to mortar constraint!',Level=12)              
     END IF
     
     IF( mcount == 1 .AND. bcount == 0 .AND. .NOT. HaveMortarDiag ) THEN
       CALL Info(Caller,'Using initial constraint matrix',Level=12)       
       Solver % Matrix % ConstraintMatrix => Solver % ConstraintMatrix
       RETURN
     END IF

     ! Now we are generating something more complex and different than last time
     IF( ASSOCIATED( Solver % Matrix % ConstraintMatrix ) ) THEN
       IF ( ListGetLogical( Solver % Values,'Apply Contact BCs', Found ) ) THEN
         CALL Info(Caller,'Remember the previous InvPerm for contact mechanics',Level=20)
         ALLOCATE( PrevInvPerm( SIZE( Solver % Matrix % ConstraintMatrix % InvPerm ) ) )
         PrevInvPerm = Solver % Matrix % ConstraintMatrix % InvPerm       
       END IF
       
       CALL Info(Caller,'Releasing previous constraint matrix',Level=12)     
       CALL FreeMatrix(Solver % Matrix % ConstraintMatrix)
       Solver % Matrix % ConstraintMatrix => NULL()
     END IF
       
     EpsVal = ListGetConstReal( Solver % Values,&
         'Minimum Projector Value', Found )
     IF(.NOT. Found ) EpsVal = 1.0d-8
     
     
     SumProjectors = ListGetLogical( Solver % Values,&
         'Mortar BCs Additive', Found )
     IF( .NOT. Found ) THEN
       IF( bcount > 1 .AND. ListGetLogical( Solver % Values, &
           'Eliminate Linear Constraints',Found ) ) THEN
         CALL Info(Caller,'Enforcing > Mortar BCs Additive < to True to enable elimination',Level=8)
         SumProjectors = .TRUE.
       END IF       
       IF( .NOT. SumProjectors .AND. ListGetLogical( Solver % Values, &
           'Apply Conforming BCs',Found ) ) THEN
         CALL Info(Caller,'Enforcing > Mortar BCs Additive < to True because of conforming BCs',Level=8)
         SumProjectors = .TRUE.
       END IF
     END IF
     EliminatedRows = 0

     CALL Info(Caller,'There are '&
         //I2S(row)//' initial rows in constraint matrices',Level=10)

     dim = Mesh % MeshDim              
     dofs = Solver % Variable % DOFs
     
     Perm => Solver % Variable % Perm
     permsize = SIZE( Perm )
     maxperm  = MAXVAL( Perm )
     AllocationsDone = .FALSE.
     arows = Solver % Matrix % NumberOfRows

     ! Create a table that shows one way how continuous nodal dofs maps to
     ! DG nodal dofs. Only one is needed since we assume reduced basis!
     IF( IsDG ) THEN
       ALLOCATE(DgSome(permsize))
       DgSome = 0

       DO i=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(i)
         DO j=1,Element % TYPE % NumberOfNodes
           k = Element % NodeIndexes(j)
           k2 = Element % DGIndexes(j)
           DgSome(k) = k2
         END DO
       END DO       
     END IF

       
     ComplexMatrix = Solver % Matrix % Complex
     ComplexSumRow = .FALSE.
     
     IF( ComplexMatrix ) THEN
       IF( MODULO( Dofs,2 ) /= 0 ) CALL Fatal(Caller,&
           'Complex matrix should have even number of components!')
     ELSE
       ! Currently complex matrix is enforced if there is an even number of 
       ! entries since it seems that we cannot rely on the flag to be set.
       ComplexMatrix = ListGetLogical( Solver % Values,'Linear System Complex',Found )
       IF( .NOT. Found ) ComplexMatrix = ( Dofs == 2*dim) 
     END IF

     IF( ComplexMatrix ) THEN
       IF(dofs==dim .OR. dofs == 2) THEN
         cdofs = dofs
       ELSE
         CALL Fatal(Caller,'Invalid number of dofs for field: '//I2S(dofs))
       END IF
     ELSE     
       IF(dofs==dim .OR. dofs == 1) THEN
         cdofs = dofs
       ELSE IF(dofs==dim+1) THEN
         ! For contact mechanics we want to ignore the pressure. 
         IF( ListGetLogical( Solver % Values,'Apply Contact BCs',Found ) ) THEN
           cdofs = dim
         ELSE
           cdofs = dofs
         END IF
       ELSE
         CALL Fatal(Caller,'Invalid number of dofs for field: '//I2S(dofs))
       END IF
     END IF
     
     ALLOCATE( ActiveComponents(dofs), SetDefined(dofs), rsum(dofs) ) 
     
     IF( SumProjectors ) THEN
       ALLOCATE( SumPerm( dofs * permsize ) )
       SumPerm = 0
       ALLOCATE( SumCount( arows ) )
       SumCount = 0
     END IF
          
     AnyPriority = ListCheckPresentAnyBC( Model,'Projector Priority') 
     IF( AnyPriority ) THEN
       IF(.NOT. SumProjectors ) THEN
         CALL Warn(Caller,'Priority has effect only in additive mode!')
         AnyPriority = .FALSE.
       ELSE
         CALL Info(Caller,'Using priority for projector entries',Level=7)
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
     sumrow0 = 0
     k20 = 0
     
     TransposePresent = .FALSE.
     Ctmp => Solver % ConstraintMatrix

     DO constraint_ind = Model % NumberOFBCs + Model % NumberOfBodyForces + mcount,1,-1
       
       ! This is the default i.e. all components are applied mortar BCs
       ActiveComponents = .TRUE.
       ThisIsRobin = .FALSE.
       
       IF(constraint_ind > Model % NumberOfBCs + Model % NumberOfBodyForces ) THEN
         ThisIsMortar = .FALSE.
         Reorder = .FALSE.                
         SumThis = .FALSE.
         Atmp => Ctmp
         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE
         Ctmp => Ctmp % ConstraintMatrix
         IF( .NOT. ASSOCIATED( Atmp % InvPerm ) ) THEN
           IF(.NOT. AllocationsDone ) THEN
             CALL Warn(Caller,'InvPerm is expected, using identity!')
           END IF
         END IF
         CALL Info(Caller,'Adding initial constraint matrix: '&
             //I2S(constraint_ind - Model % NumberOfBCs),Level=8)         
       ELSE
         ThisIsMortar = .TRUE.
         
         ! Assume the mortar matrices refer to unordered mesh dofs
         ! and existing ConstraintMatrix to already ordered entities. 
         Reorder = ThisIsMortar

         SumThis = SumProjectors
         IsBodyForce = .FALSE.
         IF( constraint_ind > Model % NumberOfBCs ) THEN
           IsBodyForce = .TRUE.
           bc_ind = constraint_ind
         ELSE IF( AnyPriority ) THEN           
           bc_ind = BCOrdering(constraint_ind)
         ELSE
           bc_ind = constraint_ind 
         END IF
         
         MortarBC => Solver % MortarBCs(bc_ind) 
         Atmp => MortarBC % Projector

         ! Add the number of rows already populated in the constraint matrix so we can associate
         ! the single constraint to the correct entry in the constraint matrix.
         MortarBC % rowoffset = sumrow

         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE

         IF(IsBodyForce ) THEN
           BC => Model % BodyForces(bc_ind - Model % NumberOfBCs) % Values
           Priority = 0
         ELSE
           BC => Model % BCs(bc_ind) % Values         
           IF( AnyPriority ) THEN
             Priority = ListGetInteger( BC,'Projector Priority',Found)
           END IF
         END IF
           
         IntegralBC = ListGetLogical( BC,'Integral BC',Found ) 

         IF( Atmp % ProjectorType == PROJECTOR_TYPE_ROBIN ) THEN
           ThisIsRobin = .TRUE.
           Reorder = .FALSE.         
         ELSE
           IF( .NOT. ASSOCIATED( Atmp % InvPerm ) ) THEN
             CALL Fatal(Caller,'InvPerm is required for geometric projector!')
           END IF
         END IF

         IF(.NOT. AllocationsDone ) THEN
           IF( ThisIsRobin ) THEN
             CALL Info(Caller,'Adding flux constraint for BC: '//I2S(bc_ind),Level=8)
           ELSE IF( IntegralBC ) THEN
             CALL Info(Caller,'Adding integral constraint for BC: '//I2S(bc_ind),Level=8)
           ELSE             
             CALL Info(Caller,'Adding mortar projector for BC: '//I2S(bc_ind),Level=8)
           END IF
           CALL Info(Caller,'Adding projector rows: '//I2S(Atmp % NumberOfRows),Level=12)
         END IF           
         
         ! Enable that the user can for vector valued cases either set some 
         ! or skip some field components. 
         SomeSet = .FALSE.
         SomeSkip = .FALSE.
         DO i=1,cDofs
           IF( cDofs > 1 ) THEN
             str = ComponentName( Solver % Variable, i )
           ELSE
             str = TRIM(Solver % Variable % Name)
           END IF

           IF( IntegralBC ) THEN
             SetVal(i) = ListGetCReal( BC,'Integral BC '//TRIM(str),SetDof )
             IF(ThisIsRobin) SetDof = .TRUE.
           ELSE
             SetDof = ListGetLogical( BC,'Mortar BC '//TRIM(str),Found )
           END IF

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
         IF( SomeSet .AND. .NOT. ALL(SetDefined(1:cdofs)) ) THEN
           IF( SomeSkip ) THEN
             CALL Fatal(Caller,'Do not know what to do with all '//I2S(cdofs)//' components')
           ELSE
             CALL Info(Caller,'Unspecified components will not be set for BC '//I2S(bc_ind),Level=10)
             DO i=1,cDofs
               IF( .NOT. SetDefined(i) ) ActiveComponents(i) = .FALSE.
             END DO
           END IF
         END IF
       END IF

       TransposePresent = TransposePresent .OR. ASSOCIATED(Atmp % Child)
       IF( TransposePresent ) THEN
         IF(ASSOCIATED(Atmp % Child) ) THEN
           IF(.NOT. ASSOCIATED(Atmp % Child % Values)) THEN
             CALL Fatal(Caller,'Atmp has Child which has no values!')
           END IF
         END IF
         CALL Info(Caller,'Transpose matrix is present',Level=8)
       END IF

       ! If the projector is of type x_s=P*x_m then generate a constraint matrix
       ! of type [D-P]x=0 where D is diagonal unit matrix. 
       CreateSelf = ( Atmp % ProjectorType == PROJECTOR_TYPE_NODAL ) 
       
       IF( SumThis .AND. CreateSelf ) THEN
         CALL Fatal(Caller,'It is impossible to sum up nodal projectors!')
       END IF

       ComplexSumRow = ListGetLogical( Solver % Values,'Complex Sum Row ', Found )
       IF(.NOT. Found ) THEN       
         ComplexSumRow = ( dofs == 2 .AND. ComplexMatrix .AND. .NOT. CreateSelf .AND. &
             SumThis .AND. .NOT. (ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag ) )
       END IF
       
       
       ! We deal with the Robin Flux cBC's here even though they would be associated 
       ! to vector or complex valued field. 
       IF( Dofs == 1 .OR. ThisIsRobin ) THEN         
         
         IF( .NOT. ActiveComponents(1) ) THEN
           CALL Info(Caller,'Skipping component: '//I2S(1),Level=12)
           CYCLE
         END IF

         ! Number the rows. 
         IF( SumThis ) THEN
           DO i=1,Atmp % NumberOfRows                               
             ! Skip empty row
             IF( Atmp % Rows(i) >= Atmp % Rows(i+1) ) CYCLE 

             ! If the mortar boundary is not active at this round don't apply it
             IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Active ) ) THEN
                 IF( .NOT. MortarBC % Active(i) ) CYCLE
               END IF
             END IF
             
             ! Use InvPerm if it is present
             IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
               k = Atmp % InvPerm(i)
               ! Node does not have an active dof to be constrained
               IF( k == 0 ) CYCLE
             ELSE
               k = i
             END IF

             kk = k             
             IF( Reorder ) THEN
               IF(IsDG) THEN
                 kk = Perm(DgSome(k))
               ELSE
                 kk = Perm(k)
               END IF
               IF( kk == 0 ) CYCLE
             END IF
             
             NewRow = ( SumPerm(kk) == 0 )
             IF( NewRow ) THEN
               sumrow = sumrow + 1                
               SumPerm(kk) = sumrow 
             ELSE IF(.NOT. AllocationsDone ) THEN
               IF( Priority /= PrevPriority .AND. SumPerm(kk) < 0 ) THEN
                 NeglectedRows = NeglectedRows + 1
               ELSE
                 EliminatedRows = EliminatedRows + 1
               END IF
             END IF
           END DO
           CALL Info(Caller,'Number of rows: '//I2S(sumrow),Level=20)         
         END IF
         
         IF( ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag) THEN
           CALL Info(Caller,'MotarBC diag exists!',Level=30)
           IF( ASSOCIATED(Atmp % InvPerm) ) THEN
             CALL Info(Caller,'MotarBC InvPerm exists!',Level=30)
             IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN                   
               k = MAXVAL( Atmp % Cols )
               ALLOCATE( MortarBC % Perm(k) )
               MortarBC % Perm = 0
               DO k=1,SIZE(Atmp % InvPerm )
                 j = Atmp % InvPerm(k)
                 MortarBC % Perm( j ) = k
               END DO
             END IF
           END IF
         END IF


         DO i=1,Atmp % NumberOfRows                     
           
           IF( Atmp % Rows(i) >= Atmp % Rows(i+1) ) CYCLE ! skip empty rows

           ! If the mortar boundary is not active at this round don't apply it
           IF( ThisIsMortar ) THEN
             IF( ASSOCIATED( MortarBC % Active ) ) THEN
               IF( .NOT. MortarBC % Active(i) ) CYCLE
             END IF
           END IF
             
           IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
             k = Atmp % InvPerm(i)
             IF( k == 0 ) CYCLE
           ELSE
             k = i
           END IF
           
           kk = k
           IF( Reorder ) THEN
             IF(IsDg) THEN
               kk = Perm(DgSome(k)) 
             ELSE
               kk = Perm(k)
             END IF
             IF( kk == 0 ) CYCLE
           END IF
             
           IF( SumThis ) THEN             
             row = SumPerm(kk)
               
             ! Mark this for future contributions so we know this is already set
             ! and can skip this above.
             IF( AnyPriority ) THEN
               IF( row < 0 ) CYCLE
               IF( Priority /= PrevPriority ) SumPerm(kk) = -SumPerm(kk)
             END IF
             
             IF( row <= 0 ) THEN
               CALL Fatal(Caller,'Invalid row index: '//I2S(row))
             END IF
           ELSE
             sumrow = sumrow + 1
             row = sumrow
           END IF

           ! This fixes a problem for Robin type of constraints.
           ! Rather than understanding this just omits the problem...
           IF( AllocationsDone .AND. .NOT. ThisIsRobin ) THEN
             Btmp % InvPerm(row) = rowoffset + kk
           END IF
           
           wsum = 0.0_dp
           rsum = 0.0_dp
           
           valsum = 0.0_dp
           DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1             
             valsum = valsum + ABS( Atmp % Values(l) ) 
           END DO             

           DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1
             
             col = Atmp % Cols(l) 
             val = Atmp % Values(l)

             IF( ABS( val ) < EpsVal * valsum ) CYCLE
             
             IF( Reorder ) THEN               
               IF( col <= permsize ) THEN
                 IF(IsDg) THEN
                   col2 = Perm(DgSome(col))
                 ELSE
                   col2 = Perm(col)
                 END IF
                 IF( col2 == 0 ) CYCLE
               ELSE
                 CALL Fatal(Caller,'col index too large: '//I2S(col)//' vs '//I2S(permsize))
               END IF
             ELSE
               col2 = col
             END IF
               
             IF( AllocationsDone ) THEN
               ! By Default there is no scaling
               Scale = 1.0_dp
               IF( ThisIsMortar ) THEN
                 IF( CreateSelf ) THEN
                   ! We want to create [D-P] hence the negative sign
                   Scale = MortarBC % MasterScale
                   wsum = wsum + val
                 ELSE IF( ASSOCIATED( MortarBC % Perm ) ) THEN
                   ! Look if the component refers to the slave
                   IF( MortarBC % Perm( col ) > 0 ) THEN
                     Scale = MortarBC % SlaveScale 
                     wsum = wsum + val
                   ELSE
                     Scale = MortarBC % MasterScale
                   END IF
                 ELSE IF( ThisIsRobin ) THEN
                   j = MODULO(col-1,dofs)+1
                   rsum(j) = rsum(j) + val
                 ELSE
                   wsum = wsum + val
                 END IF

                 ! If we sum up to anti-periodic dof then use different sign
                 ! - except if the target is also antiperiodic.
                 IF( PerFlipActive ) THEN
                   IF(  PerFlip(col) .NEQV. PerFlip(k) ) Scale = -Scale
                 END IF                 
               END IF

               ! Add a new column index to the summed up row               
               ! At the first sweep we need to find the first unset position
               IF( SumThis ) THEN
                 k2 = Btmp % Rows(row)
                 DO WHILE( Btmp % Cols(k2) > 0 )
                   k2 = k2 + 1
                 END DO
               ELSE
                 k2 = k2 + 1
               END IF
               
               Btmp % Cols(k2) = col2
               Btmp % Values(k2) = Scale * val
               IF(ASSOCIATED(Btmp % TValues)) THEN
                 IF(ASSOCIATED(Atmp % Child)) THEN
                   Btmp % TValues(k2) = Scale * Atmp % Child % Values(l)
                 ELSE
                   Btmp % TValues(k2) = Scale * val
                 END IF
               END IF
             ELSE
               k2 = k2 + 1
               IF( SumThis ) THEN
                 SumCount(row) = SumCount(row) + 1
               END IF
             END IF
           END DO
           
           ! Add the self entry as in 'D'
           IF( CreateSelf ) THEN
             k2 = k2 + 1
             IF( AllocationsDone ) THEN
               IF(IsDG) THEN
                 Btmp % Cols(k2) = Perm( DGSome( Atmp % InvPerm(i) ) )
               ELSE
                 Btmp % Cols(k2) = Perm( Atmp % InvPerm(i) )
               END IF
               Btmp % Values(k2) = MortarBC % SlaveScale * wsum
             ELSE               
               IF( SumThis) SumCount(row) = SumCount(row) + 1
             END IF
           END IF
           
           ! Add a diagonal entry if requested. When this is done at the final stage
           ! all the hassle with the right column index is easier.
           IF( ThisIsMortar .OR. ThisIsRobin ) THEN
             diag: IF( ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag ) THEN
               IF( .NOT. HaveMortarDiag ) THEN
                 MortarDiag = MortarBC % Diag(i)
                 LumpedDiag = MortarBC % LumpedDiag
               END IF

               IF( ThisIsRobin ) THEN
                 DO j=1,dofs
                   k2 = k2 + 1
                   IF( AllocationsDone ) THEN
                     Btmp % Cols(k2) = j + arows                      
                     Btmp % Values(k2) = Btmp % Values(k2) - MortarDiag * rsum(j)
                   END IF                     
                 END DO
               
               ELSE IF( LumpedDiag ) THEN
                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   Btmp % Cols(k2) = row + arows 
                   Btmp % Values(k2) = Btmp % Values(k2) - MortarDiag * wsum
                 ELSE
                   IF( SumThis) SumCount(row) = SumCount(row) + 1
                 END IF
               ELSE
                 IF(ThisIsMortar .AND. .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN                   
                   CALL Fatal(Caller,'MortarBC % Perm required, try lumped')
                 END IF
                 
                 DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1                 
                   col = Atmp % Cols(l) 

                   IF( Reorder ) THEN
                     IF( col > permsize ) THEN
                       PRINT *,'col too large',col,permsize
                       CYCLE
                     END IF
                     IF(IsDg) THEN
                       col2 = Perm(DgSome(col))
                     ELSE
                       col2 = Perm(col)
                     END IF
                   ELSE
                     col2 = col
                   END IF                     
                   IF( col2 == 0 ) CYCLE
                     
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
                     IF( SumThis ) THEN
                       l2 = ABS( SumPerm( col2) )
                     ELSE
                       l2 = MortarBC % Perm(col)
                     END IF
                     
                     Btmp % Cols(k2) = l2 + arows + rowoffset
                     Btmp % Values(k2) = Btmp % Values(k2) - val * MortarDiag
                   ELSE
                     IF( SumThis) SumCount(row) = SumCount(row) + 1
                   END IF
                 END DO
               END IF
             END IF diag
           END IF

           IF( AllocationsDone ) THEN
             IF( IntegralBC ) THEN
               Btmp % Rhs(row) = SetVal(1)
             ELSE IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
                 Btmp % Rhs(row) = Btmp % Rhs(row) + wsum * MortarBC % rhs(i)
               END IF
             END IF

             ! If every component is uniquely summed we can compute the row indexes simply
             IF( .NOT. SumThis ) THEN
               Btmp % Rows(row+1) = k2 + 1
             END IF
           END IF
         END DO
         
       ELSE IF( ComplexSumRow ) THEN

         IF(IsDG) CALL Fatal(Caller,'DG not implemented for complex systems!')
         
         CALL Info(Caller,'Using simplified complex summing!',Level=8)
         ComplexSumRow = .TRUE.
         
         ! In case of a vector valued problem create a projector that acts on all 
         ! components of the vector. Otherwise follow the same logic.
         IF( SumThis ) THEN
           DO i=1,Atmp % NumberOfRows                        
             
             IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
               k = Atmp % InvPerm(i)
               IF( k == 0 ) CYCLE
             ELSE
               k = i
             END IF
             
             kk = k
             IF( Reorder ) THEN
               kk = Perm(k)
               IF( kk == 0 ) CYCLE
             END IF
             
             NewRow = ( SumPerm(kk) == 0 )
             IF( NewRow ) THEN
               sumrow = sumrow + 1                
               SumPerm(kk) = sumrow 
             ELSE IF(.NOT. AllocationsDone ) THEN
               EliminatedRows = EliminatedRows + 1
             END IF
           END DO
         END IF
           
         
         DO i=1,Atmp % NumberOfRows           
           
           IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
             k = Atmp % InvPerm(i)
             IF( k == 0 ) CYCLE
           ELSE
             k = i
           END IF
            
           kk = k
           IF( Reorder ) THEN
             kk = Perm(k) 
             IF( kk == 0 ) CYCLE
           END IF
             
           IF( SumThis ) THEN             
             row = SumPerm(kk)
           ELSE
             sumrow = sumrow + 1
             row = sumrow
           END IF

           ! For complex matrices 
           IF( AllocationsDone ) THEN
             Btmp % InvPerm(2*row-1) = rowoffset + 2*(kk-1)+1
             Btmp % InvPerm(2*row) = rowoffset + 2*kk
           END IF

           wsum = 0.0_dp
                        

           DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1
             
             col = Atmp % Cols(l) 
             val = Atmp % Values(l)
             
             IF( Reorder ) THEN
               col2 = Perm(col)
               IF( col2 == 0 ) CYCLE
             ELSE
               col2 = col
             END IF
               
             IF( AllocationsDone ) THEN
               ! By Default there is no scaling
               Scale = 1.0_dp
               IF( ThisIsMortar ) THEN
                 IF( ASSOCIATED( MortarBC % Perm ) ) THEN
                   ! Look if the component refers to the slave
                   IF( MortarBC % Perm( col ) > 0 ) THEN
                     Scale = MortarBC % SlaveScale 
                     wsum = wsum + val
                   ELSE
                     Scale = MortarBC % MasterScale
                   END IF
                 ELSE
                   wsum = wsum + val
                 END IF
                 
                 ! If we sum up to anti-periodic dof then use different sign
                 ! - except if the target is also antiperiodic.
                 IF( PerFlipActive ) THEN
                   IF(  PerFlip(col) .NEQV. PerFlip(k) ) Scale = -Scale
                 END IF
                 
               END IF

               ! Add a new column index to the summed up row               
               ! At the first sweep we need to find the first unset position
               ! Real part
               IF( SumThis ) THEN
                 k2 = Btmp % Rows(2*row-1)
                 DO WHILE( Btmp % Cols(k2) > 0 )
                   k2 = k2 + 1
                 END DO
               ELSE
                 k2 = k2 + 1
               END IF
                                            
               Btmp % Cols(k2) = 2 * col2 - 1
               Btmp % Values(k2) = Scale * val

               k2 = k2 + 1
               Btmp % Cols(k2) = 2 * col2
               Btmp % Values(k2) = 0.0

               ! Complex part
               IF( SumThis ) THEN
                 k2 = Btmp % Rows(2*row)
                 DO WHILE( Btmp % Cols(k2) > 0 )
                   k2 = k2 + 1
                 END DO
               ELSE
                 k2 = k2 + 1
               END IF

               Btmp % Cols(k2) = 2 * col2 - 1 
               Btmp % Values(k2) = 0.0
             
               k2 = k2 + 1
               Btmp % Cols(k2) = 2 * col2 
               Btmp % Values(k2) = Scale * val
             ELSE
               k2 = k2 + 4
               IF( SumThis ) THEN
                 SumCount(2*row-1) = SumCount(2*row-1) + 2
                 SumCount(2*row) = SumCount(2*row) + 2
               END IF
             END IF
           END DO
           
           IF( AllocationsDone ) THEN
             IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
                 Btmp % Rhs(2*row-1) = Btmp % Rhs(2*row-1) + wsum * MortarBC % rhs(i)
               END IF
             END IF
           END IF
         END DO
         
       ELSE

         IF(IsDG) CALL Fatal(Caller,'DG not implemented for vector systems!')

         
         ! dofs > 1
         ! In case of a vector valued problem create a projector that acts on all 
         ! components of the vector. Otherwise follow the same logic.
         DO i=1,Atmp % NumberOfRows           
           DO j=1,cDofs
             
             IF( .NOT. ActiveComponents(j) ) THEN
               CALL Info(Caller,'Skipping component: '//I2S(j),Level=12)
               CYCLE
             END IF
             
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

             IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Active ) ) THEN
                 IF( .NOT. MortarBC % Active(cDofs*(i-1)+j) ) CYCLE
               END IF
             END IF

             IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
               k = Atmp % InvPerm(i)
               IF( k == 0 ) CYCLE
             ELSE
               k = i
             END IF

             kk = k
             IF( Reorder ) THEN
               kk = Perm(k)
               IF( kk == 0 ) CYCLE
             END IF

             IF( SumThis ) THEN
               IF( cDofs*(k-1)+j > SIZE(SumPerm) ) THEN
                 CALL Fatal(Caller,'Index out of range')
               END IF
               NewRow = ( SumPerm(cDofs*(kk-1)+j) == 0 )
               IF( NewRow ) THEN
                 sumrow = sumrow + 1                
                 IF( Priority /= 0 ) THEN
                   ! Use negative sign to show that this has already been set by priority
                   SumPerm(cDofs*(kk-1)+j) = -sumrow 
                 ELSE
                   SumPerm(cDofs*(kk-1)+j) = sumrow 
                 END IF
               ELSE IF( Priority /= PrevPriority .AND. SumPerm(cDofs*(kk-1)+j) < 0 ) THEN
                 IF(.NOT. AllocationsDone ) THEN
                   NeglectedRows = NeglectedRows + 1
                 END IF                 
                 CYCLE
               ELSE
                 IF(.NOT. AllocationsDone ) THEN
                   EliminatedRows = EliminatedRows + 1
                 END IF
               END IF
               row = ABS( SumPerm(cDofs*(kk-1)+j) )
             ELSE
               sumrow = sumrow + 1
               row = sumrow
             END IF

             IF( AllocationsDone ) THEN
               Btmp % InvPerm(row) = rowoffset + Dofs * ( kk - 1 ) + j
             END IF

             
             wsum = 0.0_dp

             DO k=Atmp % Rows(i),Atmp % Rows(i+1)-1             

               col = Atmp % Cols(k)                

               IF( Reorder ) THEN                 
                 IF( col <= permsize ) THEN
                   col2 = Perm(col)
                   IF( col2 == 0 ) CYCLE
                 ELSE 
                   PRINT *,'col too large',col,permsize
                   CYCLE
                 END IF
               ELSE
                 col2 = col
               END IF

                 
               k2 = k2 + 1
               
               IF( AllocationsDone ) THEN
                 Scale = 1.0_dp
                 IF( ThisIsMortar ) THEN
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
                   END IF

                   ! If we sum up to anti-periodic dof then use different sign
                   ! - except if the target is also antiperiodic.
                   IF( PerFlipActive ) THEN
                     IF( PerFlip(col) .NEQV. PerFlip(k) ) Scale = -Scale
                   END IF

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
                 IF( SumThis ) THEN
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

                 IF( Reorder ) THEN
                   IF( col <= permsize ) THEN
                     col2 = Perm(col)
                     IF( col2 == 0 ) CYCLE
                   END IF
                 ELSE
                   col2 = col
                 END IF

                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   Btmp % Cols(k2) = Dofs * ( col2 - 1) + j2
                 ELSE
                   IF( SumThis ) THEN
                     SumCount(row) = SumCount(row) + 1
                   END IF
                 END IF
               END DO

               IF( CreateSelf ) THEN
                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   Btmp % Cols(k2) = Dofs * ( Perm( Atmp % InvPerm(i) ) -1 ) + j2
                 END IF
               END IF
             END IF


             IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag ) THEN
                 IF( .NOT. HaveMortarDiag ) THEN
                   MortarDiag = MortarBC % Diag(cDofs*(i-1)+j)
                   LumpedDiag = MortarBC % LumpedDiag
                 END IF

                 IF( LumpedDiag ) THEN
                   k2 = k2 + 1
                   IF( AllocationsDone ) THEN
                     Btmp % Cols(k2) = row + arows
                     Btmp % Values(k2) = -wsum * MortarDiag
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
                       Btmp % Values(k2) = -Atmp % Values(k) * MortarDiag
                     END IF
                   END DO
                 END IF
               END IF
             END IF

               
             IF( AllocationsDone ) THEN
               IF( IntegralBC ) THEN
                 Btmp % Rhs(row) = SetVal(j)
               ELSE IF( ThisIsMortar ) THEN
                 IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
                   Btmp % Rhs(row) = wsum * MortarBC % rhs(cDofs*(i-1)+j)
                 END IF
               END IF
               IF(.NOT. SumThis ) THEN
                 Btmp % Rows(row+1) = k2 + 1
               END IF
             END IF

           END DO
         END DO
       END IF ! dofs > 1
       
       IF( .NOT. SumThis ) THEN
         rowoffset = rowoffset + Arows
         IF( SumProjectors ) THEN
           CALL Info(Caller,'Not summed up size is ' &
           //I2S(sumrow)//' rows and '//I2S(k2)//' nonzeros',Level=8)
           sumrow0 = sumrow
           k20 = k2
         END IF
       END IF
         
       PrevPriority = Priority 
     END DO ! constrain_ind

     IF( k2 == 0 ) THEN
       CALL Info(Caller,'No entries in constraint matrix!',Level=8)
!      Solver % Matrix % ConstraintMatrix => NULL()
       RETURN
     END IF

     ! Allocate the united matrix of all the boundary matrices
     !-------------------------------------------------------
     IF( .NOT. AllocationsDone ) THEN
       CALL Info(Caller,'Allocating '//&
           I2S(sumrow)//' rows and '//I2S(k2)//' nonzeros',&
           Level=8)

       IF( ComplexSumRow ) THEN
         sumrow = 2 * sumrow
       END IF
       
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
         Btmp % Rows(sumrow0+1) = k20+1 
         DO i=sumrow0+2,sumrow+1
           Btmp % Rows(i) = Btmp % Rows(i-1) + SumCount(i-1)
         END DO
         SumPerm = 0
         DEALLOCATE( SumCount ) 
       END IF

       AllocationsDone = .TRUE.

       GOTO 100
     END IF
     
     CALL Info(Caller,'Used '//I2S(sumrow)//&
         ' rows and '//I2S(k2)//' nonzeros',Level=7)
          
     ! Eliminate entries
     IF( SumProjectors ) THEN
       CALL Info(Caller,'Number of eliminated rows: '//I2S(EliminatedRows),Level=6)
       IF( EliminatedRows > 0 ) CALL CRS_PackMatrix( Btmp ) 
     END IF

     IF( NeglectedRows > 0 ) THEN
       CALL Info(Caller,'Number of neglected rows: '//I2S(NeglectedRows),Level=6)
     END IF
        
     IF( InfoActive(30) ) THEN
       BLOCK
         REAL(KIND=dp), POINTER :: px(:)
         px => Btmp % Values
         CALL VectorValuesRange(px,SIZE(px),'ConstraintMatrix')
       END BLOCK
     END IF

     Solver % Matrix % ConstraintMatrix => Btmp     
     Solver % MortarBCsChanged = .FALSE.

     IF( InfoActive(20) ) THEN
       WRITE(Message,'(A,ES12.3)') 'Sum of constraint matrix entries: ',SUM(Btmp % Values)
       CALL Info(Caller,Message)
       WRITE(Message,'(A,ES12.3)') 'Sum of constraint matrix rhs: ',SUM(Btmp % Rhs)
       CALL Info(Caller,Message)
       CALL Info(Caller,'Constraint matrix cols min:'//I2S(MINVAL(Btmp%Cols)))
       CALL Info(Caller,'Constraint matrix cols max:'//I2S(MAXVAL(Btmp%Cols)))
       CALL Info(Caller,'Constraint matrix rows min:'//I2S(MINVAL(Btmp%Rows)))
       CALL Info(Caller,'Constraint matrix rows max:'//I2S(MAXVAL(Btmp%Rows)))
     END IF

     ! For contact mechanics the number of lagrange multipliers may change.
     ! Hence redistribute the old values to the new initial guess using the InvPerm
     ! to identify the correct location.
     !---------------------------------------------------------------------------------     
     IF ( ListGetLogical( Solver % Values,'Apply Contact BCs', Found ) .AND. &
         ALLOCATED( PrevInvPerm ) ) THEN
       MultName = LagrangeMultiplierName( Solver, SetUnfound = .TRUE. ) 
       Var => VariableGet(Solver % Mesh % Variables, MultName)
       IF( ASSOCIATED( Var ) ) THEN
         ALLOCATE( PrevValues( SIZE( Var % Values ) ) )
         PrevValues = Var % Values
         
         k = 0
         l = SIZE(Btmp % InvPerm) 
         Var % Values = 0.0_dp
         
         DO i=1,l
           DO j=1,SIZE(PrevInvPerm)
             IF( Btmp % InvPerm(i) == PrevInvPerm(j) ) THEN
               k = k + 1
               Var % Values(i) = PrevValues(j)
               EXIT
             END IF
           END DO
         END DO

         CALL Info(Caller,'Previous Lagrange multipliers utilized: '&
             //I2S(k)//' out of '//I2S(l),Level=8)
       END IF
     END IF

     CALL Info(Caller,'Finished creating constraint matrix',Level=12)

   END SUBROUTINE GenerateConstraintMatrix
     

   SUBROUTINE ReleaseConstraintMatrix(Solver) 
     TYPE(Solver_t) :: Solver

     CALL FreeMatrix(Solver % Matrix % ConstraintMatrix)
     Solver % Matrix % ConstraintMatrix => NULL()

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



   ! This routine changes the IP field to DG field just while the results are being written.
   !---------------------------------------------------------------------------------------
   SUBROUTINE Ip2DgSwapper( Mesh, FromVar, ToVar, ToType, ToName )

     TYPE( Mesh_t), POINTER :: Mesh
     TYPE( Variable_t), POINTER :: FromVar
     TYPE( Variable_t), POINTER, OPTIONAL :: ToVar
     INTEGER, OPTIONAL :: ToType
     CHARACTER(*), OPTIONAL :: ToName
       
     TYPE( Variable_t), POINTER :: TmpVar
     INTEGER :: TmpType
     INTEGER :: permsize,ipsize,varsize,i,j,k,n,m,e,t,allocstat,dofs
     TYPE(Element_t), POINTER :: Element
     REAL(KIND=dp) :: fip(32),fdg(32)
     LOGICAL :: DgField, NodalField, ElemField, Listed
     INTEGER, ALLOCATABLE :: NodeHits(:)
     INTEGER, POINTER :: Indexes(:)
     CHARACTER(:), ALLOCATABLE :: TmpName
     CHARACTER(*), PARAMETER :: Caller = 'Ip2DgSwapper'

     INTERFACE 
       SUBROUTINE Ip2DgFieldInElement( Mesh, Parent, nip, fip, np, fdg )
         USE Types
         IMPLICIT NONE
         TYPE(Mesh_t), POINTER :: Mesh
         TYPE(Element_t), POINTER :: Parent
         INTEGER :: nip, np
         REAL(KIND=dp) :: fip(:), fdg(:)
       END SUBROUTINE Ip2DgFieldInElement
     END INTERFACE
       
     IF( FromVar % TYPE /= Variable_on_gauss_points ) THEN
       CALL Warn(Caller,'Only IP fields can be swapped!: '//TRIM(FromVar % Name))
       RETURN
     END IF
       
     CALL Info(Caller,'Swapping variable from ip field: '//TRIM(FromVar % Name),Level=8)

     TmpType = Variable_on_nodes_on_elements
     IF( PRESENT( ToType ) ) THEN
       TmpType = ToType
     ELSE IF( PRESENT( ToVar ) ) THEN
       IF( ASSOCIATED( ToVar ) ) THEN
         TmpType = ToVar % TYPE
       END IF
     END IF
     
     DgField = ( TmpType == Variable_on_nodes_on_elements )
     NodalField = ( TmpType == Variable_on_nodes )
     ElemField = ( TmpType == Variable_on_elements )

     IF(.NOT. (DgField .OR. NodalField .OR. ElemField ) ) THEN
       CALL Fatal(Caller,'Wrong type of variable: '//I2S(TmpType))
     END IF
     
     IF( PRESENT( ToName ) ) THEN
       TmpName = TRIM(ToName)
     ELSE IF( DgField ) THEN
       TmpName = TRIM(FromVar % Name)//'_dg'
     ELSE IF( NodalField ) THEN
       TmpName = TRIM(FromVar % Name)//'_nodal'
     ELSE IF( ElemField ) THEN
       TmpName = TRIM(FromVar % Name)//'_elem'
     END IF
     CALL Info(Caller,'Projected variable is named: '//TRIM(TmpName),Level=20)

          
     TmpVar => NULL()
     IF( PRESENT( ToVar ) ) THEN
       TmpVar => ToVar
     ELSE 
       TmpVar => VariableGet( Mesh % Variables, TmpName )       
     END IF

     dofs = FromVar % Dofs
     Listed = ASSOCIATED( TmpVar ) 
     IF(.NOT. Listed ) THEN
       CALL Info(Caller,'Allocating temporal variable for projection',Level=20)
       ALLOCATE( TmpVar )       
       ! Inherit stuff from the primary field to temporal field
       TmpVar % Name = TmpName 
       TmpVar % Dofs = dofs
       TmpVar % Type = TmpType
       TmpVar % NameLen = LEN_TRIM(TmpVar % Name)
       TmpVar % Solver => FromVar % Solver
     END IF

     ! Calculate the sizes related to the primary variable
     n = Mesh % NumberOfBulkElements
     ipsize = FromVar % Perm(n+1) - FromVar % Perm(1)
     CALL Info(Caller,'Size of ip table: '//I2S(ipsize),Level=20)

     IF( DgField ) THEN
       permsize = 0
       DO t=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(t)
         n = Element % TYPE % NumberOfNodes           
         permsize = permsize + n            
       END DO
     ELSE IF( ElemField ) THEN
       permsize = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
     ELSE
       permsize = Mesh % NumberOfNodes                  
     END IF
     CALL Info(Caller,'Size of permutation table: '//I2S(permsize),Level=20)


     IF(ASSOCIATED( TmpVar % Perm ) ) THEN
       IF( SIZE( TmpVar % Perm ) < permsize ) DEALLOCATE( TmpVar % Perm )
     END IF
     IF(.NOT. ASSOCIATED( TmpVar % Perm ) ) THEN
       ALLOCATE( TmpVar % Perm(permsize), STAT=allocstat )
       IF( allocstat /= 0 ) CALL Fatal(Caller,'Allocation error for TmpVar % Perm')
     END IF
      
     ! Mark the existing permutations in the temporal variable
     TmpVar % Perm = 0 
     DO t=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(t)
       n = Element % Type % NumberOfNodes                     
       e = Element % ElementIndex
       m = FromVar % Perm(e+1) - FromVar % Perm(e)

       IF( m > 0 ) THEN
         IF( DgField ) THEN
           IF( .NOT. ASSOCIATED( Element % DGIndexes ) ) THEN
             CALL Warn(Caller,'Cannot project to DG field without DGIndexes!')
             EXIT
           END IF
           TmpVar % Perm( Element % DgIndexes ) = 1
         ELSE IF( ElemField ) THEN
           TmpVar % Perm( Element % ElementIndex ) = 1
         ELSE
           TmpVar % Perm( Element % NodeIndexes ) = 1
         END IF
       END IF
     END DO

     ! Number the permutations in the temporal variable
     j = 0
     DO i = 1, permsize
       IF( TmpVar % Perm( i ) == 0 ) CYCLE
       j = j + 1
       TmpVar % Perm( i ) = j
     END DO
     varsize = j
     CALL Info(Caller,'Size of target variable: '//I2S(varsize),Level=20)

     
     IF(ASSOCIATED( TmpVar % Values ) ) THEN
       IF( SIZE( TmpVar % Values ) < varsize * dofs ) DEALLOCATE( TmpVar % Values )
     END IF
     IF(.NOT. ASSOCIATED( TmpVar % Values ) ) THEN
       ALLOCATE( TmpVar % Values(varsize * dofs), STAT=allocstat)
       IF( allocstat /= 0 ) CALL Fatal(Caller,'Allocation error for TmpVar % Values')      
     END IF
     TmpVar % Values = 0.0_dp
          
     IF( NodalField ) THEN
       ALLOCATE( NodeHits(varsize) )
       NodeHits = 0
     END IF
     
     DO t=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(t)
       n = Element % Type % NumberOfNodes            
       e = Element % ElementIndex

       m = FromVar % Perm(e+1) - FromVar % Perm(e)

       IF( m == 0 ) CYCLE 

       IF( DgField ) THEN
         IF( ALL( TmpVar % Perm( Element % DgIndexes ) == 0 ) ) CYCLE          
       ELSE IF( ElemField ) THEN
         IF( TmpVar % Perm( t ) == 0 ) CYCLE
       ELSE
         IF( ALL( TmpVar % Perm( Element % NodeIndexes ) == 0 ) ) CYCLE                  
       END IF

       DO k=1,dofs
         DO i=1,m        
           j = FromVar % Perm(t) + i 
           fip(i) = FromVar % Values(dofs*(j-1)+k)
         END DO

         IF( ElemField ) THEN
           j = TmpVar % Perm( t ) 
           TmpVar % Values(dofs*(j-1)+k) = SUM( fip(1:m) ) / m 
           CYCLE
         END IF

         ! Solve the elemental equation involving mass matrix
         CALL Ip2DgFieldInElement( Mesh, Element, m, fip, n, fdg )

         IF( DgField ) THEN
           Indexes => Element % DgIndexes
         ELSE IF( NodalField ) THEN
           Indexes => Element % NodeIndexes
         END IF

         DO i=1,n        
           j = TmpVar % Perm( Indexes(i) ) 
           IF( j > 0 ) THEN
             IF( DgField ) THEN
               TmpVar % Values(dofs*(j-1)+k) = fdg(i)
             ELSE
               TmpVar % Values(dofs*(j-1)+k) = TmpVar % Values(dofs*(j-1)+k) + fdg(i)
               IF( k==1 ) NodeHits(j) = NodeHits(j) + 1
             END IF
           END IF
         END DO
       END DO
     END DO

     IF(.NOT. Listed ) THEN
       IF( PRESENT( ToVar ) ) THEN
         ToVar => TmpVar
       ELSE
         CALL VariableAppend( Mesh % Variables,TmpVar)
       END IF
     END IF
            
     IF( DgField ) THEN
       CALL Info(Caller,'Swapping variable from ip to dg done',Level=12)
     ELSE IF( ElemField ) THEN
       CALL Info(Caller,'Swapping variable from ip to elemental done',Level=12)
     ELSE
       DO k=1,dofs
         WHERE( NodeHits > 0 ) 
           TmpVar % Values(k::dofs) = TmpVar % Values(k::dofs) / NodeHits
         END WHERE
       END DO
       DEALLOCATE( NodeHits ) 
       CALL Info(Caller,'Swapping variable from ip to nodal done',Level=12)
     END IF
     
   END SUBROUTINE Ip2DgSwapper
   !-------------------------------------------------------------------------------------


   ! This routine changes from p-DOFs to higher-order Lagrange DOFs.
   !---------------------------------------------------------------------------------------
   SUBROUTINE p2LagrangeSwapper( Mesh, FromVar, ToVar, LagN, LagPerm, LagSize )

     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Variable_t), POINTER :: FromVar
     TYPE(Variable_t), POINTER :: ToVar
     INTEGER :: LagN
     INTEGER :: LagPerm(:)
     INTEGER :: LagSize

     INTEGER, PARAMETER :: MAX_LAGRANGE_NODES = 729
 
     TYPE(Solver_t), POINTER, SAVE :: PDefsSolver => NULL()
     TYPE(Variable_t), POINTER :: TmpVar
     TYPE(Element_t), POINTER :: Element
     INTEGER :: TmpType
     INTEGER :: varsize,i,j,k,np,nl,t,dofs
     REAL(KIND=dp), ALLOCATABLE :: pSol(:,:), lSol(:,:)
     INTEGER, ALLOCATABLE :: pIndexes(:), lIndexes(:)
     INTEGER, ALLOCATABLE :: NodeHits(:)
     LOGICAL, SAVE :: Visited = .FALSE.
     LOGICAL :: DgField, NodalField, ElemField

     CHARACTER(*), PARAMETER :: Caller = 'p2LagrangeSwapper'

     INTERFACE 
       SUBROUTINE HierarchicPToLagrange(PElement, Degree, PSol, LSol, DOFs, PSolver)
         USE Types
         IMPLICIT NONE
         TYPE(Element_t), POINTER :: PElement 
         INTEGER :: Degree                         
         REAL(KIND=dp) :: PSol(:,:)                
         REAL(KIND=dp) :: LSol(:,:)                
         INTEGER, OPTIONAL :: DOFs                 
         TYPE(Solver_t), POINTER, OPTIONAL :: PSolver 
       END SUBROUTINE HierarchicPToLagrange         
     END INTERFACE

     IF(.NOT. ASSOCIATED(FromVar) ) THEN
       CALL Fatal(Caller,'From variable is not associated!')
     END IF
     
     CALL Info(Caller,'Swapping variable from p field to Lagrange field: '//TRIM(FromVar % Name),Level=8)

     IF (.NOT. Visited) THEN
       ! Pick some p-solver in order to handle special cases
       TmpVar => Mesh % Variables
       DO WHILE(ASSOCIATED(TmpVar))
         IF (TmpVar % Valid) THEN 
           IF (ASSOCIATED(TmpVar % Solver)) THEN
             IF (ALLOCATED(TmpVar % Solver % Def_Dofs)) THEN
               IF (ANY(TmpVar % Solver % Def_Dofs(:,:,6)>0)) THEN
                 PDefsSolver => TmpVar % Solver  
                 EXIT
               END IF
             END IF
           END IF
         END IF
         TmpVar => TmpVar % Next
       END DO
       Visited = .TRUE.
     END IF
     
     ! We can only map p-variables and nodal variables!
     TmpType = FromVar % TYPE
     DgField = ( TmpType == Variable_on_nodes_on_elements )
     ElemField = ( TmpType == Variable_on_elements )
     IF(DgField .OR. ElemField ) THEN
       CALL Warn(Caller,'Wrong type of variable: '//I2S(TmpType))
       RETURN
     END IF
               
     dofs = FromVar % Dofs

     ToVar % Name = FromVar % Name 
     ToVar % Dofs = dofs
     ToVar % TYPE = TmpType
     ToVar % NameLen = FromVar % NameLen
     ToVar % Solver => FromVar % Solver
     ToVar % Perm => NULL()

     IF( ASSOCIATED( ToVar % Values ) ) THEN
       IF( SIZE(ToVar % Values) < dofs * LagSize ) THEN
         DEALLOCATE( ToVar % Values )
       END IF
     END IF
     IF( .NOT. ASSOCIATED( ToVar % Values ) ) THEN
       ALLOCATE( ToVar % Values( dofs * LagSize ) )       
     END IF
     ToVar % Values = 0.0_dp
          
     ALLOCATE( pSol(dofs,MAX_LAGRANGE_NODES), lSol(dofs,MAX_LAGRANGE_NODES), &
         pIndexes(MAX_LAGRANGE_NODES), lIndexes(MAX_LAGRANGE_NODES) ) 
     pSol = 0.0_dp
     lSol = 0.0_dp
          
     DO t=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(t)
       
       IF (ASSOCIATED(FromVar % Solver)) THEN
         np = mGetElementDOFs(pIndexes, Element, FromVar % Solver)
       ELSE
         IF (.NOT. ASSOCIATED(PDefsSolver)) THEN
           np = mGetElementDOFs(pIndexes, Element)
         ELSE
           np = mGetElementDOFs(pIndexes, Element, PDefsSolver)
         END IF
       END IF

       ! If all corner nodes are not active, then there is no possibility for interpolation
       IF( ASSOCIATED(FromVar % Perm)) THEN
         k = Element % TYPE % ElementCode / 100
         IF(k>=5 .AND. k<=7) k=k-1
         IF( ANY(FromVar % Perm(pIndexes(1:k)) == 0) ) CYCLE
       END IF

       ! Copy from a global p-solution to local one (pSol) 
       DO i=1,np
         j = pIndexes(i)
         IF( ASSOCIATED(FromVar % Perm)) j = FromVar % Perm(j)         
         DO k=1,dofs           
           IF(j>0) THEN
             pSol(k,i) = FromVar % Values(dofs*(j-1)+k)
           ELSE
             pSol(k,i) = 0.0_dp
           END IF
         END DO
       END DO       

       IF( ASSOCIATED( FromVar % Solver ) ) THEN      
         CALL HierarchicPToLagrange(Element, LagN, PSol, LSol, DOFs, FromVar % Solver)
       ELSE
         CALL HierarchicPToLagrange(Element, LagN, PSol, LSol, DOFs, PDefsSolver)
       END IF
       
       ! Copy from the local Lagrange solution (lSol) to global one
       nl = GetLagrangeIndexes( Mesh, LagN, Element, lIndexes ) 
       DO i=1,nl
         j = lIndexes(i)
         IF(j<=0 .OR. j>SIZE(LagPerm)) THEN
           PRINT *,'Index error:',i,nl,j,SIZE(LagPerm)
           CYCLE
         END IF
         j = LagPerm(j)
         IF(j==0) CYCLE
         DO k=1,dofs           
           ToVar % Values(dofs*(j-1)+k) = lSol(k,i)
         END DO
       END DO
     END DO

     ! There are just gentle reminders that we could also map discontinuous fields to L-elements
     !IF( NodalField ) THEN
     !  ALLOCATE( NodeHits(varsize) )
     !  NodeHits = 0
     !END IF               
     !DO k=1,dofs
     !  WHERE( NodeHits > 0 ) 
     !    TmpVar % Values(k::dofs) = TmpVar % Values(k::dofs) / NodeHits
     !  END WHERE
     !END DO
     !DEALLOCATE( NodeHits ) 

     CALL Info(Caller,'Swapping variable from p to Lagrange done',Level=12)
     
   END SUBROUTINE P2LagrangeSwapper
   !-------------------------------------------------------------------------------------


   

  
   ! Generic evaluation of field value at given point of element.
   ! The idea is that the field value to be evaluated may be nodal, elemental,
   ! dg, or gauss point field. Perhaps even edge or face element field.
   !-------------------------------------------------------------------------------------
   FUNCTION EvalFieldAtElem( Mesh, Var, Element, Basis, dofi, eigeni, imVal, GotVal ) RESULT ( Val )

     INTERFACE 
       SUBROUTINE Ip2DgFieldInElement( Mesh, Parent, nip, fip, np, fdg )
         USE Types
         IMPLICIT NONE
         TYPE(Mesh_t), POINTER :: Mesh
         TYPE(Element_t), POINTER :: Parent
         INTEGER :: nip, np
         REAL(KIND=dp) :: fip(:), fdg(:)
       END SUBROUTINE Ip2DgFieldInElement
     END INTERFACE

     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Variable_t), POINTER :: Var
     TYPE(Element_t), POINTER :: Element
     REAL(KIND=dp), POINTER :: Basis(:)
     INTEGER, OPTIONAL :: dofi
     INTEGER, OPTIONAL :: eigeni
     LOGICAL, OPTIONAL :: GotVal
     REAL(KIND=dp), OPTIONAL :: ImVal
     REAL(KIND=dp) :: Val

     COMPLEX(KIND=dp) :: cVal
     INTEGER :: VIndexes(100)
     REAL(KIND=dp) :: vValues(100),fValues(100)
     REAL(KIND=dp), POINTER :: rValues(:)
     COMPLEX(KIND=dp), POINTER :: cValues(:)
     INTEGER :: dofs,n,i,j,k,l,i1,i2,nip     
     LOGICAL :: GotIt, IsEigen


     IF( PRESENT(GotVal) ) GotVal = .FALSE.
     GotIt = .FALSE.
     Val = 0.0_dp
     cVal = 0.0_dp

     n = Element % Type % NumberOfNodes

     rValues => Var % Values
     dofs = Var % Dofs
     IF( dofs > 1 ) THEN
       IF(.NOT. PRESENT( dofi ) ) THEN
         CALL Fatal('EvalFieldAtElem','We need "dofi" as argument!')
       END IF
     END IF

     IF( PRESENT( EigenI ) ) THEN       
       IF( dofs > 1 ) THEN
         cValues => Var % EigenVectors(eigeni,dofi :: dofs )
       ELSE
         cValues => Var % EigenVectors(eigeni,:)
       END IF
       IsEigen = .TRUE.
     ELSE
       IF( dofs > 1 ) THEN
         rValues => Var % Values(dofi::dofs)
       ELSE
         rValues => Var % Values
       END IF
       IsEigen = .FALSE.
     END IF

     IF( Var % TYPE == Variable_on_nodes ) THEN
       VIndexes(1:n) = Var % Perm( Element % NodeIndexes ) 

       IF( ALL( VIndexes(1:n) /= 0 ) ) THEN         
         IF( IsEigen ) THEN
           cVal = SUM( Basis(1:n) * cvalues( VIndexes(1:n) ) )
         ELSE
           Val = SUM( Basis(1:n) * rvalues( VIndexes(1:n) ) )
         END IF
         GotIt = .TRUE.
       END IF

     ELSE IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN
       VIndexes(1:n) = Var % Perm( Element % DGIndexes ) 

       IF( ALL( VIndexes(1:n) /= 0 ) ) THEN
         IF( isEigen ) THEN
           cVal = SUM( Basis(1:n) * cValues( VIndexes(1:n) ) )
         ELSE
           Val = SUM( Basis(1:n) * rValues( VIndexes(1:n) ) )
         END IF
         GotIt = .TRUE.
       END IF

     ELSE IF( Var % TYPE == Variable_on_elements ) THEN
       j = Var % Perm( Element % ElementIndex )

       IF( j > 0 ) THEN
         IF( isEigen ) THEN
           cVal = rValues( j ) 
         ELSE
           Val = rValues( j ) 
         END IF
         GotIt = .TRUE.
       END IF

     ELSE IF( Var % Type == Variable_on_gauss_points ) THEN              
       i1 = Var % Perm(Element % ElementIndex)       
       i2 = Var % Perm(Element % ElementIndex+1)-1
       nip = i2-i1
       IF(nip > 0) THEN
         IF( IsEigen ) THEN
           VValues = REAL( cValues(i1:i2-1) )
         ELSE
           VValues(1:nip) = rValues(i1:i2-1)
         END IF
         CALL Ip2DgFieldInElement( Mesh, Element, nip, VValues, n, FValues )         

         Val = SUM( Basis(1:n) * FValues(1:n) )

         ! We don't have complex operator for Ip2Dg yet!
         cVal = Val
         GotIt = .TRUE.
       END IF
     END IF

     IF( PRESENT( GotVal ) ) THEN
       GotVal = GotIt
     ELSE
       CALL Warn('EvalFieldAtElem','Could not evaluate field at element!')
     END IF

     IF( IsEigen ) THEN
       Val = REAL( cVal ) 
       IF( PRESENT( imVal ) ) imVal = AIMAG( cVal ) 
     END IF

     
   END FUNCTION EvalFieldAtElem
   !---------------------------------------------------------------------------------------


   
   ! When we have a transient and time-periodic system it may be
   ! beneficial to store the values and use them as an initial guess
   ! for the next round. This tabulates the system and performs the
   ! initial guess. Should be called first time when coming to new
   ! solver for a given steady-state iteration.
   !-----------------------------------------------------------------------------
   SUBROUTINE StoreCyclicSolution(Solver)
     
     TYPE CyclicSol_t
       REAL(KIND=dp), ALLOCATABLE :: PeriodicSol(:,:), PeriodicMult(:,:), &
           PeriodicNrm(:), PeriodicChange(:), dx(:), dy(:)
       LOGICAL :: DoGuess = .FALSE., PeriodicConv = .FALSE.
       INTEGER :: NConv = 0, N1st = 0       
     END TYPE CyclicSol_t
     
     TYPE(Solver_t), POINTER :: Solver     
     TYPE(Model_t), POINTER :: Model

     REAL(KIND=dp), POINTER :: PeriodicSol(:,:), PeriodicMult(:,:), &
         PeriodicNrm(:), PeriodicChange(:), dx(:), dy(:)
     LOGICAL :: DoGuess, PeriodicConv
     INTEGER :: NConv, N1st       
     
     LOGICAL :: Found
     TYPE(Variable_t), POINTER :: v
     REAL(KIND=dp), POINTER :: x(:), y(:)
     INTEGER :: n, m, Ncycle, Ntime, Nguess, Nstore, GuessMode, LGuessMode, Ntimes, Nslices, mmin
     LOGICAL :: ExportMult, ParallelTime, ParallelSlices
     TYPE(Variable_t), POINTER :: Var
     REAL(KIND=dp) :: Relax, AveErr, AveNrm, Tol
     CHARACTER(*), PARAMETER :: Caller = 'StoreCyclicSolution'
     INTEGER :: nSol     
     TYPE(CyclicSol_t), POINTER :: Sols(:) => NULL(), pSol
     CHARACTER(:), ALLOCATABLE :: MultName
     
     
     SAVE Sols

     CALL Info(Caller,'Saving and restoring cyclic solution!',Level=7)

     Model => CurrentModel
     IF(.NOT. ASSOCIATED( Sols ) ) THEN
       ALLOCATE( Sols(Model % NumberOfSolvers ) )
     END IF
     
     v => VariableGet( Solver % Mesh % Variables, 'coupled iter' )     
     IF( NINT(v % Values(1)) > 1 ) RETURN

     Ncycle = ListGetInteger( Model % Simulation,'Periodic Timesteps',Found)
     IF(.NOT. Found ) THEN
       CALL Fatal(Caller,'"Periodic Timesteps" needed to store the cyclic solution!')
     END IF       
     Ntimes = ListGetInteger( Model % Simulation,'Number of Times',Found )  
     Nslices = ListGetInteger( Model % Simulation,'Number of Slices',Found )

     ParallelSlices = ( Nslices > 1 ) 
     ParallelTime = ( Ntimes > 1 )
     IF( ParallelTime ) Ncycle = Ncycle / Ntimes
     
     v => VariableGet( Solver % Mesh % Variables, 'timestep' )
     Ntime = NINT(v % Values(1))

     ! Nothing to do before solution exists
     IF( Ntime == 1 ) RETURN
     
     x => Solver % Variable % Values 
     IF(.NOT. ASSOCIATED( x ) ) THEN
       CALL Warn(Caller,'Cannot store solution without a solution!')
       RETURN
     END IF

     IF( ASSOCIATED( Solver % Matrix ) ) THEN
       n = Solver % Matrix % NumberOfRows
     ELSE
       n = SIZE( x ) 
     END IF
     
     Relax = ListGetConstReal( Solver % Values,'Parallel Timestepping Relaxation Factor',Found ) 
     IF(.NOT. Found ) Relax = 1.0_dp
     
     GuessMode = ListGetInteger( Solver % Values,'Cyclic Guess Mode',Found ) 
     LGuessMode = ListGetInteger( Solver % Values,'Cyclic Lagrange Guess Mode',Found ) 
     IF(.NOT. Found) LGuessMode = GuessMode

     
     ExportMult = ListGetLogical( Solver % Values, 'Cyclic Lagrange Multiplier', Found )
     IF(.NOT. Found ) THEN
       ExportMult = ListGetLogical( Solver % Values, 'Export Lagrange Multiplier', Found )          
     END IF
       
     IF ( ExportMult ) THEN
       MultName = LagrangeMultiplierName( Solver, SetUnfound = .TRUE. ) 
       Var => VariableGet(Solver % Mesh % Variables, MultName)
       ExportMult = ASSOCIATED( Var ) 
       IF( ExportMult ) THEN
         y => Var % Values
         m = SIZE( y )
       END IF
     END IF
     

     nSol = Solver % SolverId
     pSol => Sols(nsol) 
     
     IF( Ntime == 2 ) THEN       
       ! allocate stuff to save vectors
       IF(.NOT. ALLOCATED( pSol % PeriodicSol ) ) THEN
         CALL Info(Caller,'Allocating for periodic solver values for Solver: '//I2S(nSol),Level=6)
         ALLOCATE( pSol % PeriodicSol(n,Ncycle), pSol % PeriodicNrm(n), pSol % PeriodicChange(n), pSol % dx(n) )
         pSol % PeriodicSol = 0.0_dp
         pSol % PeriodicNrm = 0.0_dp
         pSol % PeriodicChange = 0.0_dp
         
         IF( ExportMult ) THEN
           CALL Info(Caller,'Allocating for periodic Lagrange values of size: '//I2S(m),Level=6)
           ALLOCATE( pSol % PeriodicMult(m,Ncycle), pSol % dy(m) )
           pSol % PeriodicMult = 0.0_dp
         END IF         
       END IF
     END IF

     ! This was added afterwards to support multiple Solvers.
     ! Hence we make the original variables point to the Solver-specific variables.
     PeriodicSol => pSol % PeriodicSol
     PeriodicNrm => pSol % PeriodicNrm
     PeriodicChange => pSol % PeriodicChange     
     dx => pSol % dx
     IF( ExportMult ) THEN
       PeriodicMult => pSol % PeriodicMult
       dy => pSol % dy
     END IF
     
     DoGuess = pSol % doGuess
     PeriodicConv = pSol % PeriodicConv
     N1st = pSol % N1st
     NConv = pSol % NConv
 

     ! Both should be in [1,Ncycle]
     Nstore = MODULO( Ntime-2,Ncycle)+1
     Nguess = MODULO( Ntime-1,Ncycle)+1              

     ! Only add guessing if the values have not been tabulated.
     ! When we return for next "run" then we may use already tabulated values.  
     IF(.NOT. DoGuess .AND. Ntime > Ncycle + 1 ) THEN
       DoGuess = .TRUE.
       N1st = Ntime
     END IF
       
     ! This is the 1st iteration for which initial guess is provided.
     ! If we are a returning customer let's keep the initial guess intact.
     IF( DoGuess ) N1st = MIN( Ntime, N1st )
     
     IF( NGuess == 1 .AND. ParallelTime ) THEN
       CALL Info(Caller,'Performing parallel initial guess!')
       CALL CommunicateCyclicSolution(n,x,Solver % Variable % PrevValues)
       IF( ExportMult ) THEN
         IF(ASSOCIATED( Var % PrevValues ) ) THEN
           CALL Info(Caller,'Performing parallel initial guess for Lagrange variable!')
           CALL CommunicateCyclicSolution(m,y,Var % PrevValues)        
         ELSE
           CALL Warn(Caller,'PrevValues does not exist for Lagrange variable!')           
         END IF
       END  IF
     END IF

     
     ! Perform guess only when there is enough data
     IF( DoGuess ) THEN
       IF( GuessMode == 0 ) THEN
         dx = PeriodicSol(:,Nguess)-PeriodicSol(:,Nstore)
       END IF
       IF( ExportMult ) THEN
         IF( LGuessMode == 0 ) THEN
           dy = PeriodicMult(:,Nguess)-PeriodicMult(:,Nstore)
         END IF
       END IF
     END IF       
     
     PeriodicSol(:,Nstore) = x
     PeriodicNrm(Nstore) = Solver % Variable % Norm
     PeriodicChange(Nstore) = Solver % Variable % NonlinChange
     
     IF( ExportMult ) THEN       
       PeriodicMult(:,Nstore) = y
     END IF
     
     IF( DoGuess ) THEN
       CALL Info(Caller,'Using values from previous cycle for initial guess!')
       IF( ExportMult ) THEN
         CALL Info(Caller,'Initializing Lagrange multipliers of size: '//I2S(SIZE(y)),Level=8)
       END IF

       IF( GuessMode == 0 ) THEN
         x = x + dx
         Solver % Variable % Norm = SQRT( SUM(x(1:n)**2) / n )
       ELSE
         x = PeriodicSol(:,Nguess)
         Solver % Variable % Norm = PeriodicNrm(Nguess)
       END IF

       IF( ExportMult ) THEN
         IF( LGuessMode == 0 ) THEN
           y = y + dy 
         ELSE
           y = PeriodicMult(:,Nguess)
         END IF
       END IF
       
     END IF
     
     ! We have computed one full cycle to deduce we have converged.
     ! After having converged the third one is used for producing the results.
     !------------------------------------------------------------------------
     IF( Ntime > N1st + Ncycle ) THEN
       AveErr = SUM( PeriodicChange ) / Ncycle
       WRITE(Message,'(A,ES12.5)') 'Average cyclic error '//I2S(Ntime)//': ',AveErr
       CALL Info(Caller,Message )

       AveNrm = SUM( PeriodicNrm ) / Ncycle
       WRITE(Message,'(A,ES12.5)') 'Average cyclic norm '//I2S(Ntime)//': ',AveNrm
       CALL Info(Caller,Message )
            
       IF( ParallelTime .OR. ParallelSlices ) THEN
         AveErr = ParallelReduction( AveErr ) / ParEnv % PEs
         WRITE(Message,'(A,ES12.5)') 'Parallel cyclic error '//I2S(Ntime)//': ',AveErr
         CALL Info(Caller,Message )        
         
         AveNrm = ParallelReduction( AveNrm ) / ParEnv % PEs
         WRITE(Message,'(A,ES12.5)') 'Parallel cyclic norm '//I2S(Ntime)//': ',AveNrm
         CALL Info(Caller,Message )        
       END IF
         
       mmin = ListGetInteger( Solver % Values,'Cyclic System Min Iterations',Found)
       Tol = ListGetCReal( Solver % Values,'Cyclic System Convergence Tolerance',Found)

       IF( Found ) THEN
         ! We want to start production from the 1st periodic timestep.
         m = NINT( 1.0_dp * Ntime / Ncycle ) 
         
         IF( Nguess == 1 .AND. AveErr < Tol .AND. m >= mmin ) THEN
           PeriodicConv = .TRUE.
           CALL Info(Caller,'Cyclic convergence reached at step: '//I2S(Ntime),Level=4)         
           m = NINT( 1.0_dp * Ntime / Ncycle ) 
           CALL Info(Caller,'Cyclic convergence reached at cycle: '//I2S(m),Level=4)         
           
           ! Set marker to postprocessing solvers.
           V => VariableGet( Solver % Mesh % Variables, 'Produce' )
           V % Values = 1.0_dp   
         END IF

         ! Update counter if this is a converged solution.
         IF( PeriodicConv ) THEN
           Nconv = Nconv + 1
           IF( Nconv == Ncycle ) THEN
             V => VariableGet( Solver % Mesh % Variables, 'Finish' )
             V % Values = 1.0_dp  
           END IF
         END IF
       END IF
     END IF

     ! Store the flags in Solver-specific container.
     pSol % DoGuess = doGuess
     pSol % PeriodicConv = PeriodicConv
     pSol % N1st = N1st
     pSol % NConv = NConv

     
   CONTAINS

     ! This routine is associated to parallel timestepping.
     ! Here we communicate data among different periodic segments each
     ! of which takes certain interval of the periodic system.
     !--------------------------------------------------------------------------------
     SUBROUTINE CommunicateCyclicSolution(n,x,prevx)
       REAL(KIND=dp), POINTER :: x(:),prevx(:,:)
       INTEGER :: n
       
       INTEGER :: toproc, fromproc
       INTEGER :: mpistat(MPI_STATUS_SIZE), ierr
       REAL(KIND=dp), ALLOCATABLE :: tovals(:), fromvals(:)
       INTEGER :: rank, size, Nslices
       INTEGER :: mpitag
       INTEGER, SAVE :: VisitedTimes = 0

       VisitedTimes = VisitedTimes + 1

       CALL Info(Caller,'Communicating data between time segments!',Level=5)
       
       ! Sent data forward in time.
       ! For multislice model the offset to next/previous partition is bigger. 
       Nslices = ListGetInteger( CurrentModel % Simulation,'Number of Slices',Found )
       IF(.NOT. Found) Nslices = 1
              
       toproc = MODULO( ParEnv % MyPe + Nslices, ParEnv % PEs )
       fromproc = MODULO( ParEnv % MyPe - Nslices, ParEnv % PEs )
            
       ALLOCATE( tovals(n), fromvals(n) )

       !PRINT *,'TimeError'//I2S(ParEnv % Mype)//':',VisitedTimes, SUM(ABS(x-tovals))/SUM(ABS(x))

       tovals = x         
       CALL CheckBuffer( 2*n+n*MPI_BSEND_OVERHEAD )
       
       CALL MPI_BSEND( tovals, n, MPI_DOUBLE_PRECISION, &
           toproc, 2005, MPI_COMM_WORLD, ierr )

       CALL MPI_RECV( fromvals, n, MPI_DOUBLE_PRECISION, &
           fromproc, 2005, MPI_COMM_WORLD, mpistat, ierr )

       prevx(:,1) = fromvals 

       DEALLOCATE( tovals, fromvals ) 

     END SUBROUTINE CommunicateCyclicSolution
   
     
   END SUBROUTINE StoreCyclicSolution
   

   !----------------------------------------------------------------------
   !> This subroutine saves a projector assuming time-periodic system.
   !> There are two operation modes.
   !> a) Fetching a precomputed projector when GotProj argument is provided.
   !> b) Storing a projector when no GotProj argument is provided. 
   !----------------------------------------------------------------------
   SUBROUTINE StoreCyclicProjector(Solver,Proj,GotProj)
     TYPE ProjTable_t
       TYPE(Matrix_t), POINTER :: Proj
     END TYPE ProjTable_t
     TYPE(Solver_t), POINTER :: Solver
     TYPE(Matrix_t), POINTER :: Proj
     LOGICAL, OPTIONAL :: GotProj
     
     TYPE(Variable_t), POINTER :: v
     TYPE(Matrix_t), POINTER :: A     
     TYPE(Model_t), POINTER :: Model
     LOGICAL :: Found
     TYPE(ProjTable_t), POINTER :: ProjTable(:)
     INTEGER :: n, i, Ncycle, Ntime, Nstore, Ntimes
     LOGICAL :: SetProj
     
     SAVE ProjTable
     
     Model => CurrentModel 
     Ncycle = ListGetInteger( Model % Simulation,'Periodic Timesteps')
     Ntimes = ListGetInteger( Model % Simulation,'Number Of Times',Found )
     IF(Found ) Ncycle = Ncycle / Ntimes     

     v => VariableGet( Solver % Mesh % Variables, 'timestep' )
     Ntime = NINT(v % Values(1))

     A => Solver % Matrix
     n = A % NumberOfRows

     ! allocate space for projectors 
     IF(.NOT. ASSOCIATED( ProjTable ) ) THEN
       ALLOCATE( ProjTable(Ncycle) )
       DO i=1,Ncycle
         ProjTable(i) % Proj => NULL()
       END DO
     END IF

     ! Nstrore in [1,Ncycle]
     Nstore = MODULO( Ntime-1,Ncycle)+1

     IF( PRESENT( GotProj ) ) THEN
       ! getting projector
       IF( Ntime <= Ncycle ) THEN
         Proj => NULL()
       ELSE
         Proj => ProjTable(Nstore) % Proj
       END IF
       GotProj = ASSOCIATED( Proj ) 
       IF( InfoActive(20) ) THEN
         PRINT *,'Getting cyclic projector:',GotProj,Ntime,Nstore,Ncycle,ASSOCIATED(Proj)
       END IF
     ELSE
       ! storing projector
       SetProj = .NOT. ASSOCIATED( ProjTable(Nstore) % Proj )       
       IF( SetProj ) ProjTable(Nstore) % Proj => Proj
       IF( InfoActive(20) ) THEN
         PRINT *,'Setting cyclic projector:',SetProj,Ntime,Nstore,Ncycle,ASSOCIATED(Proj)
       END IF
     END IF
         
   END SUBROUTINE StoreCyclicProjector


   ! Set the friction for a vector field in an implicit manner by copying matrix rows of
   ! the normal component to matrix rows of the tangential component multiplied by friction
   ! coefficient and direction vector. If the direction is the velocity/displacement field
   ! this means that the amplitude is treated implicitly while the force direction is
   ! explicit. This should be called after assembling the finite element matrices but
   ! before setting the Dirichlet BCs. 
   !---------------------------------------------------------------------------------------
   SUBROUTINE SetImplicitFriction(Model, Solver, FrictionName, DirectionName )
     TYPE(Model_t) :: Model                       !< Complele model structure
     TYPE(Solver_t) :: Solver                     !< Solver for which to set the friction
     CHARACTER(LEN=*) :: FrictionName             !< Name of the friction coefficient
     CHARACTER(LEN=*), OPTIONAL :: DirectionName  !< Name of the optional direction field that may be used to set the direction
         
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Matrix_t), POINTER :: A
     TYPE(Variable_t), POINTER :: FlowVar     
     TYPE(Nodes_t), SAVE :: Nodes
     TYPE(ValueList_t), POINTER :: BC
     REAL(KIND=dp), POINTER :: Values(:)
     INTEGER, POINTER :: FlowPerm(:)
     LOGICAL, ALLOCATABLE :: NodeDone(:)
     REAL(KIND=dp) :: Coeff, Normal(3), LocalNormal(3), LocalT1(3), LocalT2(3), NTT(3,3)
     TYPE(Element_t), POINTER :: Element
     INTEGER, POINTER :: NodeIndexes(:)
     INTEGER :: i,j,k,k2,k3,l,l2,l3,n,t,nb
     INTEGER :: dofN, dofT1, dofT2, bc_id, dim, dofs
     LOGICAL :: Rotated, Found, ExcludePressure
     REAL(KIND=dp), POINTER :: VeloDir(:,:)
     REAL(KIND=dp) :: VeloCoeff(3),AbsVeloCoeff
     CHARACTER(*), PARAMETER :: Caller = 'SetImplicitFriction'
     INTEGER :: Indexes(100)
     REAL(KIND=dp) :: NodeCoeff(27)
     
     IF(.NOT. ListCheckPresentAnyBC( Model, FrictionName ) ) RETURN

     CALL Info(Caller,'Setting fluid friction for boundaries on matrix level!',Level=7)

     Mesh => Model % Mesh
     FlowVar => Solver % Variable 
     FlowPerm => FlowVar % Perm
     A => Solver % Matrix        

     dofs = FlowVar % dofs
     dim = CoordinateSystemDimension()
               
     ALLOCATE( NodeDone( SIZE( FlowPerm ) ) )
     NodeDone = .FALSE.

     ExcludePressure = ListGetLogical( Solver % Values,'Implicit Friction Exclude Pressure',Found)
     
     DO t = Mesh % NumberOfBulkElements+1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(t)

       DO bc_id = 1,Model % NumberOfBCs
         IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) EXIT
       END DO
       IF( bc_id > Model % NumberOfBCs ) CYCLE      
       BC => Model % BCs(bc_id) % Values

       IF( .NOT. ListCheckPresent( BC, FrictionName ) ) CYCLE

       ! Ok we have an active element for friction
       Model % CurrentElement => Element
       NodeIndexes => Element % NodeIndexes
       n = Element % TYPE % NumberOfNodes

       nb = mGetElementDOFs( Indexes, Element, Solver )
       
       ! Normal vector may be needed if this is not all normal-tangential nodes
       CALL CopyElementNodesFromMesh( Nodes, Mesh, n, NodeIndexes)
       Normal = NormalVector( Element, Nodes )

       NodeCoeff(1:n) = ListGetReal( BC, FrictionName, n, NodeIndexes, UnfoundFatal = .TRUE.)
       
       DO i = 1, nb ! n
         j = Indexes(i) ! Nodeindexes(i) 

         IF( NodeDone( j ) ) CYCLE
         IF( FlowPerm( j ) == 0 ) CYCLE

         NodeDone( j ) = .TRUE.         

         IF( i <= n ) THEN
           Coeff = NodeCoeff(i)
         ELSE
           ! Ok, this is a little dirty and only works for constants. 
           Coeff = SUM(NodeCoeff(1:n)) / n
         END IF
           
         ! There is no point of setting too small friction coefficient
         IF(ABS(Coeff) < 1.0d-10) CYCLE

         Rotated = GetSolutionRotation(NTT, j )
         IF( Rotated ) THEN
           ! Get the normal-tangential coordinate system
           LocalNormal = NTT(:,1)
           LocalT1 = NTT(:,2)
           LocalT2 = NTT(:,3)
           DofN = 1
           DofT1 = 2
           DofT2 = 3
         ELSE
           ! Use the standard cartesian coordinate system.
           ! We need to know to which direction the normal is associated to.
           LocalNormal = Normal
           DofN = 1
           DO k=2,dim
             IF( ABS( Normal(k) ) > ABS( Normal(DofN) ) ) DofN = k
           END DO
           IF( ABS(Normal(dofN)) < 0.99 ) THEN
             CALL Warn(Caller,'No normal-tangential system for implicit friction!')            
           END IF           
           IF( DofN == 1 ) THEN
             DofT1 = 2
           ELSE
             DofT1 = 1
           END IF
           DofT2 = 6 - DofT1 - DofN           
           CALL TangentDirections(LocalNormal, LocalT1, LocalT2 )
         END IF

         VeloCoeff = 0.0_dp           

         Found = .FALSE.
         IF( PRESENT( DirectionName ) ) THEN
           VeloDir => ListGetConstRealArray( BC, DirectionName, Found)
         END IF
         IF( Found ) THEN
           VeloCoeff(DofT1) = SUM( VeloDir(1:dim,1) * LocalT1 )
           IF( dim == 3 ) VeloCoeff(DofT2) = SUM( VeloDir(1:dim,1) * LocalT2 )
         ELSE
           k = FlowPerm(j)
           VeloCoeff(DofT1) = FlowVar % Values(Dofs*(k-1)+DofT1) 
           IF(dim==3) VeloCoeff(DofT2) = FlowVar % Values(Dofs*(k-1)+DofT2) 
         END IF

         ! Normalize coefficient to unity so that it only represents the direction of the force
         AbsVeloCoeff = SQRT( SUM( VeloCoeff**2 ) )
         IF( AbsVeloCoeff > TINY(AbsVeloCoeff) ) THEN
           VeloCoeff = VeloCoeff / AbsVeloCoeff
         ELSE
           CYCLE
         END IF
         
         ! Add the friction coefficient 
         VeloCoeff = Coeff * VeloCoeff 

         ! Matrix row associated to normal direction
         j = FlowPerm( j ) 
         k = DOFs * (j-1) + DofN 

         ! Matrix row associated to the two tangent directions
         k2 = DOFs * (j-1) + DofT1 
         A % Rhs(k2) = A % Rhs(k2) + VeloCoeff(DofT1) * A % Rhs(k)

         IF( dim == 3 ) THEN
           k3 = DOFs * (j-1) + DofT2
           A % Rhs(k3) = A % Rhs(k3) + VeloCoeff(DofT2) * A % Rhs(k)             
         END IF

         DO l = A % Rows(k),A % Rows(k+1)-1
           IF( Dofs > dim .AND. ExcludePressure ) THEN
             IF( MODULO(A % Cols(l),Dofs) == 0 ) CYCLE
           END IF
           DO l2 = A % Rows(k2), A % Rows(k2+1)-1
             IF( A % Cols(l2) == A % Cols(l) ) EXIT
           END DO
           A % Values(l2) = A % Values(l2) + VeloCoeff(DofT1) * A % Values(l)
                     
           IF( dim == 3 ) THEN
             DO l3 = A % Rows(k3), A % Rows(k3+1)-1
               IF( A % Cols(l3) == A % Cols(l) ) EXIT
             END DO
             A % Values(l3) = A % Values(l3) + VeloCoeff(DofT2) * A % Values(l)
           END IF
           
         END DO
       END DO
     END DO

     n = COUNT( NodeDone ) 
     CALL Info(Caller,'Number of friction nodes set: '//I2S(n),Level=10)

     DEALLOCATE( NodeDone )

   END SUBROUTINE SetImplicitFriction
   
END MODULE SolverUtils

!> \}
