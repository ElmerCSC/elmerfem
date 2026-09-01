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


MODULE SolverBasics

#include "../config.h"

   USE LoadMod
   USE Multigrid
   USE ElementUtils
   USE IpFieldInterface
   USE PElementVisual
   USE LumpingUtils, ONLY : ComponentStokesTheorem, ComponentCoilEnergy, BoundaryWaveFlux, &
       UpdateDependentComponents, ComponentNodalForceReduction
   USE TimeIntegrate
   USE ModelDescription
   USE MeshBasics, ONLY : CommunicateParallelSystemTag, CylinderFit, &
       DisplaceMesh, FindExtremumNodes, FindMeshEdges, GetLagrangeIndexes, &
       IntegralProjector, MakePermUsingMask
   USE MortarUtils, ONLY : PeriodicProjector, SaveProjector
   USE ParallelUtils
   USE ParallelEigenSolve
   USE MatrixAssembly
   USE MatrixScaling
   
   IMPLICIT NONE
   ! Not re-exported: the external procedure itself USEs modules that would
   ! then import its own name (see module IpFieldInterface).
   PRIVATE :: Ip2DgFieldInElement

CONTAINS

!> Initialize matrix structure and vector to zero initial value.
!------------------------------------------------------------------------------
   SUBROUTINE InitializeToZero( A, ForceVector )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: A  !< Matrix to be initialized
     REAL(KIND=dp) :: ForceVector(:)         !< vector to be initialized
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
     TYPE(ValueListEntry_t), POINTER :: ptr
     CHARACTER(LEN=MAX_NAME_LEN) :: Method
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
     ! Read into a fixed length local instead of
     !   Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )
     ! because that assignment needed a hidden length for the ALLOCATABLE
     ! deferred-length result, which gfortran keeps in static, thread shared
     ! storage. This routine runs once per element in the threaded assembly of
     ! every transient case, so it was the most exposed instance of that problem
     ! in the tree, and its failure mode was silent: a thread handed a zero length
     ! would fall through the SELECT below to CASE DEFAULT and integrate that
     ! element with Newmark-Beta instead of the requested scheme. Assigning to a
     ! fixed length local from the stored string needs no temporary, and saves an
     ! allocate and free per element as a side benefit. The CASE labels below
     ! compare blank padded, so they match exactly as they did before.
     Method = ' '
     ptr => ListFind( Solver % Values, 'Timestepping Method', GotIt )
     IF( ASSOCIATED( ptr ) ) THEN
       IF( ptr % Type /= LIST_TYPE_STRING ) THEN
         CALL Fatal('Add1stOrderTime','Invalid list type for: Timestepping Method')
       END IF
       Method = ptr % CValue
     END IF

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
     TYPE(Solver_t) :: Solver
     
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
     TYPE(Matrix_t) :: StiffMatrix  !< Global stiffness matrix.
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
       k = 0
       
       np = mGetElementDOFs(pIndexes,Element)
       np = MIN(np, n)
       DO i=1,np
         j = pIndexes(i)
         IF(j>0 .AND. j<= SIZE(NT % BoundaryReorder)) THEN
           NormalIndexes(i) = NT % BoundaryReorder(j)
           k = k+1
         END IF
       END DO
       
       IF(k>0) CALL RotateMatrix( LocalStiffMatrix, LocalForce, n, dim, NDOFs, &
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
     TYPE(Matrix_t) :: StiffMatrix  !< The global matrix structure
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
     TYPE(Element_t) :: Element
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
!> Populate values for the upper and lower limiters.
!------------------------------------------------------------------------------
   SUBROUTINE PopulateLimiterValues( Solver )
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver
!-----------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(Mesh_t), POINTER :: Mesh    
     TYPE(variable_t), POINTER :: Var 
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k,n,t,ind,dofs,dof, bf, bc, Upper, ElemFirst, ElemLast, totsize
     REAL(KIND=dp), POINTER :: FieldValues(:), ElemLimit(:)
     REAL(KIND=dp) :: val, bigval 
     INTEGER, POINTER :: FieldPerm(:), NodeIndexes(:)
     LOGICAL :: Found,AnyLimitBC, AnyLimitBF 
     TYPE(ValueList_t), POINTER :: Entity
     REAL(KIND=dp), POINTER :: LimitValues(:)
     CHARACTER(:), ALLOCATABLE :: Name, LimitName
     CHARACTER(*), PARAMETER :: Caller = 'PopulateLimiterValues'
     
     Model => CurrentModel
     Var => Solver % Variable
     Mesh => Solver % Mesh
     
     FieldValues => Var % Values
     FieldPerm => Var % Perm
     totsize = SIZE( FieldValues )
     dofs = Var % dofs
     
     n = Mesh % MaxElementNodes
     ALLOCATE( ElemLimit(n) )

     bigval = 1.0_dp / EPSILON(bigval)
     
     ! Loop through upper and lower limits     
     !------------------------------------------------------------------------
     DO Upper=0,1      
       
       ! Go through the components of the field, if many
       !-------------------------------------------------
       DO DOF = 1,dofs
         
         Name = TRIM(Var % name)
         IF ( Var % DOFs > 1 ) name = ComponentName(name,DOF)
         
         ! The keywords for the correct lower or upper limit of the variable
         !------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           LimitName = TRIM(name)//' Lower Limit'           
         ELSE
           LimitName = TRIM(name)//' Upper Limit' 
         END IF
         
         AnyLimitBC = ListCheckPresentAnyBC( Model, LimitName )
         AnyLimitBF = ListCheckPresentAnyBodyForce( Model, LimitName )

         ! If there is no active keyword then there really is nothing to do
         !----------------------------------------------------------------
         IF( .NOT. ( AnyLimitBC .OR. AnyLimitBF ) ) CYCLE
         
         CALL Info(Caller,'Populating limit: '//TRIM(LimitName),Level=10)
         
         ! Define the range of elements for which the limiters are active
         !---------------------------------------------------------------
         ElemFirst = Mesh % NumberOfBulkElements + 1           
         ElemLast = Mesh % NumberOfBulkElements 
         
         IF( AnyLimitBF ) ElemFirst = 1
         IF( AnyLimitBC ) ElemLast = Mesh % NumberOfBulkElements + &
             Mesh % NumberOfBoundaryElements 
         
         ! Check that active set vectors for limiters exist, otherwise allocate
         !---------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           IF(ASSOCIATED(Var % LowerLimit)) THEN
             IF(SIZE(Var % LowerLimit) /= totsize) THEN
               DEALLOCATE(Var % LowerLimit)
             END IF
           END IF
           IF( .NOT. ASSOCIATED(Var % LowerLimit ) ) THEN
             CALL Info(Caller,'Allocating LowerLimit for variable: '//TRIM(Name),Level=10)
             ALLOCATE( Var % LowerLimit( totsize ) )
             Var % LowerLimit = -bigval
           END IF
           LimitValues => Var % LowerLimit
         ELSE
           IF(ASSOCIATED(Var % UpperLimit)) THEN
             IF(SIZE(Var % UpperLimit) /= totsize) THEN
               DEALLOCATE(Var % UpperLimit)
             END IF
           END IF
           IF( .NOT. ASSOCIATED( Var % UpperLimit ) ) THEN
             CALL Info(Caller,'Allocating UpperLimit for variable: '//TRIM(Name),Level=10)
             ALLOCATE( Var % UpperLimit( totsize ) )
             Var % UpperLimit = bigval
           END IF
           LimitValues => Var % UpperLimit
         END IF
 
         ! In the first time set the initial set 
         !----------------------------------------------------------------------
         DO t = ElemFirst, ElemLast
           
           Element => Model % Elements(t)
           Model % CurrentElement => Element

           n = Element % TYPE % NumberOfNodes
           NodeIndexes => Element % NodeIndexes
           
           Found = .FALSE.
           IF( t > Mesh % NumberOfBulkElements ) THEN
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
           
           DO i=1,n
             j = FieldPerm( NodeIndexes(i) )
             IF( j == 0 ) CYCLE
             ind = Dofs * ( j - 1) + Dof
             
             IF(Upper==0) THEN
               LimitValues(ind) = MAX( LimitValues(ind), ElemLimit(i))
             ELSE
               LimitValues(ind) = MIN( LimitValues(ind), ElemLimit(i))
             END IF
           END DO
         END DO
       END DO
     END DO
     
     DEALLOCATE( ElemLimit ) 
     CALL Info(Caller,'All done',Level=12)
          
   END SUBROUTINE PopulateLimiterValues
     

   
!------------------------------------------------------------------------------
!> Determine soft limiters set. This is called after the solution.
!> and can therefore be active only on the 2nd nonlinear iteration round.
!------------------------------------------------------------------------------
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
    TYPE(Element_t) :: Element
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
    TYPE(Matrix_t) :: A   !< The global matrix
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
     INTEGER, INTENT(IN) :: Perm(:)              !< The permutation of the associated variable
     TYPE(Matrix_t), INTENT(INOUT) :: A         !< The coefficient matrix of the problem
     REAL(KIND=dp), INTENT(INOUT) :: F(:)           !< The RHS vector of the problem
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
          .NOT. ASSOCIATED( Solver % Variable ) ) RETURN
     IF ( .NOT. ASSOCIATED( Solver % Variable % Values ) ) RETURN

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
    TYPE(Mesh_t), TARGET  :: TopMesh,PrimaryMesh
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
      IF ( .NOT.ASSOCIATED( Mesh, PrimaryMesh ) ) THEN
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
     INTEGER, POINTER :: PeriodicPerm(:)
     LOGICAL :: Found
!------------------------------------------------------------------------------

     NT => CurrentModel % Solver % NormalTangential
     IF ( NT % NormalTangentialNOFNodes <= 0 ) RETURN

     PeriodicPerm => NULL()
     IF(ListGetLogical(CurrentModel % Solver % Values,'Apply Conforming BCs',Found)) THEN
       PeriodicPerm => CurrentModel % Mesh % PeriodicPerm
     END IF
          
     dim = CoordinateSystemDimension()
     IF ( ndofs < dim ) RETURN
     
     CALL Info('BackRotateNTSystem','Rotating n-t solution back to cartesian coordinates',Level=12)
     
     DO i=1,SIZE(NT % BoundaryReorder)
       k = NT % BoundaryReorder(i)
       IF ( k <= 0 ) CYCLE
       j = Perm(i)
       IF ( j <= 0 ) CYCLE

       IF(ASSOCIATED(PeriodicPerm)) THEN
         IF(PeriodicPerm(i) > 0) CYCLE
       END IF
       
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
    TYPE(NeighbourList_t),POINTER  :: NeighbourList(:) => NULL()
    
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

      IF( ASSOCIATED(Solver % Matrix) ) THEN
        ! Usually the neighbours are available in the parallel matrix. 
        NeighbourList => Solver % Matrix % ParallelInfo % NeighbourList
      ELSE
        ! There are some exceptions when no matrix, and hence no associated
        ! communication has been created and we have to use the communication structure
        ! of the mesh. Note that this is currently limited to scalar fields!
        NeighbourList => Solver % Mesh % ParallelInfo % NeighbourList
      END IF
      
      SELECT CASE(NormDim)

      CASE(0) 
        DO j=1,n
          IF(PassiveDof(MODULO(j-1,Dofs))) CYCLE
          IF( NeighbourList(j) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = MAX( Norm, ABS( val ) )
          totn = totn + 1
        END DO

      CASE(1)
        DO j=1,n
          IF(PassiveDof(MODULO(j-1,Dofs))) CYCLE
          IF( NeighbourList(j) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = Norm + ABS(val)
          totn = totn + 1
        END DO

      CASE(2)          
        DO j=1,n
          IF(PassiveDof(MODULO(j-1,Dofs))) CYCLE
          IF( NeighbourList(j) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = Norm + val**2
          totn = totn + 1
        END DO

      CASE DEFAULT
        DO j=1,n
          IF(PassiveDof(MODULO(j-1,Dofs))) CYCLE
          IF( NeighbourList(j) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
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
    ! Steady-state level is often sufficient if the dependence is on some other solver.
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

    IF(DoIt) THEN
      IF( LIstGetLogical( SolverParams,'Calculate Velocity Skip First Timestep',Found ) ) THEN
        dtVar => VariableGet( Solver % Mesh % Variables, 'timestep' )
        DoIt = (NINT(dtVar % Values(1)) > 1 )
      END IF
    END IF
    
    IF( DoIt ) THEN
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
        CALL Info(Caller,'Updating variable velocity: '//TRIM(str))
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
    REAL(KIND=dp), POINTER :: x0(:), dx(:)
    REAL(KIND=dp) :: Norm, PrevNorm, rNorm, bNorm, Change, PrevChange, Relaxation, tmp(1),dt, &
        Tolerance, MaxNorm, eps, Ctarget, Poffset, nsum, dpsum
    INTEGER, TARGET  ::  Dnodes(1)
    INTEGER, POINTER :: Indexes(:)
    TYPE(Variable_t), POINTER :: iterVar, VeloVar, dtVar, WeightVar
    CHARACTER(:), ALLOCATABLE :: str,SolverName,ConvergenceType
    
    LOGICAL :: Stat, ConvergenceAbsolute, Relax, RelaxBefore, DoIt, Skip, &
        SkipConstraints, ResidualMode, RelativeP, NodalNorm, DoAitkenRelax
    TYPE(Matrix_t), POINTER :: MMatrix
    REAL(KIND=dp), POINTER CONTIG :: Mx(:), Mb(:), Mr(:)
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: TmpXVec, TmpRVec, TmpRHSVec
    INTEGER :: ipar(1)
    TYPE(ValueList_t), POINTER :: SolverParams
    CHARACTER(*), PARAMETER :: Caller = 'ComputeChange'
    LOGICAL :: Parallel, SingleMesh, x0Allocated, LimitRelax
    LOGICAL, ALLOCATABLE :: LimitMask(:)
    
    SolverParams => Solver % Values
    RelativeP = .FALSE.
    SingleMesh = Solver % Mesh % SingleMesh
    LimitRelax = .FALSE.
    DoAitkenRelax = .FALSE.
    
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

      Relaxation = ListGetCReal( SolverParams,'Steady State Relaxation Factor', Stat )
      IF(.NOT. Stat) Relaxation = 1.0_dp
      Relax = ( ABS(Relaxation-1.0_dp) > EPSILON(Relaxation) )

      iterVar => VariableGet( Solver % Mesh % Variables, 'coupled iter', UnfoundFatal=.TRUE.)
      IterNo = NINT( iterVar % Values(1) )

      DoAitkenRelax = ListGetLogical( SolverParams,'Steady State Aitken Relaxation',Stat)
      IF( DoAitkenRelax ) Relax = .TRUE.
        
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
        LimitRelax = ASSOCIATED(Solver % Variable % LowerLimitActive) .OR. &
            ASSOCIATED(Solver % Variable % UpperLimitActive)
      END IF
     
      ! Steady state system has never any constraints
      SkipConstraints = .FALSE.
      
    ELSE
      IterNo = 0
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
              
      Relaxation = ListGetCReal( SolverParams,'Nonlinear System Relaxation Factor', Stat )
      IF(.NOT. Stat) Relaxation = 1.0_dp
      Relax = ( ABS(Relaxation-1.0_dp) > EPSILON(Relaxation) )

      DoAitkenRelax = ListGetLogical( SolverParams,'Nonlinear System Aitken Relaxation',Stat)
      IF( DoAitkenRelax ) Relax = .TRUE.
      
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

    IF( DoAitkenRelax ) THEN
      ! We need to store the suggested change before relaxation because Aitken
      ! relaxated needs the suggested change, not the relaxed one. 
      IF(ASSOCIATED(Solver % Variable % DeltaValues) ) THEN
        IF(SIZE(Solver % Variable % DeltaValues) /= n) THEN
          DEALLOCATE(Solver % Variable % DeltaValues)
        END IF
      END IF
      IF(.NOT. ASSOCIATED(Solver % Variable % DeltaValues )) THEN
        ALLOCATE(Solver % Variable % DeltaValues(n))
      END IF
      dx => Solver % Variable % DeltaValues

      CALL CalculateAitkenRelaxation(IterNo, Relaxation, Solver % AitkenRelax,A,x,x0,dx)

      ! Save current 'dx' and relaxation for the next round estimate.
      Solver % AitkenRelax = Relaxation
      dx = x - x0 
    END IF
    
    x0allocated = .FALSE.
    IF(Stat .AND. .NOT. SkipConstraints ) THEN
      IF (SIZE(x0) /= SIZE(x)) THEN
         CALL Info(Caller,'WARNING: Possible mismatch in length of vectors ('&
            //I2S(SIZE(x))//' vs. '//I2S(SIZE(x0))//')!',Level=10)

         ! Try to account for differing sizes of the x & x0 (effectively "Relaxation Factor=1"
         ! for the additional "x"-vector dofs."ResidualMode" -case unhandled.
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

    
    IF(LimitRelax) THEN
      ! If we do steady-state relaxation then the soft limiters might not be honored.
      ! This is a trial to fix this but still not quite there. Better to do relaxation
      ! on nonlinear system level when having soft limiters. 
      ALLOCATE(LimitMask(n))
      LimitMask = .FALSE.
      IF(ASSOCIATED(Solver % Variable % LowerLimitActive)) &
          LimitMask = Solver % Variable % LowerLimitActive 
      IF(ASSOCIATED(Solver % Variable % UpperLimitActive)) &
          LimitMask = LimitMask .OR. Solver % Variable % UpperLimitActive 
    END IF
    
    IF( ResidualMode ) THEN
      IF(LimitRelax) THEN
        CALL Fatal(Caller,'Residual mode and limited relaxation cannot be combined!')
      END IF
      IF(Relax .AND. RelaxBefore) THEN
        x(1:n) = x0(1:n) + Relaxation*x(1:n)
      ELSE
        x(1:n) = x0(1:n) + x(1:n)
      END IF
    ELSE IF(Relax .AND. RelaxBefore) THEN
      IF(LimitRelax) THEN
        WHERE(.NOT. LimitMask)
          x(1:n) = (1-Relaxation)*x0(1:n) + Relaxation*x(1:n)
        END WHERE          
      ELSE
        x(1:n) = (1-Relaxation)*x0(1:n) + Relaxation*x(1:n)
      END IF
      IF( RelativeP ) x(dofs:n:dofs) = x(dofs:n:dofs) + Poffset
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
      WRITE( Message, '(a,g15.8,g15.8,a,i0)') &
          'Norm and previous norm: ',Norm,PrevNorm,' size: ',n
      CALL Info(Caller,Message,Level=3)
      CALL NumericalError(Caller,'Norm of solution appears to be NaN')
    END IF

    IF( SteadyState ) THEN
      MaxNorm = ListGetCReal( SolverParams, &
          'Steady State Max Norm', Stat )
    ELSE
      MaxNorm = ListGetCReal( SolverParams, &
          'Nonlinear System Max Norm', Stat )
    END IF    

    IF( Stat .AND. Norm > MaxNorm ) THEN
      WRITE( Message, '(a,g15.8,g15.8)') &
          'Computed norm and given max norm: ',Norm,MaxNorm
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
      IF(LimitRelax) THEN
        WHERE(.NOT. LimitMask)
          x(1:n) = (1-Relaxation)*x0(1:n) + Relaxation*x(1:n)
        END WHERE          
      ELSE        
        x(1:n) = (1-Relaxation)*x0(1:n) + Relaxation*x(1:n)
      END IF
      IF( RelativeP ) x(dofs:n:dofs) = x(dofs:n:dofs) + Poffset
      Solver % Variable % Norm = ComputeNorm(Solver,n,x)
    END IF
    IF(LimitRelax) DEALLOCATE(LimitMask)
    
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
    ! May be useful for some floating systems where one want to impose some integral 
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
        iterVar => VariableGet( Solver % Mesh % Variables, 'coupled iter',UnfoundFatal=.TRUE.)
        IterNo = NINT( iterVar % Values(1) )
      ELSE
        iterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter',UnfoundFatal=.TRUE.)
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


    SUBROUTINE CalculateAitkenRelaxation(iter,r,r0,A,x0,x1,dx) 
      INTEGER :: iter
      TYPE(Matrix_t) :: A
      REAL(KIND=dp) :: r,r0,r1,x0(:),x1(:),dx(:)

      INTEGER :: i,j,k,n,dofs,citer,nOwned
      REAL(KIND=dp) :: aa, dd, s, ds
      INTEGER, ALLOCATABLE :: OwnedPerm(:)

      IF(iter <= 1 ) RETURN
            
      citer = ListGetInteger( Solver % Values,'Aitken Relaxation Factor After Iterations', Stat )
      citer = MAX(1,citer)
      
      n = A % NumberOfRows

      k = 0
      dofs = Solver % Variable % dofs
      IF(dofs>1) k = ListGetInteger( Solver % Values,'Aitken Relaxation Component', Stat ) 

      ! In parallel a dof may be shared by several partitions and must be
      ! accounted for only once. Loop over the dofs owned by this partition and
      ! reduce the two sums separately over the partitions. In serial the
      ! permutation is the identity and the same loop applies.
      n = MIN(n,SIZE(dx),SIZE(x0),SIZE(x1))
      nOwned = ParallelOwnedPerm(n,OwnedPerm,Matrix=A,Mesh=Solver % Mesh)

      aa = 0.0_dp
      dd = 0.0_dp
      DO j=1,nOwned
        i = OwnedPerm(j)
        IF( k > 0 ) THEN
          IF( MODULO(i-1,dofs)+1 /= k ) CYCLE
        END IF
        ds = x0(i) - x1(i) - dx(i)
        aa = aa + dx(i) * ds
        dd = dd + ds**2
      END DO

      IF( ParEnv % PEs > 1 ) THEN
        aa = ParallelReduction(aa)
        dd = ParallelReduction(dd)
      END IF

      ! Nothing to extrapolate from, the previous iterate is already at a fixed point.
      IF( dd < TINY(dd) ) THEN
        CALL Info('CalculateAitkenRelaxation',&
            'Vanishing denominator, keeping the previous relaxation factor',Level=8)
        RETURN
      END IF
      
      r1 = - r0 * aa / dd 

      WRITE(Message,'(A,ES12.3)') 'Aitken relaxation factor suggested: ',r1
      CALL Info('CalculateAitkenRelaxation', Message, Level=10) 

      IF(iter > citer ) THEN
        s = ListGetCReal( Solver % Values,'Aitken Relaxation Factor Relaxation', Stat )
        IF(Stat) r1 = s * r1 + (1-s) * r0

        s = ListGetCReal( Solver % Values,'Aitken Relaxation Factor Minimum', Stat )
        IF(Stat) r1 = MAX(s,r1)

        s = ListGetCReal( Solver % Values,'Aitken Relaxation Factor Maximum', Stat ) 
        IF(Stat) r1 = MIN(s,r1)

        s = ListGetCReal( Solver % Values,'Aitken Relaxation Factor Increase', Stat )
        IF( Stat ) THEN
          ! Don't allow for too large increase
          IF( ABS(r1) > s * ABS(r0) ) r1 = SIGN(s*r0,r1) 
        END IF

        s = ListGetCReal( Solver % Values,'Aitken Relaxation Factor Decrease', Stat )
        IF( Stat ) THEN
          ! Don't allow for too large reduction
          IF( ABS(r1) < s * ABS(r0) ) r1 = SIGN(s*r0,r1) 
        END IF
        
        WRITE(Message,'(A,ES12.3)') 'Aitken relaxation factor regulated: ',r1
        CALL Info('CalculateAitkenRelaxation', Message, Level=5) 
       
        r = r1
      END IF
      
    END SUBROUTINE CalculateAitkenRelaxation
      
    
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
    INTEGER :: AdaptOrder, AdaptNp, Np, RelOrder, BaseNp, BaseRelOrder
    REAL(KIND=dp) :: MinLim, MaxLim, MinV, MaxV, V
    LOGICAL :: UseAdapt, Found,ElementalRule
    INTEGER :: i,j,n,ElementalNp(8),prevVisited = -1
    LOGICAL :: Debug, InitDone, pRef, IsBC, prevIsBC, AdaptSplit, UseNameSpace, EdgePRef
    INTEGER :: EdgeBasisDegree
    REAL(KIND=dp) :: ElemPhi(27)
    LOGICAL :: ElemCut(8)
    TYPE(Nodes_t) :: ElemNodes

    ! NOTE: pRef, Np, RelOrder and ElemPhi are deliberately NOT in the SAVE
    ! list. They are element-dependent (written per element in the rule
    ! selection at the bottom of this function) and must be per-call locals so
    ! that concurrent threads each get their own copy.
    !
    ! pRef used to be SAVEd shared state written only inside the init block;
    ! when several threads hit their first element concurrently they all
    ! re-ran the init, which transiently toggled the shared pRef
    ! .FALSE. -> .TRUE. while other threads were already reading it for the
    ! rule selection at the bottom of this function. A p-element could then
    ! rarely be integrated with the non-p rule: slightly wrong local matrix,
    ! or - with too few points for the bubble basis - a singular block in
    ! CondensateP (the intermittent "LUDecomp: Matrix is singular" failure
    ! in e.g. Step_stokes_heat_vec / SD_Step_stokes_heat_vec).
    !
    ! Np and RelOrder are the same hazard class: they hold a solver-level
    ! default (captured once during init) that is then overwritten per element
    ! for the elemental / adaptive integration rules. The persistent default
    ! lives in the SAVEd BaseNp/BaseRelOrder; the per-call Np/RelOrder are
    ! re-seeded from those after the init block on every call. ElemPhi is
    ! pure per-element scratch (nodal values of the adaptive variable).
    !
    ! REMAINING HAZARD (not fixed here): the "Adaptive Integration Split"
    ! sub-branch below still uses SAVEd scratch with shared state - ElemNodes
    ! (POINTER components, lazy-allocated) and, in its inner BLOCK,
    ! PieceElement / IPtmp / PieceNodes. This 2D-only path is exercised for
    ! correctness by the ModelPDEipsplit test, but only serially: its solver
    ! (ModelPDEhandle / AdvDiffSolver) has no !$OMP PARALLEL assembly, so the
    ! SAVEd scratch is never actually contended there. It remains latently
    ! thread-unsafe and would race if a threaded (*Vec-style) solver ever used
    ! Adaptive Integration Split; rework it to per-call scratch before doing
    ! so. ElemNodes is left SAVE here (its POINTER components would leak on
    ! every call as a plain non-SAVE local) pending that rework.
    SAVE prevSolver, UseAdapt, MinLim, MaxLim, IntegVar, AdaptOrder, AdaptNp, BaseRelOrder, BaseNp, &
        ElementalRule, ElementalNp, prevVisited, EdgePRef, prevIsBC, AdaptSplit, ElemNodes

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
      ! Unsynchronized double-checked-lock, same anti-pattern already fixed
      ! elsewhere this session (EnsureStores, Ip2DgFieldInElement): the fast
      ! InitDone check above happens with no lock, so several threads hitting
      ! their first element of a parallel assembly concurrently (no serial
      ! warm-up call exists before e.g. StatElecSolve's BulkAssembly) can all
      ! enter this block together. Unlike the already-fixed pRef/Np/RelOrder
      ! hazard, the danger here is not that the *values* differ across
      ! threads (for a plain SIF with no adaptive/elemental/namespace keys
      ! they are deterministic) -- it is that nothing orders or flushes the
      ! writes to BaseNp/BaseRelOrder/ElementalRule/UseAdapt/MinLim/MaxLim/
      ! EdgePRef/ElementalNp/AdaptNp/AdaptOrder/AdaptSplit relative to the
      ! terminal "prevSolver => pSolver" / "prevVisited = ..." publish below.
      ! A third thread can observe InitDone flip to TRUE (via those two
      ! writes becoming visible) while an earlier field is still stale or
      ! mid-write under compiler/CPU store reordering, then read garbage into
      ! the per-element rule selection just below this block. Re-checking
      ! InitDone inside the critical section restores real double-checked
      ! locking: cheap in the common (already-initialized) case, since the
      ! unlocked fast path above still avoids the lock entirely once done.
      ! NOTE: this MUST stay an unnamed CRITICAL. A named variant
      ! (GaussPointsAdaptInit) deterministically SIGSEGVs several MPI tests
      ! (radiation*/radiator*/mgdyn_torus/pointload _np4) on the Windows
      ! MSYS2/UCRT MinGW gomp runtime -- a toolchain bug in the named-critical
      ! machinery, unrelated to the logic here (removing the name, or the
      ! directive entirely, makes them pass; the unnamed lock behaves).
      !$OMP CRITICAL
      InitDone = ASSOCIATED( pSolver, prevSolver ) .AND. &
          ( prevVisited == pSolver % TimesVisited ) .AND. (.NOT. (IsBC .NEQV. PrevIsBC) )
      IF( .NOT. InitDone ) THEN
      PrevIsBC = IsBC

      IF(IsBC ) THEN
        UseNameSpace = ListCheckPrefix( pSolver % Values,'bc gauss:')         
        IF( UseNameSpace) THEN
          CALL Info('GaussPointsAdapt','Using namespace "bc gauss:" for integration rules!',Level=10)
          CALL ListPushNamespace('bcgauss:')
        END IF
      ELSE
        UseNameSpace = ListCheckPrefix( pSolver % Values,'bulk gauss:')         
        IF( UseNameSpace) THEN
          CALL Info('GaussPointsAdapt','Using namespace "bulk gauss:" for integration rules!',Level=10)
          CALL ListPushNamespace('bulk gauss:')
        END IF
      END IF
                    
      BaseRelOrder = ListGetInteger( pSolver % Values,'Relative Integration Order',Found )
      AdaptNp = 0
      AdaptSplit = .FALSE.
      BaseNp = ListGetInteger( pSolver % Values,'Number of Integration Points',Found )

      ! Elemental explicit rule will dominate over all other rules
      GaussDef = ListGetString( pSolver % Values,'Element Integration Points',ElementalRule )
      IF( ElementalRule ) THEN
        CALL ElementalGaussRules( GaussDef )
      END IF
      
      VarName = ListGetString( pSolver % Values,'Adaptive Integration Variable',UseAdapt )
      IF( UseAdapt ) THEN
        CALL Info('GaussPointsAdapt','Using adaptive gaussian integration rules',Level=10)
        IntegVar => VariableGet( pSolver % Mesh % Variables, VarName )
        IF( .NOT. ASSOCIATED( IntegVar ) ) THEN
          CALL Fatal('GaussPointsAdapt','> Adaptive Integration Variable < does not exist')
        END IF
        IF( IntegVar % TYPE /= Variable_on_nodes ) THEN
          CALL Fatal('GaussPointsAdapt','Wrong type of integration variable!')
        END IF
        AdaptSplit = ListGetLogical( pSolver % Values,'Adaptive Integration Split',Found )          
        IF(AdaptSplit .AND. Element % TYPE % ElementCode > 500 ) THEN
          CALL Warn('GaussPointsAdapt','Adaptive Integration Split only implemented in 2D!')
        END IF
        
        IF( AdaptSplit ) THEN
          MinLim = ListGetCReal( pSolver % Values,'Adaptive Integration Split Limit', Found )
          MaxLim = MinLim
        ELSE
          MinLim = ListGetCReal( pSolver % Values,'Adaptive Integration Lower Limit' )
          MaxLim = ListGetCReal( pSolver % Values,'Adaptive Integration Upper Limit' )
        END IF
        AdaptNp = ListGetInteger( pSolver % Values,'Adaptive Integration Points',Found )
        IF(.NOT. Found ) THEN
          AdaptOrder = ListGetInteger( pSolver % Values,'Adaptive Integration Order',Found )        
        END IF
        IF(.NOT. Found ) AdaptOrder = 1
        !PRINT *,'Adaptive Integration Strategy:',MinLim,MaxLim,AdaptOrder,AdaptSplit,AdaptNp
      END IF

      IF(UseNameSpace) THEN
        IF(IsBC) THEN
          CALL ListPopNamespace('bc gauss:')
        ELSE
          CALL ListPopNamespace('bulk gauss:')
        END IF
      END IF
      
      
      pRef = .FALSE.

      EdgeBasisDegree = 0
      IF( PRESENT(EdgeBasis) ) THEN
        IF( EdgeBasis ) THEN

          CALL EdgeElementStyle(pSolver % Values, pRef, BasisDegree=EdgeBasisDegree)
          ! Store the solver-constant edge-basis reference-element flag for
          ! the per-call pRef selection below (single same-value store; the
          ! local pRef is used within this init block only).
          EdgePRef = pRef

          ! If elemental rule has not been given then use special edge basis rules
          ! to overrule any other rule for the gauss points. 
          IF(.NOT. ElementalRule ) THEN
            ! We can alter between the two explicit rules of edge basis using relative integration order.
            IF( BaseRelOrder /= 0 ) THEN
              IF( BaseRelOrder == 1 .AND. EdgeBasisDegree == 1 ) THEN
                EdgeBasisDegree = 2
              ELSE IF( BaseRelOrder == -1 .AND. EdgeBasisDegree == 2 ) THEN
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
      END IF

      prevSolver => pSolver
      prevVisited = pSolver % TimesVisited
      END IF
      !$OMP END CRITICAL
    END IF

    ! Select the reference-element style for THIS call/element. Computed as a
    ! local on every call - never cached in shared SAVE state - to stay
    ! race-free when concurrent threads re-run the init block above (see the
    ! NOTE at the SAVE statement). isActivePElement is only a flag check plus
    ! a small array scan, so this costs next to nothing per element.
    IF( PRESENT(EdgeBasis) ) THEN
      pRef = .FALSE.
      IF( EdgeBasis ) pRef = EdgePRef
    ELSE IF( PRESENT(pReferenceElement) ) THEN
      pRef = pReferenceElement
    ELSE
      ! Apart from edge elements we may have p-elements defined.
      ! If not specified check from the current solver.
      pRef = isActivePElement(Element,pSolver)
    END IF

    ! Re-seed the per-call rule parameters from the solver-level defaults
    ! captured during init. These may be overwritten per element just below;
    ! keeping them as call-locals (seeded from the SAVEd base each call)
    ! avoids threads clobbering each other's rule selection.
    Np = BaseNp
    RelOrder = BaseRelOrder

    IF( ElementalRule ) THEN
      ! Elemental explicit rule always has the prevalance
      Np = ElementalNp( Element % TYPE % ElementCode / 100 )
    ELSE IF( UseAdapt ) THEN
      RelOrder = 0
      Np = 0

      ! Count number of corner nodes only 
      n = Element % TYPE % ElementCode / 100
      IF(n > 4 .AND. n < 8) n = n-1
      
      BLOCK
        INTEGER :: nn,j1,j2
        REAL(KIND=dp) :: r

        ElemPhi(1:n) = 0.0_dp
        
        DO i=1,n
          j = Element % NodeIndexes(i)
          IF ( j>0 .AND. j<=SIZE(IntegVar % Perm) ) THEN
            j = IntegVar % Perm(j)
            IF ( j>0 ) THEN
              ElemPhi(i) = IntegVar % Values(j)
            END IF
          ELSE IF( ASSOCIATED( Solver % CutInterp ) ) THEN
            nn = SIZE(IntegVar % Perm)
            j1 = IntegVar % Perm(Solver % Mesh % Edges(j-nn) % NodeIndexes(1))
            j2 = IntegVar % Perm(Solver % Mesh % Edges(j-nn) % NodeIndexes(2))
            IF(j1 > 0 .AND. j2 > 0) THEN
              r = Solver % CutInterp(j-nn)
              ElemPhi(i) = r*IntegVar % Values(j1) + (1-r)*IntegVar % Values(j2)
            END IF
          END IF
        END DO
        
        MinV = MINVAL(ElemPhi(1:n))           
        MaxV = MAXVAL(ElemPhi(1:n))
      END BLOCK

        
      IF( .NOT. ( MaxV < MinLim .OR. MinV > MaxLim ) ) THEN

        IF( AdaptSplit ) THEN
          ! THREADING: this adaptive-split branch is SERIAL-ONLY. ElemNodes
          ! (below) and the PieceElement/IPtmp/PieceNodes SAVEd scratch in the
          ! inner BLOCK are shared across calls, so calling this from an OMP
          ! parallel assembly loop would race. Currently safe only because its
          ! solver (ModelPDEhandle/AdvDiffSolver) assembles serially. Rework
          ! this scratch to per-call locals before using adaptive-split with a
          ! threaded (*Vec-style) solver. See the NOTE at the SAVE statement.
          ElemPhi(1:n) = ElemPhi(1:n) - MinLim
          IF(.NOT. ASSOCIATED(ElemNodes % x)) THEN
            ALLOCATE(ElemNodes % x(2*n),ElemNodes % y(2*n), ElemNodes % z(2*n))
          END IF
          ElemNodes % x(1:n) = pSolver % Mesh % Nodes % x(Element % NodeIndexes(1:n) )
          ElemNodes % y(1:n) = pSolver % Mesh % Nodes % y(Element % NodeIndexes(1:n) )
          ElemNodes % z(1:n) = pSolver % Mesh % Nodes % z(Element % NodeIndexes(1:n) )
            
          ElemCut = .FALSE.
          CALL CutSingleElement(Element, ElemNodes, ElemPhi, ElemCut )

          IF(COUNT(ElemCut(1:2*n)) > 1) THEN
            BLOCK
              LOGICAL :: IsCut, IsMore, stat
              INTEGER :: SgnNode, CutCnt, LocalInds(4), m, t
              TYPE(Element_t), TARGET, SAVE :: PieceElement
              TYPE(Element_t), POINTER :: pElement
              REAL(KIND=dp) :: Ssum0, Ssum, u, v, w, x, y, z, Basis(4), detJ
              REAL(KIND=dp), ALLOCATABLE, SAVE :: IPtmp(:,:)
              TYPE( GaussIntegrationPoints_t ) :: IP
              TYPE( Nodes_t ), SAVE :: PieceNodes
              
              IF(.NOT. ASSOCIATED(PieceElement % NodeIndexes)) THEN
                ALLOCATE(PieceElement % NodeIndexes(4),IpTmp(4,20), &
                    PieceNodes % x(4), PieceNodes % y(4), PieceNodes % z(n) )
              END IF
              
              IntegStuff = GaussPoints( Element, PReferenceElement = .FALSE. )
              Ssum0 = SUM(IntegStuff % s(1:IntegStuff % n))

              i = 0
              DO CutCnt=1,10           
                CALL SplitSingleElement(Element, ElemCut, ElemNodes, CutCnt, &
                    IsCut, IsMore, LocalInds, SgnNode )
                IF(.NOT. IsCut) THEN
                  PRINT *,'ElemCut: ',ElemCut(1:n),CutCnt,IsCut,IsMore
                  CALL Warn('GaussPointsAdapt','This should be cut?')
                  RETURN
                END IF
                  
                m = COUNT(LocalInds > 0)
                IF(m<3 .OR. m>4) CALL Fatal('GaussPointsAdapt','This is neither triangle or quad?')
                
                PieceElement % TYPE => GetElementType( 101*m )   
                IP = GaussPoints( PieceElement, PReferenceElement = .FALSE. )

                PieceNodes % x(1:m) = ElemNodes % x(LocalInds(1:m))
                PieceNodes % y(1:m) = ElemNodes % y(LocalInds(1:m))
                PieceNodes % z(1:m) = ElemNodes % z(LocalInds(1:m))
                
                DO t=1,IP % n        
                  i = i+1
                  stat = ElementInfo( PieceElement, PieceNodes, &
                      IP % u(t), IP % v(t), IP % w(t), detJ, Basis )
                  
                  x = SUM(Basis(1:m)*PieceNodes % x(1:m))
                  y = SUM(Basis(1:m)*PieceNodes % y(1:m))
                  z = SUM(Basis(1:m)*PieceNodes % z(1:m))

                  pElement => PieceElement
                  CALL GlobalToLocal( u, v, w, x, y, z, pElement, ElemNodes ) 

                  ! We need temporal space since IP and IntegStuff refer to the same arrays!
                  IpTmp(1,i) = u
                  IpTmp(2,i) = v
                  IpTmp(3,i) = w
                  IpTmp(4,i) = detJ * IP % s(t)
                END DO

                ! No more element needed to split the reference element.
                IF(.NOT. IsMore) EXIT                                
              END DO

              IntegStuff % n = i
              ssum = SUM(IPtmp(4,1:i))

              IntegStuff % u(1:i) = IPtmp(1,1:i) 
              IntegStuff % v(1:i) = IPtmp(2,1:i) 
              IntegStuff % w(1:i) = IPtmp(3,1:i) 
              IntegStuff % s(1:i) = IPtmp(4,1:i) * (ssum0 / ssum )
              
              IF(i>0) THEN
                !PRINT *,'detJ orig:',CutCnt,i,ssum0,ssum,IP % n
              END IF
              !PRINT *,'AdaptiveElement:',Element % ElementIndex, MinV, MaxV, AdaptSplit, IntegStuff % n

              RETURN
              
            END BLOCK

          END IF
        ELSE
          RelOrder = AdaptOrder
          Np = AdaptNp
        END IF
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
            
      IterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter',UnfoundFatal=.TRUE.)
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

    TYPE(Solver_t) :: Solver
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

    TYPE(Solver_t) :: Solver
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

  

!------------------------------------------------------------------------------
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



!-----------------------------------------------------------------------------------
!> Assemble a right-hand-side related to adjoint problem using library functioality.
!-----------------------------------------------------------------------------------
  SUBROUTINE AssembleAdjointRhs( Solver, AdjointRhsName )
    
    TYPE(Solver_t) :: Solver
    CHARACTER(*) :: AdjointRhsName
    !------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Indexes(:)
    INTEGER :: i,j,k,n,t,istat
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: FORCE(:), Load(:), Basis(:)
    REAL(KIND=dp), POINTER :: rhs(:)
    REAL(KIND=dp) :: detJ, val, s
    LOGICAL :: Stat, Found
    TYPE(ValueList_t), POINTER :: vList
    INTEGER :: bf_id, bc_id, body_id, dofs

    CALL Info('AssembleAdjointRhs','Integrating source by name: '//TRIM(AdjointRhsName),Level=10)

    IF(.NOT. ( ListCheckPrefixAnyBC( CurrentModel, AdjointRhsName ) .OR. &
        ListCheckPrefixAnyBodyForce( CurrentModel, AdjointRhsName ) ) ) THEN
      CALL Fatal('AssembleAdjointRhs','Could not find any source with name: '//TRIM(AdjointRhsName))      
    END IF
        
    
    Mesh => Solver % Mesh

    dofs = Solver % Variable % Dofs 
    N = Mesh % MaxElementNodes 
    ALLOCATE( Basis(n), FORCE(dofs*N), Load(dofs*N), &
        Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        STAT=istat)

    rhs => Solver % Matrix % rhsAdjoint
    IF(.NOT. ASSOCIATED(rhs)) THEN
      CALL Fatal('AssembleAdjointRhs','RhsAdjoint not allocated!')
    END IF
    rhs = 0.0_dp
        
    DO t=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 
      Element => Mesh % Elements(t)
      CurrentModel % CurrentElement => Element
      n = Element % TYPE % NumberOfNodes

      vList => NULL()
      body_id = Element % BodyId
      IF( body_id > 0) THEN
        bf_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
            'Body Force', Found )
        IF(bf_id > 0) THEN
          vList => CurrentModel % BodyForces(bf_id) % Values
        END IF
      END IF
      
      IF(.NOT. ASSOCIATED(vList) ) THEN
        IF( t > Mesh % NumberOfBulkElements + 1 ) THEN
          DO bc_id=1,CurrentModel % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) THEN
              vList => CurrentModel % BCs(bc_id) % Values
              EXIT
            END IF
          END DO
        END IF
        IF(.NOT. ASSOCIATED(vList)) CYCLE
      END IF
      
      Indexes => Element % NodeIndexes
      Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
      Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
      Nodes % z(1:n) = Mesh % Nodes % z(Indexes)
      
      Load = 0.0_dp
      FORCE = 0.0_dp
      
      IF(dofs==1) THEN
        Load(1:n) = ListGetReal( vList, AdjointRhsName, n, Indexes, Found ) 
      ELSE
        Found = .FALSE.
        DO i=1,dofs
          Load(i:dofs*n:dofs) = ListGetReal( vList, TRIM(AdjointRhsName)//' '//I2S(i), &
              n, Indexes, Stat )
          Found = Found .OR. Stat          
        END DO        
      END IF
      IF(.NOT. Found) CYCLE
                          
      IP = GaussPoints( Element )

      DO k=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
            IP % W(k),  detJ, Basis )        
        s = IP % s(k) * DetJ 

        DO i=1,dofs
          val = SUM(Basis(1:n) * Load(i:dofs*n:dofs))
          FORCE(i:dofs*n:dofs) = FORCE(i:dofs*n:dofs) + val * s * Basis(1:n) 
        END DO
      END DO
      
      CALL UpdateGlobalForce( Solver % Matrix % RhsAdjoint, &
          Force, n, dofs, Solver % Variable % Perm(Indexes) )
    END DO

    DEALLOCATE( Basis, FORCE, Nodes % x, Nodes % y, Nodes % z)
    
  END SUBROUTINE AssembleAdjointRhs



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
    TYPE(Variable_t), POINTER :: dtVar
    CHARACTER(*), PARAMETER :: Caller = 'FsiCouplingAssembly'
    
    
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
    TYPE(Mesh_t) :: Mesh
    TYPE(Element_t), INTENT(IN) :: Element
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
      TYPE(Element_t) :: Element
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

     TYPE(Mesh_t), TARGET :: Mesh
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


     TYPE(Mesh_t), TARGET :: Mesh
     TYPE(Variable_t), POINTER :: Var
     TYPE(Element_t), TARGET :: Element
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
     IF(.NOT. Found) Relax = 1.0_dp
     
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
    Basis = 0.0_dp
    ElementNodes % x = 0.0_dp
    ElementNodes % y = 0.0_dp
    ElementNodes % z = 0.0_dp
    LocalIndexes = 0
    
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

END MODULE SolverBasics

!> \}
