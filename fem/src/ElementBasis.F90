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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Oct 1996
! *
! ******************************************************************************/

!--------------------------------------------------------------------------------
!>  Module defining element type and operations. The most basic FEM routines
!>  are here, handling the basis functions, global derivatives, etc...
!--------------------------------------------------------------------------------
!> \ingroup ElmerLib
!> \{

#include "../config.h"

!--------------------------------------------------------------------------------
!>  Element basis functions: type setup, polynomial primitives, basis dispatchers,
!>  CrossProduct, GetEdgeMap, ElementDiameter.
!--------------------------------------------------------------------------------
MODULE ElementBasis
   USE Messages
   USE Integration
   USE LinearAlgebra
   USE CoordinateSystems
   USE PElementMaps
   USE PElementBase
   USE H1Basis
   USE Lists
!$ USE omp_lib ! Include module conditionally (for omp_get_max_threads below)

   IMPLICIT NONE

   INTEGER, PARAMETER, PRIVATE :: MaxDeg  = 4, MaxDeg3 = MaxDeg**3, &
                                   MaxDeg2 = MaxDeg**2
   INTEGER, PARAMETER :: MAX_ELEMENT_NODES = 256
   LOGICAL, PRIVATE :: TypeListInitialized = .FALSE.
   TYPE(ElementType_t), POINTER :: ElementTypeList

CONTAINS

    SUBROUTINE SwapRefElemNodes(p)
!------------------------------------------------------------------------------
      LOGICAL :: p
!------------------------------------------------------------------------------
      INTEGER :: n
      TYPE(ElementType_t), POINTER :: et
!------------------------------------------------------------------------------
      
      et => ElementTypeList
      DO WHILE(ASSOCIATED(et))
        n = et % NumberOfNodes

        ! Single node does not really have much options here...
        IF( et % ElementCode < 200 ) THEN
          CONTINUE
        ELSE IF( p .AND. ALLOCATED(et % NodeU) ) THEN
          IF ( .NOT.ALLOCATED(et % P_NodeU) ) THEN
            ALLOCATE(et % P_NodeU(n), et % P_NodeV(n), et % P_NodeW(n))
            CALL GetRefPElementNodes( et,  et % P_NodeU, et % P_NodeV, et % P_NodeW )
          END IF
          et % NodeU = et % P_NodeU
          et % NodeV = et % P_NodeV
          et % NodeW = et % P_NodeW
        ELSE IF ( ALLOCATED(et % N_NodeU) ) THEN
          et % NodeU = et % N_NodeU
          et % NodeV = et % N_NodeV
          et % NodeW = et % N_NodeW
        END IF
        et => et % NextElementType
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE SwapRefElemNodes
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Add an element description to global list of element types.
!------------------------------------------------------------------------------
   SUBROUTINE AddElementDescription( element,BasisTerms )
!------------------------------------------------------------------------------
      INTEGER, DIMENSION(:) :: BasisTerms  !< List of terms in the basis function that should be included for this element type. 
                                           ! BasisTerms(i) is an integer from 1-27 according to the list below.
      TYPE(ElementType_t), TARGET :: element !< Structure holding element type description
!------------------------------------------------------------------------------
!     Local variables
!------------------------------------------------------------------------------
      TYPE(ElementType_t), POINTER :: temp

      INTEGER, DIMENSION(MaxDeg3) :: s
      INTEGER :: i,j,k,l,m,n,upow,vpow,wpow,i1,i2,ii(9),jj

      REAL(KIND=dp) :: u,v,w,r
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: A, B
!------------------------------------------------------------------------------

!     PRINT*,'Adding element type: ', element % ElementCode

      n = element % NumberOfNodes
      element % NumberOfEdges = 0
      element % NumberOfFaces = 0
      element % BasisFunctionDegree = 0
      NULLIFY( element % BasisFunctions )

      IF ( element % ElementCode >= 200 ) THEN

      ALLOCATE( A(n,n) )

!------------------------------------------------------------------------------
!     1D bar elements
!------------------------------------------------------------------------------
      IF ( element % DIMENSION == 1 ) THEN

         DO i = 1,n
           u = element % NodeU(i)
           DO j = 1,n
             k = BasisTerms(j) - 1
             upow = k
             IF ( u==0 .AND. upow == 0 ) THEN
                A(i,j) = 1
             ELSE
                A(i,j) = u**upow
             END IF
             element % BasisFunctionDegree = MAX(element % BasisFunctionDegree,upow) 
           END DO
         END DO

!        ALLOCATE( element % BasisFunctions(MaxDeg,MaxDeg) )

!------------------------------------------------------------------------------
!     2D surface elements
!------------------------------------------------------------------------------
      ELSE IF ( element % DIMENSION == 2 ) THEN

         DO i = 1,n
            u = element % NodeU(i)
            v = element % NodeV(i)
            DO j = 1,n
              k = BasisTerms(j) - 1
              vpow = k / MaxDeg 
              upow = MOD(k,MaxDeg)

              IF ( upow == 0 ) THEN
                 A(i,j) = 1
              ELSE
                 A(i,j) = u**upow
              END IF

              IF ( vpow /= 0 ) THEN
                 A(i,j) = A(i,j) * v**vpow
              END IF

              element % BasisFunctionDegree = MAX(element % BasisFunctionDegree,upow) 
              element % BasisFunctionDegree = MAX(element % BasisFunctionDegree,vpow) 
            END DO
         END DO

!        ALLOCATE( element % BasisFunctions(MaxDeg2,MaxDeg2) )

!------------------------------------------------------------------------------
!     3D volume elements
!------------------------------------------------------------------------------
      ELSE

         DO i = 1,n
            u = element % NodeU(i)
            v = element % NodeV(i)
            w = element % NodeW(i)
            DO j = 1,n
              k = BasisTerms(j) - 1
              upow = MOD( k,MaxDeg )
              wpow = k / MaxDeg2
              vpow = MOD( k / MaxDeg, MaxDeg )

              IF ( upow == 0 ) THEN
                 A(i,j) = 1
              ELSE
                 A(i,j) = u**upow
              END IF

              IF ( vpow /= 0 ) THEN
                 A(i,j) = A(i,j) * v**vpow
              END IF

              IF ( wpow /= 0 ) THEN
                 A(i,j) = A(i,j) * w**wpow
              END IF

              element % BasisFunctionDegree = MAX(element % BasisFunctionDegree,upow) 
              element % BasisFunctionDegree = MAX(element % BasisFunctionDegree,vpow) 
              element % BasisFunctionDegree = MAX(element % BasisFunctionDegree,wpow) 
            END DO
         END DO

!        ALLOCATE( element % BasisFunctions(MaxDeg3,MaxDeg3) )
      END IF

!------------------------------------------------------------------------------
!     Compute the coefficients of the basis function terms
!------------------------------------------------------------------------------
      CALL InvertMatrix( A,n )

      IF ( Element % ElementCode == 202 ) THEN
         ALLOCATE( Element % BasisFunctions(14) )
      ELSE
         ALLOCATE( Element % BasisFunctions(n) )
      END IF

      upow = 0
      vpow = 0
      wpow = 0

      DO i = 1,n
        Element % BasisFunctions(i) % n = n
        ALLOCATE( Element % BasisFunctions(i) % p(n) )
        ALLOCATE( Element % BasisFunctions(i) % q(n) )
        ALLOCATE( Element % BasisFunctions(i) % r(n) )
        ALLOCATE( Element % BasisFunctions(i) % Coeff(n) )

        DO j = 1,n
          k = BasisTerms(j) - 1

          SELECT CASE( Element % DIMENSION ) 
          CASE(1)
             upow = k
          CASE(2)
             vpow = k / MaxDeg 
             upow = MOD(k,MaxDeg)
          CASE(3)
             upow = MOD( k,MaxDeg )
             wpow = k / MaxDeg2
             vpow = MOD( k / MaxDeg, MaxDeg )
           END SELECT

           Element % BasisFunctions(i) % p(j) = upow
           Element % BasisFunctions(i) % q(j) = vpow
           Element % BasisFunctions(i) % r(j) = wpow
           Element % BasisFunctions(i) % Coeff(j) = A(j,i)
        END DO
      END DO

      DEALLOCATE( A )

      IF ( Element % ElementCode == 202 ) THEN
         ALLOCATE( A(14,14) )
         A = 0
         CALL Compute1DPBasis( A,14 )

         DO i=3,14
            ALLOCATE( Element % BasisFunctions(i) % p(i) )
            ALLOCATE( Element % BasisFunctions(i) % q(i) )
            ALLOCATE( Element % BasisFunctions(i) % r(i) )
            ALLOCATE( Element % BasisFunctions(i) % Coeff(i) )

            k = 0
            DO j=1,i
               IF ( A(i,j) /= 0.0d0 ) THEN
                  k = k + 1
                  Element % BasisFunctions(i) % p(k) = j-1
                  Element % BasisFunctions(i) % q(k) = 0
                  Element % BasisFunctions(i) % r(k) = 0
                  Element % BasisFunctions(i) % Coeff(k) = A(i,j)
               END IF
            END DO
            Element % BasisFunctions(i) % n = k
         END DO
         DEALLOCATE( A )
      END IF

!------------------------------------------------------------------------------

      SELECT CASE( Element % ElementCode / 100 )
        CASE(2) 
           Element % NumberOfEdges = 1
        CASE(3) 
           Element % NumberOfFaces = 1
           Element % NumberOfEdges = 3
        CASE(4) 
           Element % NumberOfFaces = 1
           Element % NumberOfEdges = 4
        CASE(5) 
           Element % NumberOfFaces = 4
           Element % NumberOfEdges = 6
        CASE(6) 
           Element % NumberOfFaces = 5
           Element % NumberOfEdges = 8
        CASE(7) 
           Element % NumberOfFaces = 5
           Element % NumberOfEdges = 9
        CASE(8) 
           Element % NumberOfFaces = 6
           Element % NumberOfEdges = 12
      END SELECT

      END IF ! type >= 200

!------------------------------------------------------------------------------
!     And finally add the element description to the global list of types
!------------------------------------------------------------------------------
      IF ( .NOT.TypeListInitialized ) THEN
        ALLOCATE( ElementTypeList )
        ElementTypeList = element
        TypeListInitialized = .TRUE.
        NULLIFY( ElementTypeList % NextElementType )
      ELSE
        ALLOCATE( temp )
        temp = element
        temp % NextElementType => ElementTypeList
        ElementTypeList => temp
      END IF

!------------------------------------------------------------------------------
     ! Allocate reference basis cache on the list node. The last dimension is
     ! per OpenMP thread so the shared per-type cache can be filled lock-free
     ! during threaded assembly (see ElementInfo in ElemInfo.F90).
     n = ElementTypeList % NumberOfNodes
     BLOCK
       INTEGER :: nthr
       nthr = 1
       !$ nthr = omp_get_max_threads()
       ALLOCATE( ElementTypeList % BasisCacheU(ELEM_BASIS_CACHE_SIZE, nthr), &
                 ElementTypeList % BasisCacheV(ELEM_BASIS_CACHE_SIZE, nthr), &
                 ElementTypeList % BasisCacheW(ELEM_BASIS_CACHE_SIZE, nthr), &
                 ElementTypeList % BasisCache(ELEM_BASIS_CACHE_SIZE, n, nthr), &
                 ElementTypeList % dBasisCache(ELEM_BASIS_CACHE_SIZE, n, 3, nthr), &
                 ElementTypeList % BasisCacheCount(nthr) )
       ElementTypeList % BasisCacheCount = 0
     END BLOCK
!------------------------------------------------------------------------------

CONTAINS


!------------------------------------------------------------------------------
!> Subroutine to compute 1D P-basis from Legendre polynomials.
!------------------------------------------------------------------------------
   SUBROUTINE Compute1DPBasis( Basis,n )
!------------------------------------------------------------------------------
     INTEGER :: n
     REAL(KIND=dp) :: Basis(:,:)
!------------------------------------------------------------------------------
     REAL(KIND=dp)   :: s,P(n+1),Q(n),P0(n),P1(n+1)
     INTEGER :: i,j,k,np,info

!------------------------------------------------------------------------------

     IF ( n <= 1 ) THEN
        Basis(1,1)     = 1.0d0
        RETURN
     END IF
!------------------------------------------------------------------------------
! Compute coefficients of n:th Legendre polynomial from the recurrence:
!
! (i+1)P_{i+1}(x) = (2i+1)*x*P_i(x) - i*P_{i-1}(x), P_{0} = 1; P_{1} = x;
!
! CAVEAT: Computed coefficients inaccurate for n > ~15
!------------------------------------------------------------------------------
     P = 0
     P0 = 0
     P1 = 0
     P0(1) = 1
     P1(1) = 1
     P1(2) = 0

     Basis(1,1) =  0.5d0
     Basis(1,2) = -0.5d0

     Basis(2,1) =  0.5d0
     Basis(2,2) =  0.5d0

     DO k=2,n
       IF ( k > 2 ) THEN
          s = SQRT( (2.0d0*(k-1)-1) / 2.0d0 )
          DO j=1,k-1
             Basis(k,k-j+1) = s * P0(j) / (k-j)
             Basis(k,1) = Basis(k,1) - s * P0(j)*(-1)**(j+1) / (k-j)
          END DO
       END IF

       i = k - 1
       P(1:i+1) = (2*i+1) * P1(1:i+1)  / (i+1)
       P(3:i+2) = P(3:i+2) - i*P0(1:i) / (i+1)
       P0(1:i+1) = P1(1:i+1)
       P1(1:i+2) = P(1:i+2)
     END DO
!--------------------------------------------------------------------------
 END SUBROUTINE Compute1DPBasis
!--------------------------------------------------------------------------

   END SUBROUTINE AddElementDescription 
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Read the element description input file and add the element types to a
!>   global list. The file is assumed to be found under the name
!>        $ELMER_HOME/lib/elements.def
!>   This is the first routine the user of the element utilities should call
!>   in his/her code.
!------------------------------------------------------------------------------
   SUBROUTINE InitializeElementDescriptions()
!------------------------------------------------------------------------------
!     Local variables
!------------------------------------------------------------------------------
      CHARACTER(LEN=:), ALLOCATABLE :: tstr, str,elmer_home

      INTEGER :: k, n
      INTEGER, DIMENSION(MaxDeg3) :: BasisTerms

      TYPE(ElementType_t) :: element

      LOGICAL :: gotit, fexist
!------------------------------------------------------------------------------
!     PRINT*,' '
!     PRINT*,'----------------------------------------------'
!     PRINT*,'Reading element definition file: elements.def'
!     PRINT*,'----------------------------------------------'

      !
      ! Add connectivity element types:
      ! -------------------------------
      BasisTerms = 0
      element % GaussPoints  = 0
      element % GaussPoints0 = 0
      element % GaussPoints2 = 0
      element % StabilizationMK = 0
      DO k=3,64
        element % NumberOfNodes = k
        element % ElementCode = 100 + k
        CALL AddElementDescription( element,BasisTerms )
      END DO

      ! then the rest of them....
      !--------------------------
      ALLOCATE(CHARACTER(MAX_PATH_LEN)::elmer_home)

      tstr = 'ELMER_LIB'
      CALL envir( tstr,elmer_home,k ) 
      
      fexist = .FALSE.
      IF (  k > 0 ) THEN
         tstr = elmer_home(1:k) // '/elements.def'
         INQUIRE(FILE=TRIM(tstr), EXIST=fexist)
      END IF
      IF (.NOT. fexist) THEN
        tstr = 'ELMER_HOME'
        CALL envir( tstr,elmer_home,k ) 
        IF ( k > 0 ) THEN
           tstr = elmer_home(1:k)//'/share/elmersolver/lib/elements.def'
           INQUIRE(FILE=TRIM(tstr), EXIST=fexist)
        END IF
        IF ((.NOT. fexist) .AND. k > 0) THEN
           tstr = elmer_home(1:k)//'/elements.def'
           INQUIRE(FILE=TRIM(tstr), EXIST=fexist)
        END IF
     END IF
     IF (.NOT. fexist) THEN
        CALL GetSolverHome(elmer_home, n)
        tstr = elmer_home(1:n)//'/lib/elements.def'
        INQUIRE(FILE=TRIM(tstr), EXIST=fexist)
     END IF
     IF (.NOT. fexist) THEN
        CALL Fatal('InitializeElementDescriptions','elements.def not found')
     END IF

      OPEN( 1,FILE=TRIM(tstr), STATUS='OLD' )

      ALLOCATE(CHARACTER(MAX_STRING_LEN)::str)
      DO WHILE( ReadAndTrim(1,str) )

        IF ( SEQL(str, 'element') ) THEN

          BasisTerms = 0

          gotit = .FALSE.
          DO WHILE( ReadAndTrim(1,str) )

            IF ( SEQL(str, 'dimension') ) THEN
              READ( str(10:), * ) element % DIMENSION

            ELSE IF ( SEQL(str, 'code') ) THEN
              READ( str(5:), * ) element % ElementCode

            ELSE IF ( SEQL(str, 'nodes') ) THEN
              READ( str(6:), * ) element % NumberOfNodes

            ELSE IF ( SEQL(str, 'node u') ) THEN
              ALLOCATE( element % NodeU(element % NumberOfNodes) )
              READ( str(7:), * ) (element % NodeU(k),k=1,element % NumberOfNodes)

            ELSE IF ( SEQL(str, 'node v') ) THEN
              ALLOCATE( element % NodeV(element % NumberOfNodes) )
              READ( str(7:), * ) (element % NodeV(k),k=1,element % NumberOfNodes)

            ELSE IF ( SEQL(str, 'node w') ) THEN
              ALLOCATE( element % NodeW(element % NumberOfNodes ) )
              READ( str(7:), * ) (element % NodeW(k),k=1,element % NumberOfNodes)

            ELSE IF ( SEQL(str, 'basis') ) THEN
              READ( str(6:), * ) (BasisTerms(k),k=1,element % NumberOfNodes)

            ELSE IF ( SEQL(str, 'stabilization') ) THEN
              READ( str(14:), * ) element % StabilizationMK

            ELSE IF ( SEQL(str, 'gauss points') ) THEN

              Element % GaussPoints2 = 0
              READ( str(13:), *,END=10 ) element % GaussPoints,&
                  element % GaussPoints2, element % GaussPoints0 

10            CONTINUE

              IF ( Element % GaussPoints2 <= 0 ) &
                   Element % GaussPoints2 = Element % GaussPoints

              IF ( Element % GaussPoints0 <= 0 ) &
                   Element % GaussPoints0 = Element % GaussPoints
             
            ELSE IF ( str == 'end element' ) THEN
              gotit = .TRUE.
              EXIT
            END IF
          END DO

          IF ( gotit ) THEN
            Element % StabilizationMK = 0.0d0
            IF ( .NOT.ALLOCATED( element % NodeV ) ) THEN
              ALLOCATE( element % NodeV(element % NumberOfNodes) )
              element % NodeV = 0.0d0
            END IF

            IF ( .NOT.ALLOCATED( element % NodeW ) ) THEN
              ALLOCATE( element % NodeW(element % NumberOfNodes) )
              element % NodeW = 0.0d0
            END IF

            CALL AddElementDescription( element,BasisTerms )
            IF ( ALLOCATED( element % NodeU ) ) DEALLOCATE( element % NodeU )
            IF ( ALLOCATED( element % NodeV ) ) DEALLOCATE( element % NodeV )
            IF ( ALLOCATED( element % NodeW ) ) DEALLOCATE( element % NodeW )
          ELSE
            IF ( ALLOCATED( element % NodeU ) ) DEALLOCATE( element % NodeU )
            IF ( ALLOCATED( element % NodeV ) ) DEALLOCATE( element % NodeV )
            IF ( ALLOCATED( element % NodeW ) ) DEALLOCATE( element % NodeW )
          END IF
        END IF
      END DO

      CLOSE(1)
!------------------------------------------------------------------------------
   END SUBROUTINE InitializeElementDescriptions
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Given element type code return pointer to the corresponding element type
!>    structure.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>   Given element structure return value of a quantity x given at element nodes
!>   at local coordinate point u inside the element. Element basis functions are
!>   used to compute the value. This is for 1D elements, and shouldn't probably
!>   be called directly by the user but through the wrapper routine
!>   InterpolateInElement.
!------------------------------------------------------------------------------
   FUNCTION InterpolateInElement1D( element,x,u ) RESULT(y)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element  !< element structure
     REAL(KIND=dp) :: u          !< Point at which to evaluate the value
     REAL(KIND=dp), DIMENSION(:) :: x  !< Nodal values of the quantity whose value we want to know
     REAL(KIND=dp) :: y                !< value of the quantity y = x(u)
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: s
     INTEGER :: i,j,k,n
     TYPE(ElementType_t), POINTER :: elt
     REAL(KIND=dp), POINTER :: Coeff(:)
     INTEGER, POINTER :: p(:)
     TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
!------------------------------------------------------------------------------

     elt => element % TYPE
     k = Elt % NumberOfNodes
     BasisFunctions => elt % BasisFunctions

     y = 0.0d0
     DO n=1,k
       IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          Coeff => BasisFunctions(n) % Coeff

          s = 0.0d0
          DO i=1,BasisFunctions(n) % n
            IF (p(i)==0) THEN
              s = s + Coeff(i)
            ELSE
              s = s + Coeff(i) * u**p(i)
            END if
          END DO
          y = y + s * x(n)
       END IF
     END DO
   END FUNCTION InterpolateInElement1D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE NodalBasisFunctions1D( y,element,u )
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element  !< element structure
     REAL(KIND=dp) :: u          !< Point at which to evaluate the value
     REAL(KIND=dp) :: y(:)       !< value of the quantity y = x(u)

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: s
     INTEGER :: i,n
     TYPE(ElementType_t), POINTER :: elt
     REAL(KIND=dp), POINTER :: Coeff(:)
     INTEGER, POINTER :: p(:)
     TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
!------------------------------------------------------------------------------

     elt => element % TYPE
     BasisFunctions => elt % BasisFunctions

     DO n=1,Elt % NumberOfNodes
       p => BasisFunctions(n) % p
       Coeff => BasisFunctions(n) % Coeff

       s = 0.0d0
       DO i=1,BasisFunctions(n) % n
         IF (p(i)==0) THEN
           s = s + Coeff(i)
         ELSE
           s = s + Coeff(i) * u**p(i)
         END if
       END DO
       y(n) = s
     END DO
   END SUBROUTINE NodalBasisFunctions1D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of the first partial derivative with
!>   respect to local coordinate of a quantity x given at element nodes at local
!>   coordinate point u inside the element. Element basis functions are used to
!>   compute the value. 
!------------------------------------------------------------------------------
   FUNCTION FirstDerivative1D( element,x,u ) RESULT(y)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element         !< element structure
     REAL(KIND=dp) :: u                 !< Point at which to evaluate the partial derivative
     REAL(KIND=dp), DIMENSION(:) :: x   !< Nodal values of the quantity whose partial derivative we want to know
     REAL(KIND=dp) :: y                 !< value of the quantity y = @x/@u
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,n,l
     TYPE(ElementType_t), POINTER :: elt
     REAL(KIND=dp) :: s
     REAL(KIND=dp), POINTER :: Coeff(:)
     INTEGER, POINTER :: p(:)
     TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)

     elt => element % TYPE
     k = Elt % NumberOfNodes
     BasisFunctions => elt % BasisFunctions

     y = 0.0d0
     DO n=1,k
       IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          Coeff => BasisFunctions(n) % Coeff

          s = 0.0d0
          DO i=1,BasisFunctions(n) % n
             IF ( p(i) >= 1 ) THEN 
                s = s + p(i) * Coeff(i) * u**(p(i)-1)
             END IF
          END DO
          y = y + s * x(n)
       END IF
     END DO
   END FUNCTION FirstDerivative1D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE NodalFirstDerivatives1D( y,element,u )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: u          !< Point at which to evaluate the partial derivative
     REAL(KIND=dp) :: y(:,:)     !< value of the quantity y = @x/@u
     TYPE(Element_t) :: element  !< element structure
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(ElementType_t), POINTER :: elt
     INTEGER :: i,n
     REAL(KIND=dp) :: s

     REAL(KIND=dp), POINTER :: Coeff(:)
     INTEGER, POINTER :: p(:)
     TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)

     elt => element % TYPE
     BasisFunctions => elt % BasisFunctions

     DO n=1, Elt % NumberOfNodes
        p => BasisFunctions(n) % p
        Coeff => BasisFunctions(n) % Coeff

        s = 0.0d0
        DO i=1,BasisFunctions(n) % n
           IF (p(i)>=1) s = s + p(i)*Coeff(i)*u**(p(i)-1)
        END DO
        y(n,1) = s
     END DO
   END SUBROUTINE NodalFirstDerivatives1D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of the second partial derivative with
!>   respect to local coordinate of a quantity x given at element nodes at local
!>   coordinate point u inside the element. Element basis functions are used to
!>   compute the value. 
!------------------------------------------------------------------------------
   FUNCTION SecondDerivatives1D( element,x,u ) RESULT(y)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element          !< element structure
     REAL(KIND=dp) :: u                  !< Point at which to evaluate the partial derivative
     REAL(KIND=dp), DIMENSION(:) :: x    !< Nodal values of the quantity whose partial derivative we want to know
     REAL(KIND=dp) :: y                  !< value of the quantity y = @x/@u
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: usum
     INTEGER :: i,j,k,n
     TYPE(ElementType_t), POINTER :: elt
     INTEGER, POINTER :: p(:),q(:)
     REAL(KIND=dp), POINTER :: Coeff(:)
     REAL(KIND=dp) :: s
     TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)

     elt => element % TYPE
     k = Elt % NumberOfNodes
     BasisFunctions => elt % BasisFunctions

     y = 0.0d0
     DO n=1,k
       IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          Coeff => BasisFunctions(n) % Coeff

          s = 0.0d0
          DO i=1,BasisFunctions(n) % n
             IF ( p(i) >= 2 ) THEN
                s = s + p(i) * (p(i)-1) * Coeff(i) * u**(p(i)-2)
             END IF
          END DO
          y = y + s * x(n)
       END IF
     END DO
   END FUNCTION SecondDerivatives1D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return the value of a quantity x known at element nodes
!>   at local coordinate point (u,v) inside the element. Element basis functions
!>   are used to compute the value. This is for 2D elements, and shouldn't probably
!>   be called directly by the user but through the wrapper routine
!>   InterpolateInElement.
!------------------------------------------------------------------------------
   FUNCTION InterpolateInElement2D( element,x,u,v ) RESULT(y)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element          !< element structure
     REAL(KIND=dp) :: u                  !< u at the point where the quantity is evaluated
     REAL(KIND=dp) :: v                  !< v at the point where the quantity is evaluated
     REAL(KIND=dp), DIMENSION(:) :: x    !< Nodal values of the quantity
     REAL(KIND=dp) :: y                  !< The value of the quantity y = x(u,v)
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: s,t

      INTEGER :: i,j,k,m,n

      TYPE(ElementType_t),POINTER :: elt
      REAL(KIND=dp), POINTER :: Coeff(:)
      INTEGER, POINTER :: p(:),q(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
!------------------------------------------------------------------------------

      elt => element % TYPE
      BasisFunctions => elt % BasisFunctions

      y = 0.0d0
      DO n = 1,elt % NumberOfNodes
        IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          q => BasisFunctions(n) % q
          Coeff => BasisFunctions(n) % Coeff

          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
             s = s + Coeff(i) * u**p(i) * v**q(i)
          END DO
          y = y + s*x(n)
        END IF
      END DO

   END FUNCTION InterpolateInElement2D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE NodalBasisFunctions2D( y,element,u,v )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: y(:)       !< The values of the reference element basis
     TYPE(Element_t) :: element  !< element structure
     REAL(KIND=dp) :: u          !< Point at which to evaluate the value
     REAL(KIND=dp) :: v          !< Point at which to evaluate the value
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: s
     INTEGER :: i,n
     TYPE(ElementType_t), POINTER :: elt
     REAL(KIND=dp), POINTER :: Coeff(:)
     INTEGER, POINTER :: p(:),q(:)
     TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: ult(0:6), vlt(0:6)

     elt => element % TYPE
     BasisFunctions => elt % BasisFunctions

     ult(0) = 1
     ult(1) = u

     vlt(0) = 1
     vlt(1) = v

     DO i=2,elt % BasisFunctionDegree
       ult(i) = u**i
       vlt(i) = v**i
     END DO

     DO n=1,Elt % NumberOfNodes
       p => BasisFunctions(n) % p
       q => BasisFunctions(n) % q
       Coeff => BasisFunctions(n) % Coeff

       s = 0.0d0
       DO i=1,BasisFunctions(n) % n
          s = s + Coeff(i)*ult(p(i))*vlt(q(i))
       END DO
       y(n) = s
     END DO
   END SUBROUTINE NodalBasisFunctions2D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return the value of the first partial derivative with
!>   respect to local coordinate u of a quantity x given at element nodes at local
!>   coordinate point u,v inside the element. Element basis functions are used to
!>   compute the value. 
!------------------------------------------------------------------------------
   FUNCTION FirstDerivativeInU2D( element,x,u,v ) RESULT(y)
!------------------------------------------------------------------------------
      TYPE(Element_t) :: element        !< element structure
      REAL(KIND=dp) :: u,v              !< Point at which to evaluate the partial derivative
      REAL(KIND=dp), DIMENSION(:) :: x  !< Nodal values of the quantity to differentiate
      REAL(KIND=dp) :: y                !< value of the quantity y = @x(u,v)/@u
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: s,t
      TYPE(ElementType_t),POINTER :: elt
      REAL(KIND=dp), POINTER :: Coeff(:)
      INTEGER, POINTER :: p(:),q(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
      INTEGER :: i,j,k,m,n

      elt => element % TYPE
      BasisFunctions => elt % BasisFunctions

      y = 0.0d0
      DO n = 1,elt % NumberOfNodes
        IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          q => BasisFunctions(n) % q
          Coeff => BasisFunctions(n) % Coeff

          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
             IF ( p(i) >= 1 ) THEN
               s = s + p(i) * Coeff(i) * u**(p(i)-1) * v**q(i)
            END IF
          END DO
          y = y + s*x(n)
        END IF
      END DO

   END FUNCTION FirstDerivativeInU2D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of the first partial derivative with
!>   respect to local coordinate v of i quantity x given at element nodes at local
!>   coordinate point u,v inside the element. Element basis functions are used to
!>   compute the value. 
!------------------------------------------------------------------------------
   FUNCTION FirstDerivativeInV2D( element,x,u,v ) RESULT(y)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element        !< element structure
     REAL(KIND=dp) :: u,v              !< Point at which to evaluate the partial derivative
     REAL(KIND=dp), DIMENSION(:) :: x  !< Nodal values of the quantity to differentiate
     REAL(KIND=dp) :: y                !< value of the quantity y = @x(u,v)/@u
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: s,t
      TYPE(ElementType_t),POINTER :: elt
      REAL(KIND=dp), POINTER :: Coeff(:)
      INTEGER, POINTER :: p(:),q(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)

      INTEGER :: i,j,k,m,n

      elt => element % TYPE
      BasisFunctions => elt % BasisFunctions

      y = 0.0d0
      DO n = 1,elt % NumberOfNodes
        IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          q => BasisFunctions(n) % q
          Coeff => BasisFunctions(n) % Coeff

          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
             IF ( q(i) >= 1  ) THEN
                s = s + q(i) * Coeff(i) * u**p(i) * v**(q(i)-1)
             END IF
          END DO
          y = y + s*x(n)
        END IF
      END DO

   END FUNCTION FirstDerivativeInV2D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE NodalFirstDerivatives2D( y,element,u,v )
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element        !< element structure
     REAL(KIND=dp) :: u,v              !< Point at which to evaluate the partial derivative
     REAL(KIND=dp) :: y(:,:)           !< value of the quantity y = @x(u,v)/@u
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: s,t
      TYPE(ElementType_t),POINTER :: elt
      REAL(KIND=dp), POINTER :: Coeff(:)
      INTEGER, POINTER :: p(:),q(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)

      INTEGER :: i,n

      REAL(KIND=dp) :: ult(0:6), vlt(0:6)
 
      elt => element % TYPE
      BasisFunctions => elt % BasisFunctions

      ult(0) = 1
      ult(1) = u

      vlt(0) = 1
      vlt(1) = v

      DO i=2,elt % BasisFunctionDegree
        ult(i) = u**i
        vlt(i) = v**i
      END DO


      DO n = 1,elt % NumberOfNodes
        p => BasisFunctions(n) % p
        q => BasisFunctions(n) % q
        Coeff => BasisFunctions(n) % Coeff

        s = 0.0d0
        t = 0.0d0
        DO i = 1,BasisFunctions(n) % n
          IF (p(i)>=1) s = s + p(i)*Coeff(i)*ult(p(i)-1)*vlt(q(i))
          IF (q(i)>=1) t = t + q(i)*Coeff(i)*ult(p(i))*vlt(q(i)-1)
        END DO
        y(n,1) = s
        y(n,2) = t
      END DO

   END SUBROUTINE NodalFirstDerivatives2D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given an element structure return the second partial derivatives of 
!>   a quantity x given at the element nodes with respect to the local coordinates
!>   u,v of the element. The element basis functions are used to compute the value. 
!------------------------------------------------------------------------------
   FUNCTION SecondDerivatives2D( element,x,u,v ) RESULT(ddx)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element        !< Element structure
     REAL(KIND=dp) :: u,v              !< Point at which to evaluate the partial derivatives
     REAL(KIND=dp), DIMENSION(:) :: x  !< The nodal values of the quantity to differentiate
     REAL(KIND=dp), DIMENSION (2,2) :: ddx !< The second partial derivatives of x
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      TYPE(ElementType_t),POINTER :: elt
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
      REAL(KIND=dp) :: s,t
      INTEGER, POINTER :: p(:),q(:)
      REAL(KIND=dp), POINTER :: Coeff(:)
      INTEGER :: i,j,k,n,m

!------------------------------------------------------------------------------
      elt => element % TYPE
      k = elt % NumberOfNodes
      BasisFunctions => elt % BasisFunctions

      ddx = 0.0d0
      DO n = 1,k
        IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          q => BasisFunctions(n) % q
          Coeff => BasisFunctions(n) % Coeff
!------------------------------------------------------------------------------
!         @^2x/@u^2
!------------------------------------------------------------------------------
          s = 0.0d0
          DO i = 1, BasisFunctions(n) % n
             IF ( p(i) >= 2 ) THEN
                s = s + p(i) * (p(i)-1) * Coeff(i) * u**(p(i)-2) * v**q(i)
             END IF
          END DO
          ddx(1,1) = ddx(1,1) + s*x(n)

!------------------------------------------------------------------------------
!         @^2x/@u@v
!------------------------------------------------------------------------------
          s = 0.0d0
          DO i = 1, BasisFunctions(n) % n
              IF ( p(i) >= 1 .AND. q(i) >= 1 ) THEN
                 s = s + p(i) * q(i) * Coeff(i) * u**(p(i)-1) * v**(q(i)-1)
              END IF
          END DO
          ddx(1,2) = ddx(1,2) + s*x(n)

!------------------------------------------------------------------------------
!         @^2x/@v^2
!------------------------------------------------------------------------------
          s = 0.0d0
          DO i = 1, BasisFunctions(n) % n
             IF ( q(i) >= 2 ) THEN
                s = s + q(i) * (q(i)-1) * Coeff(i) * u**p(i) * v**(q(i)-2)
             END IF
          END DO
          ddx(2,2) = ddx(2,2) + s*x(n)
        END IF
      END DO

      ddx(2,1) = ddx(1,2)

   END FUNCTION SecondDerivatives2D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of a quantity x given at element nodes
!>   at local coordinate point (u,v,w) inside the element. Element basis functions
!>   are used to compute the value. This is for 3D elements, and shouldn't probably
!>   be called directly by the user but through the wrapper routine
!>   InterpolateInElement.
!------------------------------------------------------------------------------
   FUNCTION InterpolateInElement3D( element,x,u,v,w ) RESULT(y)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element        !< element structure
     REAL(KIND=dp) :: u,v,w            !< Point at which to evaluate the partial derivative
     REAL(KIND=dp), DIMENSION(:) :: x  !< Nodal values of the quantity to differentiate
     REAL(KIND=dp) :: y                !< value of the quantity y = x(u,v,w)
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      TYPE(ElementType_t),POINTER :: elt
      INTEGER :: i,j,k,l,n,m
      REAL(KIND=dp) :: s,t
      INTEGER, POINTER :: p(:),q(:), r(:)
      REAL(KIND=dp), POINTER :: Coeff(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
!------------------------------------------------------------------------------

      elt => element % TYPE
      l = elt % BasisFunctionDegree
      BasisFunctions => elt % BasisFunctions

      IF ( Elt % ElementCode == 605 ) THEN
        s = 0.0d0
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        y = 0.0d0
        DO n=1,5
          IF(x(n)==0) CYCLE
          SELECT CASE(n)
          CASE(1)
            y = y + x(1)*((1-u)*(1-v) - w + u*v*w * s) / 4
          CASE(2)
            y = y + x(2)*((1+u)*(1-v) - w - u*v*w * s) / 4
          CASE(3)
            y = y + x(3)*((1+u)*(1+v) - w + u*v*w * s) / 4
          CASE(4)
            y = y + x(4)*((1-u)*(1+v) - w - u*v*w * s) / 4
          CASE(5)
            y = y + x(5)*w
          END SELECT
        END DO
        RETURN
      ELSE IF ( Elt % ElementCode == 613 ) THEN
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        y = 0.0d0
        DO n=1,13
          IF(x(n)==0) CYCLE
          SELECT CASE(n)
          CASE(1)
            y = y + x(1)  * (-u-v-1) * ( (1-u) * (1-v) - w + u*v*w * s ) / 4
          CASE(2)
            y = y + x(2)  * ( u-v-1) * ( (1+u) * (1-v) - w - u*v*w * s ) / 4
          CASE(3)
            y = y + x(3)  * ( u+v-1) * ( (1+u) * (1+v) - w + u*v*w * s ) / 4
          CASE(4)
            y = y + x(4)  * (-u+v-1) * ( (1-u) * (1+v) - w - u*v*w * s ) / 4
          CASE(5)
            y = y + x(5)  * w*(2*w-1)
          CASE(6)
            y = y + x(6)  * (1+u-w)*(1-u-w)*(1-v-w) * s / 2
          CASE(7)
            y = y + x(7)  * (1+v-w)*(1-v-w)*(1+u-w) * s / 2
          CASE(8)
            y = y + x(8)  * (1+u-w)*(1-u-w)*(1+v-w) * s / 2
          CASE(9)
            y = y + x(9)  * (1+v-w)*(1-v-w)*(1-u-w) * s / 2
          CASE(10)
            y = y + x(10) * w * (1-u-w) * (1-v-w) * s
          CASE(11)
            y = y + x(11) * w * (1+u-w) * (1-v-w) * s
          CASE(12)
            y = y + x(12) * w * (1+u-w) * (1+v-w) * s
          CASE(13)
            y = y + x(13) * w * (1-u-w) * (1+v-w) * s
          END SELECT
        END DO
        RETURN
      END IF

      y = 0.0d0
      DO n = 1,elt % NumberOfNodes
        IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          q => BasisFunctions(n) % q
          r => BasisFunctions(n) % r
          Coeff => BasisFunctions(n) % Coeff

          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
             s = s + Coeff(i) * u**p(i) * v**q(i) * w**r(i)
          END DO
          y = y + s*x(n)
        END IF
      END DO
!------------------------------------------------------------------------------
   END FUNCTION InterpolateInElement3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE NodalBasisFunctions3D( y,element,u,v,w )
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element        !< element structure
     REAL(KIND=dp) :: u,v,w            !< Point at which to evaluate the basis functions
     REAL(KIND=dp) :: y(:)             !< The values of the basis functions
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: s

     INTEGER :: i,n

     TYPE(ElementType_t), POINTER :: elt

     REAL(KIND=dp), POINTER :: Coeff(:)
     INTEGER, POINTER :: p(:),q(:),r(:)
     TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: ult(0:6), vlt(0:6), wlt(0:6)

     elt => element % TYPE
     BasisFunctions => elt % BasisFunctions
 
     ult(0) = 1
     ult(1) = u

     vlt(0) = 1
     vlt(1) = v

     wlt(0) = 1
     wlt(1) = w

     DO i=2,elt % BasisFunctionDegree
       ult(i) = u**i
       vlt(i) = v**i
       wlt(i) = w**i
     END DO

     DO n=1,Elt % NumberOfNodes
       p => BasisFunctions(n) % p
       q => BasisFunctions(n) % q
       r => BasisFunctions(n) % r
       Coeff => BasisFunctions(n) % Coeff

       s = 0.0d0
       DO i=1,BasisFunctions(n) % n
          s = s + Coeff(i)*ult(p(i))*vlt(q(i))*wlt(r(i))
       END DO
       y(n) = s
     END DO
   END SUBROUTINE NodalBasisFunctions3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>   Given element structure return value of the first partial derivative with
!>   respect to local coordinate u of a quantity x given at element nodes at
!>   local coordinate point u,v,w inside the element. Element basis functions
!>   are used to compute the value. 
!------------------------------------------------------------------------------
   FUNCTION FirstDerivativeInU3D( element,x,u,v,w ) RESULT(y)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element        !< element structure
     REAL(KIND=dp) :: u,v,w            !< Point at which to evaluate the partial derivative
     REAL(KIND=dp), DIMENSION(:) :: x  !< Nodal values of the quantity to be derivated
     REAL(KIND=dp) :: y                !< value of the quantity y =  @x(u,v,w)/@u
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      TYPE(ElementType_t),POINTER :: elt
      INTEGER :: i,j,k,l,n,m
      REAL(KIND=dp) :: s,t
      INTEGER, POINTER :: p(:),q(:), r(:)
      REAL(KIND=dp), POINTER :: Coeff(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
!------------------------------------------------------------------------------
      elt => element % TYPE
      l = elt % BasisFunctionDegree
      BasisFunctions => elt % BasisFunctions

      IF ( Elt % ElementCode == 605 ) THEN
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        y = 0.0d0
        DO n=1,5
          IF(x(n)==0) CYCLE
          SELECT CASE(n)
          CASE(1)
            y = y + x(1) * ( -(1-v) + v*w * s ) / 4
          CASE(2)
            y = y + x(2) * (  (1-v) - v*w * s ) / 4
          CASE(3)
            y = y + x(3) * (  (1+v) + v*w * s ) / 4
          CASE(4)
            y = y + x(4) * ( -(1+v) - v*w * s ) / 4
          CASE(5)
            CONTINUE
          END SELECT
        END DO
        RETURN

      ELSE IF ( Elt % ElementCode == 613 ) THEN
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        y = 0.0d0
        DO n=1,13
          IF(x(n)==0) CYCLE
          SELECT CASE(n)
          CASE(1)
             y = y + x(1) * (-((1-u)*(1-v)-w+u*v*w*s)+(-u-v-1) * (-(1-v)+v*w*s))/4
          CASE(2)
             y = y + x(2)  * ( ((1+u)*(1-v)-w-u*v*w*s)+( u-v-1) * ( (1-v)-v*w*s))/4
          CASE(3)
             y = y + x(3)  * ( ((1+u)*(1+v)-w+u*v*w*s)+( u+v-1) * ( (1+v)+v*w*s))/4
          CASE(4)
             y = y + x(4)  * (-((1-u)*(1+v)-w-u*v*w*s)+(-u+v-1) * (-(1+v)-v*w*s))/4
          CASE(5)
             CONTINUE
          CASE(6)
             y = y + x(6)  * (  (1-u-w)*(1-v-w) - (1+u-w)*(1-v-w) ) * s / 2
          CASE(7)
             y = y + x(7)  * (  (1+v-w)*(1-v-w) ) * s / 2
          CASE(8)
             y = y + x(8)  * (  (1-u-w)*(1+v-w) - (1+u-w)*(1+v-w) ) * s / 2
          CASE(9)
             y = y + x(9)  * ( -(1+v-w)*(1-v-w) ) * s / 2
          CASE(10)
             y = y - x(10) * w * (1-v-w) * s
          CASE(11)
             y = y + x(11) * w * (1-v-w) * s
          CASE(12)
             y = y + x(12) * w * (1+v-w) * s
          CASE(13)
             y = y - x(13) * w * (1+v-w) * s
          END SELECT
        END DO
        RETURN
      END IF

      y = 0.0d0
      DO n = 1,elt % NumberOfNodes
        IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          q => BasisFunctions(n) % q
          r => BasisFunctions(n) % r
          Coeff => BasisFunctions(n) % Coeff

          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
             IF ( p(i) >= 1  ) THEN
                s = s + p(i) * Coeff(i) * u**(p(i)-1) * v**q(i) * w**r(i)
             END IF
          END DO
          y = y + s*x(n)
        END IF
      END DO
!------------------------------------------------------------------------------
   END FUNCTION FirstDerivativeInU3D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of the first partial derivative with
!>   respect to local coordinate v of a quantity x given at element nodes at
!>   local coordinate point u,v,w inside the element. Element basis functions
!>   are used to compute the value. 
!------------------------------------------------------------------------------
   FUNCTION FirstDerivativeInV3D( element,x,u,v,w ) RESULT(y)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element        !< element structure
     REAL(KIND=dp) :: u,v,w            !< Point at which to evaluate the partial derivative
     REAL(KIND=dp), DIMENSION(:) :: x  !< Nodal values of the quantity to be derivated
     REAL(KIND=dp) :: y                !< value of the quantity y =  @x(u,v,w)/@v
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      TYPE(ElementType_t),POINTER :: elt
      INTEGER :: i,j,k,l,n,m
      REAL(KIND=dp) :: s,t
      INTEGER, POINTER :: p(:),q(:), r(:)
      REAL(KIND=dp), POINTER :: Coeff(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
!------------------------------------------------------------------------------
      elt => element % TYPE
      l = elt % BasisFunctionDegree
      BasisFunctions => elt % BasisFunctions

      IF ( Elt % ElementCode == 605 ) THEN
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        y = 0.0d0
        DO n=1,5
          IF(x(n)==0) CYCLE
          SELECT CASE(n)
          CASE(1)
            y = y + x(1) * ( -(1-u) + u*w * s ) / 4
          CASE(2)
            y = y + x(2) * ( -(1+u) - u*w * s ) / 4
          CASE(3)
            y = y + x(3) * (  (1+u) + u*w * s ) / 4
          CASE(4)
            y = y + x(4) * (  (1-u) - u*w * s ) / 4
          CASE(5)
            CONTINUE
          END SELECT
        END DO
        RETURN
      ELSE IF ( Elt % ElementCode == 613 ) THEN
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        y = 0.0d0
        DO n=1,13
          IF(x(n)==0) CYCLE
          SELECT CASE(n)
          CASE(1)
            y = y + x(1)  * ( -( (1-u) * (1-v) - w + u*v*w * s ) +  &
                (-u-v-1) * ( -(1-u) + u*w * s ) ) / 4
          CASE(2)
            y = y + x(2)  * ( -( (1+u) * (1-v) - w - u*v*w * s ) + &
                ( u-v-1) * ( -(1+u) - u*w * s ) ) / 4
          CASE(3)
            y = y + x(3)  * (  ( (1+u) * (1+v) - w + u*v*w * s ) + &
                ( u+v-1) * (  (1+u) + u*w * s ) ) / 4
          CASE(4)
            y = y + x(4)  * (  ( (1-u) * (1+v) - w - u*v*w * s ) + &
                (-u+v-1) * (  (1-u) - u*w * s ) ) / 4
          CASE(5)
            CONTINUE
          CASE(6)
            y = y - x(6)  *  (1+u-w)*(1-u-w) * s / 2
          CASE(7)
            y = y + x(7)  * ( (1-v-w)*(1+u-w) - (1+v-w)*(1+u-w) ) * s / 2
          CASE(8)
            y = y + x(8)  *  (1+u-w)*(1-u-w) * s / 2
          CASE(9)
            y = y + x(9)  * ( (1-v-w)*(1-u-w) - (1+v-w)*(1-u-w) ) * s / 2
          CASE(10)
            y = y - x(10) *  w * (1-u-w) * s
          CASE(11)
            y = y - x(11) *  w * (1+u-w) * s
          CASE(12)
            y = y + x(12) *  w * (1+u-w) * s
          CASE(13)
            y = y + x(13) *  w * (1-u-w) * s
          END SELECT
        END DO
        RETURN
      END IF

      y = 0.0d0
      DO n = 1,elt % NumberOfNodes
        IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          q => BasisFunctions(n) % q
          r => BasisFunctions(n) % r
          Coeff => BasisFunctions(n) % Coeff

          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
             IF ( q(i) >= 1  ) THEN
                s = s + q(i) * Coeff(i) * u**p(i) * v**(q(i)-1) * w**r(i)
             END IF
          END DO
          y = y + s*x(n)
        END IF
      END DO
   END FUNCTION FirstDerivativeInV3D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of the first partial derivatives with
!>   respect to local coordinate w of a quantity x given at element nodes at
!>   local coordinate point u,v,w inside the element. Element basis functions
!>   are used to compute the value. 
!------------------------------------------------------------------------------
   FUNCTION FirstDerivativeInW3D( element,x,u,v,w ) RESULT(y)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element        !< element structure
     REAL(KIND=dp) :: u,v,w            !< Point at which to evaluate the partial derivative
     REAL(KIND=dp), DIMENSION(:) :: x  !< Nodal values of the quantity to be derivated
     REAL(KIND=dp) :: y                !< value of the quantity y =  @x(u,v,w)/@w
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      TYPE(ElementType_t),POINTER :: elt
      INTEGER :: i,j,k,l,n,m
      REAL(KIND=dp) :: s,t
      INTEGER, POINTER :: p(:),q(:), r(:)
      REAL(KIND=dp), POINTER :: Coeff(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
!------------------------------------------------------------------------------
      elt => element % TYPE
      l = elt % BasisFunctionDegree
      BasisFunctions => elt % BasisFunctions

      IF ( Elt % ElementCode == 605 ) THEN
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        y = 0.0d0
        DO n=1,5
          IF(x(n)==0) CYCLE
          SELECT CASE(n)
          CASE(1)
            y = y + x(1) * ( -1 + u*v*s**2 ) / 4
          CASE(2)
            y = y + x(2) * ( -1 - u*v*s**2 ) / 4
          CASE(3)
            y = y + x(3) * ( -1 + u*v*s**2 ) / 4
          CASE(4)
            y = y + x(4) * ( -1 - u*v*s**2 ) / 4
          CASE(5)
            y = y + x(5)
          END SELECT
        END DO
        RETURN
      ELSE IF ( Elt % ElementCode == 613 ) THEN
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        y = 0.0d0
        DO n=1,13
          IF(x(n)==0) CYCLE
          SELECT CASE(n)
          CASE(1)
            y = y + x(1)  * (-u-v-1) * ( -1 + u*v*s**2 ) / 4
          CASE(2)
            y = y + x(2)  * ( u-v-1) * ( -1 - u*v*s**2 ) / 4
          CASE(3)
            y = y + x(3)  * ( u+v-1) * ( -1 + u*v*s**2 ) / 4
          CASE(4)
            y = y + x(4)  * (-u+v-1) * ( -1 - u*v*s**2 ) / 4
          CASE(5)
            y = y + x(5)  * (4*w-1)
          CASE(6)
            y = y + x(6)  * ( ( -(1-u-w)*(1-v-w) - (1+u-w)*(1-v-w) - (1+u-w)*(1-u-w) ) * s + &
                ( 1+u-w)*(1-u-w)*(1-v-w) * s**2 ) / 2
          CASE(7)
            y = y + x(7)  * ( ( -(1-v-w)*(1+u-w) - (1+v-w)*(1+u-w) - (1+v-w)*(1-v-w) ) * s + &
                ( 1+v-w)*(1-v-w)*(1+u-w) * s**2 ) / 2
          CASE(8)
            y = y + x(8)  * ( ( -(1-u-w)*(1+v-w) - (1+u-w)*(1+v-w) - (1+u-w)*(1-u-w) ) * s + &
                ( 1+u-w)*(1-u-w)*(1+v-w) * s**2 ) / 2
          CASE(9)
            y = y + x(9)  * ( ( -(1-v-w)*(1-u-w) - (1+v-w)*(1-u-w) - (1+v-w)*(1-v-w) ) * s + &
                ( 1+v-w)*(1-v-w)*(1-u-w) * s**2 ) / 2
          CASE(10)
            y = y + x(10) * ( ( (1-u-w) * (1-v-w) - w * (1-v-w) - w * (1-u-w) ) * s  + &
                w * (1-u-w) * (1-v-w) * s**2 )
          CASE(11)
            y = y + x(11) * ( ( (1+u-w) * (1-v-w) - w * (1-v-w) - w * (1+u-w) ) * s  + &
                w * (1+u-w) * (1-v-w) * s**2 )
          CASE(12)
            y = y + x(12) * ( ( (1+u-w) * (1+v-w) - w * (1+v-w) - w * (1+u-w) ) * s  + &
                w * (1+u-w) * (1+v-w) * s**2 )
          CASE(13)
            y = y + x(13) * ( ( (1-u-w) * (1+v-w) - w * (1+v-w) - w * (1-u-w) ) * s  + &
                w * (1-u-w) * (1+v-w) * s**2 )
          END SELECT
        END DO
        RETURN
      END IF

      y = 0.0d0
      DO n = 1,elt % NumberOfNodes
        IF ( x(n) /= 0.0d0 ) THEN
          p => BasisFunctions(n) % p
          q => BasisFunctions(n) % q
          r => BasisFunctions(n) % r
          Coeff => BasisFunctions(n) % Coeff

          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
             IF ( r(i) >= 1  ) THEN
                s = s + r(i) * Coeff(i) * u**p(i) * v**q(i) * w**(r(i)-1)
             END IF
          END DO
          y = y + s*x(n)
        END IF
      END DO
!------------------------------------------------------------------------------
   END FUNCTION FirstDerivativeInW3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Return first partial derivative in u of a quantity x at point (u,v,w)
!------------------------------------------------------------------------------
   SUBROUTINE NodalFirstDerivatives3D( y,element,u,v,w )
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element        !< element structure
     REAL(KIND=dp) :: u,v,w            !< Point at which to evaluate the partial derivative
     REAL(KIND=dp) :: y(:,:)           !< value of the quantity y =  @x(u,v,w)/@u
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: s,t,z
      TYPE(ElementType_t),POINTER :: elt
      REAL(KIND=dp), POINTER :: Coeff(:)
      INTEGER, POINTER :: p(:),q(:),r(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
      INTEGER :: i,n
      REAL(KIND=dp) :: ult(0:6), vlt(0:6), wlt(0:6)
 
      elt => element % TYPE
      BasisFunctions => elt % BasisFunctions
 
      ult(0) = 1
      ult(1) = u

      vlt(0) = 1
      vlt(1) = v
 
      wlt(0) = 1
      wlt(1) = w

      DO i=2,elt % BasisFunctionDegree
        ult(i) = u**i
        vlt(i) = v**i
        wlt(i) = w**i
      END DO

      DO n = 1,elt % NumberOfNodes
        p => BasisFunctions(n) % p
        q => BasisFunctions(n) % q
        r => BasisFunctions(n) % r
        Coeff => BasisFunctions(n) % Coeff

        s = 0.0d0
        t = 0.0d0
        z = 0.0d0
        DO i = 1,BasisFunctions(n) % n
          IF (p(i)>=1) s = s + p(i)*Coeff(i)*ult(p(i)-1)*vlt(q(i))*wlt(r(i))
          IF (q(i)>=1) t = t + q(i)*Coeff(i)*ult(p(i))*vlt(q(i)-1)*wlt(r(i))
          IF (r(i)>=1) z = z + r(i)*Coeff(i)*ult(p(i))*vlt(q(i))*wlt(r(i)-1)
        END DO
        y(n,1) = s
        y(n,2) = t
        y(n,3) = z
      END DO
   END SUBROUTINE NodalFirstDerivatives3D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given the element structure return the second partial derivatives of 
!>   a quantity x given at element nodes with respect to local coordinates
!>   at a point with local coordinates (u,v,w) inside the element. Element basis
!>   functions are used to compute the value. 
!------------------------------------------------------------------------------
   FUNCTION SecondDerivatives3D( element,x,u,v,w ) RESULT(ddx)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: x(:)
!     INPUT: Nodal values of the quantity whose partial derivatives we want to know
!
!    REAL(KIND=dp) :: u,v,w
!     INPUT: Point at which to evaluate the partial derivative
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: s
!      value of the quantity s = @^2x(u,v)/@v^2
!    
!------------------------------------------------------------------------------
   !
   !  Return matrix of second partial derivatives.
   !
!------------------------------------------------------------------------------

      TYPE(Element_t) :: element

      REAL(KIND=dp), DIMENSION(:) :: x
      REAL(KIND=dp) :: u,v,w

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      TYPE(ElementType_t),POINTER :: elt
      REAL(KIND=dp), DIMENSION (3,3) :: ddx
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)

      REAL(KIND=dp), POINTER :: Coeff(:)
      INTEGER, POINTER :: p(:), q(:), r(:)

      REAL(KIND=dp) :: s,t
      INTEGER :: i,j,k,l,n,m

!------------------------------------------------------------------------------
      elt => element % TYPE
      k = elt % NumberOfNodes
      BasisFunctions => elt % BasisFunctions

      ddx = 0.0d0
      IF ( Elt % ElementCode == 605 ) THEN
        s = 0.0d0
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        ddx(1,2) = (x(1)-x(2)+x(3)-x(4))*(1+w*s)
        ddx(2,1) = ddx(1,2)

        ddx(1,3) = (x(1)-x(2)+x(3)-x(4))*v*s**2
        ddx(3,1) = ddx(1,3)

        ddx(2,3) = (x(1)-x(2)+x(3)-x(4))*u*s**2
        ddx(3,2) = ddx(2,3)
        ddx = ddx/4

        RETURN
      ELSE IF ( Elt % ElementCode == 613 ) THEN
        s = 0.0d0
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        DO n=1,13
          IF(x(n)==0) CYCLE

          t = 0
          SELECT CASE(n)
          CASE(1)
            t = t - x(1)  * (-(1-v) + v*w*s)/2
          CASE(2)
            t = t + x(2)  * ( (1-v) - v*w*s)/2
          CASE(3)
            t = t + x(3)  * ( (1+v) + v*w*s)/2
          CASE(4)
            t = t - x(4)  * (-(1+v) - v*w*s)/2
          CASE(6)
            t = t - x(6)  *  (1-v-w) * s
          CASE(8)
            t = t - x(8)  *  (1+v-w) * s
          END SELECT
          ddx(1,1) = ddx(1,1) + t

          t = 0
          SELECT CASE(n)
          CASE(1)
            t = t + x(1)  *  ((1-u) - u*w*s)/4
            t = t + x(1)  *  ((1-v) - v*w*s)/4
            t = t + x(1)  *  (-u-v-1)*(1+w*s)/4
          CASE(2)
            t = t + x(2)  *  (-(1+u) - u*w*s)/4
            t = t + x(2)  *  ( -(1-v) + v*w*s)/4
            t = t + x(2)  *  ( u-v-1)*(-1-w*s)/4
          CASE(3)
            t = t + x(3)  *  ( (1+u) + u*w*s)/4
            t = t + x(3)  *  ( (1+v) + v*w*s)/4
            t = t + x(3)  *  ( u+v-1)*(1+w*s)/4
          CASE(4)
            t = t + x(4)  *  ( -(1-u) + u*w*s)/4
            t = t + x(4)  *  (-(1+v) - v*w*s)/4
            t = t + x(4)  *  (-u+v-1)*(-1-w*s)/4
          CASE(5)
            CONTINUE
          CASE(6)
            t = t - x(6)  * (1-u-w)*s/2
            t = t + x(6)  * (1+u-w)*s/2
          CASE(7)
            t = t + x(7)  * (1-v-w)*s/2
            t = t - x(7)  * (1+v-w)*s/2
          CASE(8)
            t = t + x(8)  * (1-u-w)*s/2
            t = t - x(8)  * (1+u-w)*s/2
          CASE(9)
            t = t - x(9)  * (1-v-w)*s/2
            t = t + x(9)  * (1+v-w)*s/2
          CASE(10)
            t = t + x(10) *  w*s
          CASE(11)
            t = t - x(11) *  w*s
          CASE(12)
            t = t + x(12) *  w*s
          CASE(13)
            t = t - x(13) *  w*s
          END SELECT
          ddx(1,2) = ddx(1,2) + t

          t = 0
          SELECT CASE(n) 
          CASE(1)
            t = t - x(1)  * (-1 + u*v*s**2) / 4
            t = t + x(1)  * (-u-v-1) * (v*s**2) / 4
          CASE(2)
            t = t + x(2)  * (-1 - u*v*s**2) / 4
            t = t + x(2)  * ( u-v-1) * (-v*s**2) / 4
          CASE(3)
            t = t + x(3)  * (-1 + u*v*s**2) / 4
            t = t + x(3)  * ( u+v-1) * (v*s**2) / 4
          CASE(4)
            t = t - x(4)  * (-1 - u*v*s**2) / 4
            t = t + x(4)  * (-u+v-1) * (-v*s**2) / 4
          CASE(5)
            CONTINUE
          CASE(6)
            t = t - x(6)  * (1-v-w) * s / 2
            t = t - x(6)  * (1-u-w) * s / 2
            t = t + x(6)  * (1-u-w)*(1-v-w) * s**2 / 2
            t = t + x(6)  * (1-v-w) * s / 2
            t = t + x(6)  * (1+u-w) * s / 2
            t = t - x(6)  * (1+u-w)*(1-v-w) * s**2 / 2
          CASE(7)
            t = t - x(7)  * (1-v-w) * s / 2
            t = t - x(7)  * (1+v-w) * s / 2
            t = t + x(7)  * (1+v-w)*(1-v-w) * s**2 / 2
          CASE(8)
            t = t - x(8)  * (1+v-w) * s / 2
            t = t - x(8)  * (1-u-w) * s / 2
            t = t + x(8)  * (1-u-w)*(1+v-w) * s**2 / 2
            t = t + x(8)  * (1+v-w) * s / 2
            t = t + x(8)  * (1+u-w) * s / 2
            t = t - x(8)  * (1+u-w)*(1+v-w) * s**2 / 2
          CASE(9)
            t = t + x(9)  * (1-v-w) * s / 2
            t = t + x(9)  * (1+v-w) * s / 2
            t = t - x(9)  * (1+v-w)*(1-v-w) * s**2 / 2
          CASE(10)
            t = t + x(10) * w * s
            t = t - x(10) * (1-v-w) * s**2
          CASE(11)
            t = t - x(11) * w * s
            t = t + x(11) * (1-v-w) * s**2
          CASE(12)
            t = t - x(12) * w * s
            t = t + x(12) * (1+v-w) * s**2
          CASE(13)
            t = t + x(13) * w * s
            t = t - x(13) * (1+v-w) * s**2
          END SELECT
          ddx(1,3) = ddx(1,3) + t

          t = 0
          SELECT CASE(n)
          CASE(1)
            t = t - x(1)  * (-(1-u) + u*w*s)/2
          CASE(2)
            t = t - x(2)  * (-(1+u) - u*w*s)/2
          CASE(3)
            t = t + x(3)  * ( (1+u) + u*w*s)/2
          CASE(4)
            t = t + x(4)  * ( (1-u) - u*w*s)/2
          CASE(7)
            t = t - x(7)  * (1+u-w)*s
          CASE(9)
            t = t - x(9)  * (1-u-w)*s
          CASE(6,8,10,11,12,13)
          END SELECT
          ddx(2,2) = ddx(2,2) + t

          t = 0
          SELECT CASE(n)
          CASE(1)
            t = t - x(1)  * (-1 + u*v*s**2) / 4
            t = t + x(1)  * (-u-v-1) * (u*s**2) / 4
          CASE(2)
            t = t - x(2)  * (-1 - u*v*s**2) / 4
            t = t + x(2)  * ( u-v-1) * (-u*s**2) / 4
          CASE(3)
            t = t + x(3)  * (-1 + u*v*s**2) / 4
            t = t + x(3)  * ( u+v-1) * (u*s**2) / 4
          CASE(4)
            t = t + x(4)  * (-1 - u*v*s**2) / 4
            t = t + x(4)  * (-u+v-1) * (-u*s**2) / 4
          CASE(5)
            CONTINUE
          CASE(6)
            t = t + x(6)  * (1-u-w) * s / 2
            t = t + x(6)  * (1+u-w) * s / 2
            t = t - x(6)  * (1+u-w)*(1-u-w) * s**2 / 2
          CASE(7)
            t = t - x(7)  * (1+u-w) * s / 2
            t = t - x(7)  * (1-v-w) * s / 2
            t = t + x(7)  * (1-v-w)*(1+u-w) * s**2 / 2
            t = t + x(7)  * (1+u-w) * s / 2
            t = t + x(7)  * (1+v-w) * s / 2
            t = t - x(7)  * (1+v-w)*(1+u-w) * s**2 / 2
          CASE(8)
            t = t - x(8)  * (1-u-w) * s / 2
            t = t - x(8)  * (1+u-w) * s / 2
            t = t + x(8)  * (1+u-w)*(1-u-w) * s**2 / 2
          CASE(9)
            t = t - x(9)  * (1-u-w) * s / 2
            t = t - x(9)  * (1-v-w) * s / 2
            t = t + x(9)  * (1-v-w)*(1-u-w) * s**2 / 2
            t = t + x(9)  * (1-u-w) * s / 2
            t = t + x(9)  * (1+v-w) * s / 2
            t = t - x(9)  * (1+v-w)*(1-u-w) * s**2 / 2
          CASE(10)
            t = t + x(10) * w * s
            t = t - x(10) * (1-u-w) * s**2
          CASE(11)
            t = t + x(11) * w * s
            t = t - x(11) * (1+u-w) * s**2
          CASE(12)
            t = t - x(12) * w * s
            t = t + x(12) * (1+u-w) * s**2
          CASE(13)
            t = t - x(13) * w * s
            t = t + x(13) * (1-u-w) * s**2
          END SELECT
          ddx(2,3) = ddx(2,3) + t

          t = 0
          SELECT CASE(n)
          CASE(1)
            t = t + x(1)  * (-u-v-1) * ( u*v*2*s**3) / 4
          CASE(2)
            t = t + x(2)  * ( u-v-1) * (-u*v*2*s**3) / 4
          CASE(3)
            t = t + x(3)  * ( u+v-1) * ( u*v*2*s**3) / 4
          CASE(4)
            t = t + x(4)  * (-u+v-1) * (-u*v*2*s**3) / 4
          CASE(5)
            t = t + x(5) * 4
          CASE(6)
            t = t + x(6)  * (1-v-w) * s / 2
            t = t + x(6)  * (1-u-w) * s / 2
            t = t - x(6)  * (1-u-w)*(1-v-w) * s**2 / 2
            t = t + x(6)  * (1-v-w) * s / 2
            t = t + x(6)  * (1+u-w) * s / 2
            t = t - x(6)  * (1+u-w)*(1-v-w) * s**2 / 2
            t = t + x(6)  * (1-u-w) * s / 2
            t = t + x(6)  * (1+u-w) * s / 2
            t = t - x(6)  * (1+u-w)*(1-u-w) * s**2 / 2
            t = t - x(6)  * (1-u-w)*(1-v-w) * s**2 / 2
            t = t - x(6)  * (1+u-w)*(1-v-w) * s**2 / 2
            t = t - x(6)  * (1+u-w)*(1-u-w) * s**2 / 2
            t = t + x(6)  * (1+u-w)*(1-u-w)*(1-v-w) * 2*s**3 / 2
          CASE(7)
            t = t + x(7)  * (1+u-w) * s / 2
            t = t + x(7)  * (1-v-w) * s / 2
            t = t - x(7)  * (1-v-w)*(1+u-w) * s**2 / 2
            t = t + x(7)  * (1+u-w) * s / 2
            t = t + x(7)  * (1+v-w) * s / 2
            t = t - x(7)  * (1+v-w)*(1+u-w) * s**2 / 2
            t = t + x(7)  * (1-v-w) * s / 2
            t = t + x(7)  * (1+v-w) * s / 2
            t = t - x(7)  * (1+v-w)*(1-v-w) * s**2 / 2
            t = t - x(7)  * (1-v-w)*(1+u-w) * s**2 / 2
            t = t - x(7)  * (1+v-w)*(1+u-w) * s**2 / 2
            t = t - x(7)  * (1+v-w)*(1-v-w) * s**2 / 2
            t = t + x(7)  * (1+v-w)*(1-v-w)*(1+u-w) * 2*s**3 / 2
          CASE(8)
            t = t + x(8)  * (1+v-w) * s / 2
            t = t + x(8)  * (1-u-w) * s / 2
            t = t - x(8)  * (1-u-w)*(1+v-w) * s**2 / 2
            t = t + x(8)  * (1+v-w) * s / 2
            t = t + x(8)  * (1+u-w) * s / 2
            t = t - x(8)  * (1+u-w)*(1+v-w) * s**2 / 2
            t = t + x(8)  * (1-u-w) * s / 2
            t = t + x(8)  * (1+u-w) * s / 2
            t = t - x(8)  * (1+u-w)*(1-u-w) * s**2 / 2
            t = t - x(8)  * (1-u-w)*(1+v-w) * s**2 / 2
            t = t - x(8)  * (1+u-w)*(1+v-w) * s**2 / 2
            t = t - x(8)  * (1+u-w)*(1-u-w) * s**2 / 2
            t = t + x(8)  * (1+u-w)*(1-u-w)*(1+v-w) * 2*s**3 / 2
          CASE(9)
            t = t + x(9)  * (1-u-w) * s / 2
            t = t + x(9)  * (1-v-w) * s / 2
            t = t - x(9)  * (1-v-w)*(1-u-w) * s**2 / 2
            t = t + x(9)  * (1-u-w) * s / 2
            t = t + x(9)  * (1+v-w) * s / 2
            t = t - x(9)  * (1+v-w)*(1-u-w) * s**2 / 2
            t = t + x(9)  * (1-v-w) * s / 2
            t = t + x(9)  * (1+v-w) * s / 2
            t = t - x(9)  * (1+v-w)*(1-v-w) * s**2 / 2
            t = t - x(9)  * (1-v-w)*(1-u-w) * s**2 / 2
            t = t - x(9)  * (1+v-w)*(1-u-w) * s**2 / 2
            t = t - x(9)  * (1+v-w)*(1-v-w) * s**2 / 2
            t = t + x(9)  * (1+v-w)*(1-v-w)*(1-u-w) * 2*s**3 / 2
          CASE(10)
            t = t + x(10) * w * s
            t = t - x(10) * (1-v-w) * s**2
            t = t + x(10) * w * s
            t = t - x(10) * (1-u-w) * s**2
            t = t - x(10) * (1-v-w) * s**2
            t = t - x(10) * (1-u-w) * s**2
            t = t + x(10) * (1-u-w) * (1-v-w) * 2*s**3
          CASE(11)
            t = t + x(11) * w * s
            t = t - x(11) * (1-v-w) * s**2
            t = t + x(11) * w * s
            t = t - x(11) * (1+u-w) * s**2
            t = t - x(11) * (1-v-w) * s**2
            t = t - x(11) * (1+u-w) * s**2
            t = t + x(11) * (1+u-w) * (1-v-w) * 2*s**3
          CASE(12)
            t = t + x(12) * w * s
            t = t - x(12) * (1+v-w) * s**2
            t = t + x(12) * w * s
            t = t - x(12) * (1+u-w) * s**2
            t = t - x(12) * (1+v-w) * s**2
            t = t - x(12) * (1+u-w) * s**2
            t = t + x(12) * (1+u-w) * (1+v-w) * 2*s**3
          CASE(13)
            t = t + x(13) * w*s
            t = t - x(13) * (1+v-w) * s**2
            t = t + x(13) * w*s
            t = t - x(13) * (1-u-w) * s**2
            t = t - x(13) * (1+v-w) * s**2
            t = t - x(13) * (1-u-w) * s**2
            t = t + x(13) * (1-u-w) * (1+v-w) * 2*s**3
          END SELECT
          ddx(3,3) = ddx(3,3) + t
        END DO
        ddx(2,1) = ddx(1,2)
        ddx(3,1) = ddx(1,3)
        ddx(3,2) = ddx(2,3)
        RETURN

      END IF

      DO n = 1,k
        IF ( x(n) /= 0.0d0 ) THEN
          p => elt % BasisFunctions(n) % p
          q => elt % BasisFunctions(n) % q
          r => elt % BasisFunctions(n) % r
          Coeff => elt % BasisFunctions(n) % Coeff
!------------------------------------------------------------------------------
!         @^2x/@u^2
!------------------------------------------------------------------------------
          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
             IF ( p(i) >= 2 ) THEN
                s = s + p(i) * (p(i)-1) * Coeff(i) * u**(p(i)-2) * v**q(i) * w**r(i)
             END IF
          END DO
          ddx(1,1) = ddx(1,1) + s*x(n)

!------------------------------------------------------------------------------
!         @^2x/@u@v
!------------------------------------------------------------------------------
          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
              IF (  p(i) >= 1 .AND. q(i) >= 1 ) THEN
                 s = s + p(i) * q(i) * Coeff(i) * u**(p(i)-1) * v**(q(i)-1) * w**r(i)
              END IF
          END DO
          ddx(1,2) = ddx(1,2) + s*x(n)

!------------------------------------------------------------------------------
!         @^2x/@u@w
!------------------------------------------------------------------------------
          s = 0.0d0
          DO i = 2,k
              IF (  p(i) >= 1 .AND. r(i) >= 1 ) THEN
                 s = s + p(i) * r(i) * Coeff(i) * u**(p(i)-1) * v**q(i) * w**(r(i)-1)
              END IF
          END DO
          ddx(1,3) = ddx(1,3) + s*x(n)

!------------------------------------------------------------------------------
!         @^2x/@v^2
!------------------------------------------------------------------------------
          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
             IF ( q(i) >= 2 ) THEN
                s = s + q(i) * (q(i)-1) * Coeff(i) * u**p(i) * v**(q(i)-2) * w**r(i)
             END IF
          END DO
          ddx(2,2) = ddx(2,2) + s*x(n)

!------------------------------------------------------------------------------
!         @^2x/@v@w
!------------------------------------------------------------------------------
          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
              IF (  q(i) >= 1 .AND. r(i) >= 1 ) THEN
                 s = s + q(i) * r(i) * Coeff(i) * u**p(i) * v**(q(i)-1) * w**(r(i)-1)
              END IF
          END DO
          ddx(2,3) = ddx(2,3) + s*x(n)

!------------------------------------------------------------------------------
!         @^2x/@w^2
!------------------------------------------------------------------------------
          s = 0.0d0
          DO i = 1,BasisFunctions(n) % n
             IF ( r(i) >= 2 ) THEN
                s = s + r(i) * (r(i)-1) * Coeff(i) * u**p(i) * v**q(i) * w**(r(i)-2)
             END IF
          END DO
          ddx(3,3) = ddx(3,3) + s*x(n)

        END IF
      END DO

      ddx(2,1) = ddx(1,2)
      ddx(3,1) = ddx(1,3)
      ddx(3,2) = ddx(2,3)

   END FUNCTION SecondDerivatives3D
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  Return the values of the reference element basis functions. In the case of
!>  p-element, the values of the lowest-order basis functions corresponding 
!>  to the background mesh are returned.
!------------------------------------------------------------------------------
   SUBROUTINE NodalBasisFunctions( n, Basis, element, u, v, w, USolver)
!------------------------------------------------------------------------------
     INTEGER :: n                 !< The number of (background) element nodes
     REAL(KIND=dp) :: Basis(:)    !< The values of reference element basis
     TYPE(Element_t) :: element   !< The element structure
     REAL(KIND=dp) :: u,v,w       !< The coordinates of the reference element point
     TYPE(Solver_t), POINTER, OPTIONAL :: USolver
!------------------------------------------------------------------------------
     INTEGER   :: i, q, dim, elemcode
     REAL(KIND=dp) :: NodalBasis(n)
     LOGICAL :: pElem
     
     dim = Element % TYPE % DIMENSION
     elemcode = element % Type % ElementCode
     pElem = isActivePElement( Element, USolver ) 
     
     ! Fast path for all standard (non-pyramid) elements and P-elements.
     IF( elemcode/100 /= 6 .AND. ( pelem .OR. elemcode/100 >= MODULO(elemcode,100) ) ) THEN
       SELECT CASE(elemcode/100)
       CASE( 2 )
         CALL LineNodalPBasisAll(u, Basis )
       CASE( 3 )
         IF( pElem ) THEN
           CALL TriangleNodalPBasisAll(u, v, Basis)
         ELSE
           CALL TriangleNodalLBasisAll(u, v, Basis)
         END IF
       CASE( 4 )
         CALL QuadNodalPBasisAll(u, v, Basis )
       CASE( 5 )
         IF( pElem ) THEN
           CALL TetraNodalPBasisAll(u, v, w, Basis)
         ELSE
           CALL TetraNodalLBasisAll(u, v, w, Basis)
         END IF
       CASE( 7 )
         IF( pElem ) THEN
           CALL WedgeNodalPBasisAll(u, v, w, Basis)
         ELSE
           CALL WedgeNodalLBasisAll(u, v, w, Basis)
         END IF
       CASE( 8 )
         CALL BrickNodalPBasisAll(u,v,w,Basis)
       END SELECT
       RETURN
     END IF
     
     IF ( pElem ) THEN
       SELECT CASE(elemcode / 100 )
       CASE(2)
         CALL NodalBasisFunctions1D( Basis, element, u )
       CASE(3)
         DO q=1,n
           Basis(q) = TriangleNodalPBasis(q, u, v)
         END DO
       CASE(4) 
         DO q=1,n
           Basis(q) = QuadNodalPBasis(q, u, v)
         END DO
       CASE(5)
         DO q=1,n
           Basis(q) = TetraNodalPBasis(q, u, v, w)
         END DO
       CASE(6) 
         DO q=1,n
           Basis(q) = PyramidNodalPBasis(q, u, v, w)
         END DO
       CASE(7)
         DO q=1,n
           Basis(q) = WedgeNodalPBasis(q, u, v, w)
         END DO
       CASE(8) 
         DO q=1,n
           Basis(q) = BrickNodalPBasis(q, u, v, w)
         END DO
       END SELECT
     ELSE
       SELECT CASE( dim )
       CASE(1)
         CALL NodalBasisFunctions1D( Basis, element, u )
       CASE(2)
         CALL NodalBasisFunctions2D( Basis, element, u,v )
       CASE(3)
         IF ( elemcode/100==6 ) THEN
           NodalBasis=0
           DO q=1,n
             NodalBasis(q)  = 1.0d0
             Basis(q) = InterpolateInElement3D( element, NodalBasis, u,v,w )
             NodalBasis(q)  = 0.0d0
           END DO
         ELSE
           CALL NodalBasisFunctions3D( Basis, element, u,v,w )
         END IF
       END SELECT
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE NodalBasisFunctions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  Return the gradient of the reference element basis functions, with the
!>  gradient taken with respect to the reference element coordinates. In the case
!>  of p-element, the gradients of the lowest-order basis functions corresponding 
!>  to the background mesh are returned.
!------------------------------------------------------------------------------
   SUBROUTINE NodalFirstDerivatives( n, dLBasisdx, element, u, v, w, USolver )
!------------------------------------------------------------------------------
     INTEGER :: n                    !< The number of (background) element nodes
     REAL(KIND=dp) :: dLBasisdx(:,:) !< The gradient of reference element basis functions
     TYPE(Element_t) :: element      !< The element structure
     REAL(KIND=dp) :: u,v,w          !< The coordinates of the reference element point
     TYPE(Solver_t), POINTER, OPTIONAL :: USolver
!------------------------------------------------------------------------------
     INTEGER   :: i, q, dim, elemcode
     REAL(KIND=dp) :: NodalBasis(n)
     LOGICAL :: pElem
!------------------------------------------------------------------------------
     dim = Element % TYPE % DIMENSION
     elemcode = element % TYPE % ElementCode
     pElem = isActivePElement( Element, USolver ) 
     
     ! Fast path for all standard (non-pyramid) elements and P-elements.
     IF( elemcode/100 /= 6 .AND. ( pelem .OR. elemcode/100 >= MODULO(elemcode,100) ) ) THEN
       SELECT CASE(elemcode/100)
       CASE( 2 )
         CALL dLineNodalPBasisAll(u, dLBasisdx )
       CASE( 3 )
         IF( pElem ) THEN
           CALL dTriangleNodalPBasisAll(u, v, dLBasisdx)
         ELSE
           CALL dTriangleNodalLBasisAll(u, v, dLBasisdx)
         END IF
       CASE( 4 )
         CALL dQuadNodalPBasisAll(u, v, dLBasisdx )
       CASE( 5 )
         IF( pElem ) THEN
           CALL dTetraNodalPBasisAll(u, v, w, dLBasisdx)
         ELSE
           CALL dTetraNodalLBasisAll(u, v, w, dLBasisdx)
         END IF
       CASE( 7 )
         IF( pElem ) THEN
           CALL dWedgeNodalPBasisAll(u, v, w, dLBasisdx)
         ELSE
           CALL dWedgeNodalLBasisAll(u, v, w, dLBasisdx)
         END IF
       CASE( 8 )
         CALL dBrickNodalPBasisAll(u,v,w,dLBasisdx)
       END SELECT
       RETURN
     END IF
     
     IF ( IsActivePElement(Element, USolver ) ) THEN
       SELECT CASE(elemcode / 100 )
       CASE(2)
         CALL NodalFirstDerivatives1D( dLBasisdx, element, u )
       CASE(3)
         DO q=1,n
           dLBasisdx(q,1:2) = dTriangleNodalPBasis(q, u, v)
         END DO
       CASE(4)
         DO q=1,n
           dLBasisdx(q,1:2) = dQuadNodalPBasis(q, u, v)
         END DO
       CASE(5)
         DO q=1,n
           dLBasisdx(q,1:3) = dTetraNodalPBasis(q, u, v, w)
         END DO
       CASE( 6 )
         DO q=1,n
           dLBasisdx(q,1:3) = dPyramidNodalPBasis(q, u, v, w)
         END DO
       CASE( 7 ) 
         DO q=1,n
           dLBasisdx(q,1:3) = dWedgeNodalPBasis(q, u, v, w)
         END DO
       CASE( 8 ) 
         DO q=1,n
           dLBasisdx(q,1:3) = dBrickNodalPBasis(q, u, v, w)
         END DO
       END SELECT
     ELSE
       SELECT CASE(dim)
       CASE(1)
         CALL NodalFirstDerivatives1D( dLBasisdx, element, u )
       CASE(2)
         CALL NodalFirstDerivatives2D( dLBasisdx, element, u,v )
       CASE(3)
         IF ( elemcode / 100 == 6 ) THEN
           NodalBasis=0
           DO q=1,n
             NodalBasis(q)  = 1.0d0
             dLBasisdx(q,1) = FirstDerivativeInU3D(element,NodalBasis,u,v,w)
             dLBasisdx(q,2) = FirstDerivativeInV3D(element,NodalBasis,u,v,w)
             dLBasisdx(q,3) = FirstDerivativeInW3D(element,NodalBasis,u,v,w)
             NodalBasis(q)  = 0.0d0
           END DO
         ELSE
           CALL NodalFirstDerivatives3D( dLBasisdx, element, u,v,w )
         END IF
       END SELECT
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE NodalFirstDerivatives
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Return basis function degrees
!------------------------------------------------------------------------------
   SUBROUTINE ElementBasisDegree( Element, BasisDegree, USolver )
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Element_t), TARGET :: Element   !< Element structure
     INTEGER :: BasisDegree(:)            !< Degree of each basis function in Basis(:) vector. 
     TYPE(Solver_t), TARGET, OPTIONAL :: USolver
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: t,s
     LOGICAL :: invert, degrees
     INTEGER :: i, j, k, l, q, p, f, n, nb, dim, cdim, locali, localj,  &
          tmp(4), direction(4), BDOFs, BodyId

     TYPE(Solver_t), POINTER :: pSolver

     LOGICAL :: SerendipityPBasis

     TYPE(Element_t) :: Bubble
     TYPE(Element_t), POINTER :: Edge, Face, Parent
!------------------------------------------------------------------------------

     IF( PRESENT( USolver ) ) THEN
       pSolver => USolver
     ELSE
       pSolver => CurrentModel % Solver
     END IF

     n    = Element % TYPE % NumberOfNodes
     dim  = Element % TYPE % DIMENSION
     cdim = CoordinateSystemDimension()


     BasisDegree = 0
     BasisDegree(1:n) = Element % Type % BasisFunctionDegree

     IF ( isActivePElement(element) ) THEN

       BodyId = Element % BodyId
       IF (BodyId==0 .AND. ASSOCIATED(Element % BoundaryInfo)) THEN
         Parent => Element % PDefs % LocalParent         
         IF(ASSOCIATED(Parent)) BodyId = Parent % BodyId
         IF( BodyId == 0 ) THEN
           Parent => Element % BoundaryInfo % Left
           IF( ASSOCIATED(Parent)) BodyId = Parent % BodyId
         END IF
         IF(BodyId == 0) THEN
           Parent => Element % BoundaryInfo % Right
           IF( ASSOCIATED(Parent)) BodyId = Parent % BodyId
         END IF
       END IF
       IF (BodyId==0) THEN
         CALL Warn('ElementBasisDegree', 'Element '//I2S(Element % ElementIndex)//' of type '//&
             I2S(Element % TYPE % ElementCode)//' has 0 BodyId, assuming index 1')
         BodyId = 1
       END IF

       ! Check for need of P basis degrees and set degree of
       ! linear basis if vector asked:
       ! ---------------------------------------------------
       BasisDegree(1:n) = 1
       q = n

       SerendipityPBasis = Element % PDefs % Serendipity
!------------------------------------------------------------------------------
     SELECT CASE( Element % TYPE % ElementCode ) 
!------------------------------------------------------------------------------

     ! P element code for line element:
     ! --------------------------------
     CASE(202)
        ! Bubbles of line element
        p = pSolver % Def_Dofs(2,BodyId,6)
        nb = pSolver % Def_Dofs(2,BodyId,5)
        BDOFs = MAX(GetBubbleDOFs(Element, p), nb)

        IF (BDOFs > 0) THEN
           p = getEffectiveBubbleP(element,p,bdofs)
           ! For each bubble in line element get value of basis function
           DO i=1, BDOFs
             IF (q >= SIZE(BasisDegree)) CYCLE
             q = q + 1
             BasisDegree(q) = 1+i
           END DO
        END IF

!------------------------------------------------------------------------------
! P element code for edges and bubbles of triangle
     CASE(303)
        ! Edges of triangle
        IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge calculate the value of edge basis function
           DO i=1,3
              Edge => CurrentModel % Solver % Mesh % Edges( Element % EdgeIndexes(i) )

              ! For each dof in edge get value of p basis function 
              DO k=1,Edge % BDOFs
                 IF (q >= SIZE(BasisDegree)) CYCLE
                 q = q + 1
                 BasisDegree(q) = 1+k
              END DO
           END DO 
        END IF

        ! Bubbles of p triangle      

        p = pSolver % Def_Dofs(3,BodyId,6)
        nb = pSolver % Def_Dofs(3,BodyId,5)
        BDOFs = MAX(GetBubbleDOFs(Element, p), nb)
        IF ( BDOFs > 0 ) THEN
           ! Get element p
           p = getEffectiveBubbleP(element,p,bdofs)

           DO i = 0,p-3
              DO j = 0,p-i-3
                 IF ( q >= SIZE(BasisDegree) ) CYCLE
                 q = q + 1
                 BasisDegree(q) = 3+i+j
              END DO
           END DO
        END IF
!------------------------------------------------------------------------------
! P element code for quadrilateral edges and bubbles 
     CASE(404)
        ! Edges of p quadrilateral
        IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge begin node calculate values of edge functions 
           DO i=1,4
              Edge => CurrentModel % Solver % Mesh % Edges( Element % EdgeIndexes(i) )
              ! For each DOF in edge calculate value of p basis function
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(BasisDegree) ) CYCLE
                 q = q + 1
                 BasisDegree(q) = 1+k
              END DO              
           END DO         
        END IF

        ! Bubbles of p quadrilateral
        p = pSolver % Def_Dofs(4,BodyId,6)
        nb = pSolver % Def_Dofs(4,BodyId,5)
        BDOFs = MAX(GetBubbleDOFs(Element, p), nb)
        IF ( BDOFs > 0 ) THEN
          ! Get element P
           p = getEffectiveBubbleP(element,p,bdofs)
          
           IF(SerendipityPBasis) THEN
             DO i=2,(p-2)
               DO j=2,(p-i)
                  IF ( q >= SIZE(BasisDegree) ) CYCLE
                  q = q + 1
                  BasisDegree(q) = i+j
               END DO
             END DO
           ELSE
             DO i=0,p-2
                DO j=0,p-2
                   IF ( q >= SIZE(BasisDegree) ) CYCLE
                   q = q + 1
                   BasisDegree(q) = 2+i+j
                END DO
             END DO
           END IF
        END IF
!------------------------------------------------------------------------------
! P element code for tetrahedron edges, faces and bubbles
     CASE(504) 
        ! Edges of p tetrahedron
        IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN   
           ! For each edge calculate value of edge functions
           DO i=1,6
              Edge => CurrentModel % Solver % Mesh % Edges (Element % EdgeIndexes(i))

              ! Do not solve edge DOFS if there is not any
              IF (Edge % BDOFs <= 0) CYCLE

              ! For each DOF in edge calculate value of edge functions 
              ! and their derivatives for edge=i, i=k+1
              DO k=1, Edge % BDOFs
                 IF (q >= SIZE(BasisDegree)) CYCLE
                 q = q + 1
                 BasisDegree(q) = 1+k
              END DO
           END DO
        END IF

        ! Faces of p tetrahedron
        IF ( ASSOCIATED( Element % FaceIndexes )) THEN
           ! For each face calculate value of face functions
           DO F=1,4
              Face => CurrentModel % Solver % Mesh % Faces (Element % FaceIndexes(F))

              ! Do not solve face DOFs if there is not any
              IF (Face % BDOFs <= 0) CYCLE

              ! Get face p 
              p = Face % PDefs % P

              ! For each DOF in face calculate value of face functions and 
              ! their derivatives for face=F and index pairs 
              ! i,j=0,..,p-3, i+j=0,..,p-3
              DO i=0,p-3
                 DO j=0,p-i-3
                    IF (q >= SIZE(BasisDegree)) CYCLE
                    q = q + 1 
                    BasisDegree(q) = 3+i+j
                 END DO
              END DO
           END DO
        END IF

        ! Bubbles of p tetrahedron
        p = pSolver % Def_Dofs(5,BodyId,6)
        nb = pSolver % Def_Dofs(5,BodyId,5)
        BDOFs = MAX(GetBubbleDOFs(Element, p), nb)
        IF ( BDOFs > 0 ) THEN
           p = getEffectiveBubbleP(element,p,bdofs)

           DO i=0,p-4
              DO j=0,p-i-4
                 DO k=0,p-i-j-4
                    IF (q >= SIZE(BasisDegree)) CYCLE
                    q = q + 1
                    BasisDegree(q) = 4+i+j+k
                 END DO
              END DO
           END DO
           
        END IF
!------------------------------------------------------------------------------
! P element code for pyramid edges, faces and bubbles
     CASE(605)
        ! Edges of P Pyramid
        IF (ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge in wedge, calculate values of edge functions
           DO i=1,8
              Edge => CurrentModel % Solver % Mesh % Edges( Element % EdgeIndexes(i) )

              ! Do not solve edge dofs, if there is not any
              IF (Edge % BDOFs <= 0) CYCLE
              
              ! For each DOF in edge calculate values of edge functions
              ! and their derivatives for edge=i and i=k+1
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(BasisDegree) ) CYCLE
                 q = q + 1
                 BasisDegree(q) = 1+k
              END DO
           END DO
        END IF
        
        ! Faces of P Pyramid
        IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
           ! For each face in pyramid, calculate values of face functions
           DO F=1,5
              Face => CurrentModel % Solver % Mesh % Faces( Element % FaceIndexes(F) )

              ! Do not solve face dofs, if there is not any
              IF ( Face % BDOFs <= 0) CYCLE
              
              ! Get face p
              p = Face % PDefs % P 
              
              ! Handle triangle and square faces separately
              SELECT CASE(F)
              CASE (1)
                 ! For each face calculate values of functions from index
                 ! pairs i,j=2,..,p-2 i+j=4,..,p
!                DO i=2,p-2
!                   DO j=2,p-i
                 DO i=0,p-2
                    DO j=0,p-2
                       IF ( q >= SIZE(BasisDegree) ) CYCLE
                       q = q + 1
                       BasisDegree(q) = 2+i+j
                    END DO
                 END DO

              CASE (2,3,4,5)
                 ! For each face calculate values of functions from index
                 ! pairs i,j=0,..,p-3 i+j=0,..,p-3
                 DO i=0,p-3
                    DO j=0,p-i-3
                       IF ( q >= SIZE(BasisDegree) ) CYCLE
                       q = q + 1
                       BasisDegree(q) = 3+i+j
                    END DO
                 END DO
              END SELECT    
           END DO
        END IF

        ! Bubbles of P Pyramid
        p = pSolver % Def_Dofs(6,BodyId,6)
        nb = pSolver % Def_Dofs(6,BodyId,5)
        BDOFs = MAX(GetBubbleDOFs(Element, p), nb)
        IF (BDOFs > 0) THEN 
           ! Get element p
           p = getEffectiveBubbleP(element,p,bdofs)

           ! Calculate value of bubble functions from indexes
           ! i,j,k=0,..,p-3 i+j+k=0,..,p-3
           DO i=0,p-3
              DO j=0,p-i-3
                 DO k=0,p-i-j-3
                    IF ( q >= SIZE(BasisDegree)) CYCLE
                    q = q + 1
                    BasisDegree(q) = 3+i+j+k
                 END DO
              END DO
           END DO
        END IF
        
!------------------------------------------------------------------------------
! P element code for wedge edges, faces and bubbles
     CASE(706)
        ! Edges of P Wedge
        IF (ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge in wedge, calculate values of edge functions
           DO i=1,9
              Edge => CurrentModel % Solver % Mesh % Edges( Element % EdgeIndexes(i) )

              ! Do not solve edge dofs, if there is not any
              IF (Edge % BDOFs <= 0) CYCLE
              
              ! For each DOF in edge calculate values of edge functions
              ! and their derivatives for edge=i and i=k+1
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(BasisDegree) ) CYCLE
                 q = q + 1
                 BasisDegree(q) = 1+k
              END DO
           END DO
        END IF

        ! Faces of P Wedge 
        IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
           ! For each face in wedge, calculate values of face functions
           DO F=1,5
              Face => CurrentModel % Solver % Mesh % Faces( Element % FaceIndexes(F) )

              ! Do not solve face dofs, if there is not any
              IF ( Face % BDOFs <= 0) CYCLE

              p = Face % PDefs % P 
              
              ! Handle triangle and square faces separately
              SELECT CASE(F)
              CASE (1,2)
                 ! For each face calculate values of functions from index
                 ! pairs i,j=0,..,p-3 i+j=0,..,p-3
                 DO i=0,p-3
                    DO j=0,p-i-3
                       IF ( q >= SIZE(BasisDegree) ) CYCLE
                       q = q + 1
                       BasisDegree(q) = 3+i+j
                    END DO
                 END DO
              CASE (3,4,5)
                 ! For each face calculate values of functions from index
                 ! pairs i,j=2,..,p-2 i+j=4,..,p
                 IF(SerendipityPBasis) THEN
                   DO i=2,p-2
                      DO j=2,p-i
                         IF ( q >= SIZE(BasisDegree) ) CYCLE
                         q = q + 1
                         BasisDegree(q) = i+j
                      END DO
                   END DO
                 ELSE
                   DO i=0,p-2
                      DO j=0,p-2
                         IF ( q >= SIZE(BasisDegree) ) CYCLE
                         q = q + 1
                         BasisDegree(q) = 2+i+j
                      END DO
                   END DO
                 END IF
              END SELECT
                           
           END DO
        END IF

        ! Bubbles of P Wedge
        p = pSolver % Def_Dofs(7,BodyId,6)
        nb = pSolver % Def_Dofs(7,BodyId,5)
        BDOFs = MAX(GetBubbleDOFs(Element, p), nb)
        IF ( BDOFs > 0 ) THEN
           ! Get p from element
           p = getEffectiveBubbleP(element,p,bdofs)

           ! For each bubble calculate value of basis function and its derivative
           ! for index pairs i,j=0,..,p-5 k=2,..,p-3 i+j+k=2,..,p-3
           IF(SerendipityPBasis) THEN
             DO i=0,p-5
                DO j=0,p-5-i
                   DO k=2,p-3-i-j
                      IF ( q >= SIZE(BasisDegree) ) CYCLE
                      q = q + 1
                      BasisDegree(q) = 3+i+j+k
                   END DO
                END DO
             END DO
           ELSE
             DO i=0,p-3
                DO j=0,p-i-3
                   DO k=0,p-2
                      IF ( q >= SIZE(BasisDegree) ) CYCLE
                      q = q + 1
                      BasisDegree(q) = 3+i+j+k
                   END DO
                END DO
             END DO
           END IF
        END IF

!------------------------------------------------------------------------------
! P element code for brick edges, faces and bubbles
     CASE(808) 
        ! Edges of P brick
        IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge in brick, calculate values of edge functions 
           DO i=1,12
              Edge => CurrentModel % Solver % Mesh % Edges( Element % EdgeIndexes(i) )

              ! Do not solve edge dofs, if there is not any
              IF (Edge % BDOFs <= 0) CYCLE
              
              ! For each DOF in edge calculate values of edge functions
              ! and their derivatives for edge=i and i=k+1
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(BasisDegree) ) CYCLE
                 q = q + 1
                 BasisDegree(q) = 1+k
              END DO
           END DO 
        END IF

        ! Faces of P brick
        IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
           ! For each face in brick, calculate values of face functions
           DO F=1,6
              Face => CurrentModel % Solver % Mesh % Faces( Element % FaceIndexes(F) )
                          
              ! Do not calculate face values if no dofs
              IF (Face % BDOFs <= 0) CYCLE
              
              ! Get p for face
              p = Face % PDefs % P

              ! For each face calculate values of functions from index
              ! pairs i,j=2,..,p-2 i+j=4,..,p
              IF(SerendipityPBasis) THEN
                DO i=2,p-2
                   DO j=2,p-i
                      IF ( q >= SIZE(BasisDegree) ) CYCLE
                      q = q + 1
                      BasisDegree(q) = i+j
                   END DO
                 END DO
              ELSE
                DO i=0,p-2
                   DO j=0,p-2
                      IF ( q >= SIZE(BasisDegree) ) CYCLE
                      q = q + 1
                      BasisDegree(q) = 2+i+j
                   END DO
                END DO
              END IF
           END DO
        END IF

        ! Bubbles of p brick
        p = pSolver % Def_Dofs(7,BodyId,6)
        nb = pSolver % Def_Dofs(7,BodyId,5)
        BDOFs = MAX(GetBubbleDOFs(Element, p), nb)
        IF ( BDOFs > 0 ) THEN
           ! Get p from bubble DOFs 
           p = getEffectiveBubbleP(element,p,bdofs)

           ! For each bubble calculate value of basis function and its derivative
           ! for index pairs i,j,k=2,..,p-4, i+j+k=6,..,p
           IF(SerendipityPBasis) THEN
             DO i=2,p-4
               DO j=2,p-i-2
                 DO k=2,p-i-j
                   IF ( q >= SIZE(BasisDegree) ) CYCLE
                   q = q + 1
                   BasisDegree(q) = i+j+k
                 END DO
               END DO
             END DO
           ELSE
             DO i=0,p-2
                DO j=0,p-2
                   DO k=0,p-2
                      IF ( q >= SIZE(BasisDegree) ) CYCLE
                      q = q + 1
                      BasisDegree(q) = 2+i+j+k
                   END DO
                END DO
             END DO
           END IF
        END IF

     END SELECT
     END IF ! P element flag check
!------------------------------------------------------------------------------
   END SUBROUTINE ElementBasisDegree
     FUNCTION CrossProduct( v1, v2 ) RESULT( v3 )
!------------------------------------------------------------------------------
       IMPLICIT NONE
       REAL(KIND=dp) :: v1(3), v2(3), v3(3)
       v3(1) =  v1(2)*v2(3) - v1(3)*v2(2)
       v3(2) = -v1(1)*v2(3) + v1(3)*v2(1)
       v3(3) =  v1(1)*v2(2) - v1(2)*v2(1)
!------------------------------------------------------------------------------
     END FUNCTION CrossProduct
   FUNCTION InterpolateInElement( elm,f,u,v,w,Basis ) RESULT(val)
!------------------------------------------------------------------------------
!
!  DESCRIPTION:
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: f(:)
!     INPUT: Nodal values of the quantity whose value we want to know
!
!    REAL(KIND=dp) :: u,v,w
!     INPUT: Point at which to evaluate the value
!
!    REAL(KIND=dp), OPTIONAL :: Basis(:)
!      INPUT: Values of the basis functions at the point u,v,w can be given here,
!      if known, otherwise the will be computed from the definition
!                 
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: y
!       value of the quantity y = x(u,v,w)
!    
!------------------------------------------------------------------------------

     TYPE(Element_t) :: elm
     REAL(KIND=dp) :: u,v,w
     REAL(KIND=dp) :: f(:)
     REAL(KIND=dp), OPTIONAL :: Basis(:)

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: val
     INTEGER :: n

     IF ( PRESENT( Basis ) ) THEN
!------------------------------------------------------------------------------
!      Basis function values given, just sum the result ...
!------------------------------------------------------------------------------
       n = elm % TYPE % NumberOfNodes
       val = SUM( f(1:n)*Basis(1:n) )
     ELSE
!------------------------------------------------------------------------------
!      ... otherwise compute from the definition.
!------------------------------------------------------------------------------
       SELECT CASE (elm % TYPE % DIMENSION)
         CASE (0)
           val = f(1)
         CASE (1)
           val = InterpolateInElement1D( elm,f,u )
         CASE (2)
           val = InterpolateInElement2D( elm,f,u,v )
         CASE (3)
           val = InterpolateInElement3D( elm,f,u,v,w )
       END SELECT
     END IF
  
   END FUNCTION InterpolateInElement
 FUNCTION GetEdgeMap( ElementFamily ) RESULT(EdgeMap)
!------------------------------------------------------------------------------
    INTEGER :: ElementFamily
    INTEGER, POINTER :: EdgeMap(:,:)

    INTEGER, TARGET :: Point(1,1)
    INTEGER, TARGET :: Line(1,2)
    INTEGER, TARGET :: Triangle(3,2)
    INTEGER, TARGET :: Quad(4,2)
    INTEGER, TARGET :: Tetra(6,2)
    INTEGER, TARGET :: Pyramid(8,2)
    INTEGER, TARGET :: Wedge(9,2)
    INTEGER, TARGET :: Brick(12,2)

    LOGICAL :: Initialized(8) = .FALSE.
  
    SAVE Line, Triangle, Wedge, Brick, Tetra, Quad, Pyramid, Initialized
    !$OMP THREADPRIVATE(Line, Triangle, Wedge, Brick, Tetra, Quad, Pyramid, Initialized)

    SELECT CASE(ElementFamily)
    CASE(1)
      EdgeMap => Point
    CASE(2)
      EdgeMap => Line
    CASE(3)
      EdgeMap => Triangle
    CASE(4) 
      EdgeMap => Quad
    CASE(5) 
      EdgeMap => Tetra
    CASE(6) 
      EdgeMap => Pyramid
    CASE(7) 
      EdgeMap => Wedge
    CASE(8) 
      EdgeMap => Brick
    CASE DEFAULT
      WRITE( Message,'(A,I0,A)') 'Element family ',ElementFamily,' is not known!'
      CALL Fatal( 'GetEdgeMap', Message )
    END SELECT
 
    IF ( .NOT. Initialized(ElementFamily) ) THEN
       Initialized(ElementFamily) = .TRUE.
       SELECT CASE(ElementFamily)
       CASE(1)
         EdgeMap(1,1) = 1

       CASE(2)
         EdgeMap(1,:) = [ 1,2 ]

       CASE(3)
         EdgeMap(1,:) = [ 1,2 ]
         EdgeMap(2,:) = [ 2,3 ]
         EdgeMap(3,:) = [ 3,1 ]

       CASE(4)
         EdgeMap(1,:) = [ 1,2 ]
         EdgeMap(2,:) = [ 2,3 ]
         EdgeMap(3,:) = [ 3,4 ]
         EdgeMap(4,:) = [ 4,1 ]

       CASE(5)
         EdgeMap(1,:) = [ 1,2 ]
         EdgeMap(2,:) = [ 2,3 ]
         EdgeMap(3,:) = [ 3,1 ]
         EdgeMap(4,:) = [ 1,4 ]
         EdgeMap(5,:) = [ 2,4 ]
         EdgeMap(6,:) = [ 3,4 ]

       CASE(6)
         EdgeMap(1,:) = [ 1,2 ]
         EdgeMap(2,:) = [ 2,3 ]
         EdgeMap(3,:) = [ 4,3 ]
         EdgeMap(4,:) = [ 1,4 ]
         EdgeMap(5,:) = [ 1,5 ]
         EdgeMap(6,:) = [ 2,5 ]
         EdgeMap(7,:) = [ 3,5 ]
         EdgeMap(8,:) = [ 4,5 ]
 
       CASE(7)
         EdgeMap(1,:) = [ 1,2 ]
         EdgeMap(2,:) = [ 2,3 ]
         EdgeMap(3,:) = [ 3,1 ]
         EdgeMap(4,:) = [ 4,5 ]
         EdgeMap(5,:) = [ 5,6 ]
         EdgeMap(6,:) = [ 6,4 ]
         EdgeMap(7,:) = [ 1,4 ]
         EdgeMap(8,:) = [ 2,5 ]
         EdgeMap(9,:) = [ 3,6 ]

       CASE(8)
         EdgeMap(1,:)  = [ 1,2 ]
         EdgeMap(2,:)  = [ 2,3 ]
         EdgeMap(3,:)  = [ 4,3 ]
         EdgeMap(4,:)  = [ 1,4 ]
         EdgeMap(5,:)  = [ 5,6 ]
         EdgeMap(6,:)  = [ 6,7 ]
         EdgeMap(7,:)  = [ 8,7 ]
         EdgeMap(8,:)  = [ 5,8 ]
         EdgeMap(9,:)  = [ 1,5 ]
         EdgeMap(10,:) = [ 2,6 ]
         EdgeMap(11,:) = [ 3,7 ]
         EdgeMap(12,:) = [ 4,8 ]
       END SELECT
     END IF
!------------------------------------------------------------------------------
  END FUNCTION GetEdgeMap
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Figure out element diameter parameter for stabilization.
!------------------------------------------------------------------------------
   FUNCTION ElementDiameter( elm, nodes, UseLongEdge ) RESULT(hK)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: elm  !< element structure
     TYPE(Nodes_t) :: nodes  !< Nodal coordinate arrays of the element
     LOGICAL, OPTIONAL :: UseLongEdge  !< Use the longest edge to determine the diameter.
     REAL(KIND=dp) :: hK     !< hK
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp), DIMENSION(:), POINTER :: X,Y,Z
     INTEGER :: i,j,k,Family
     INTEGER, POINTER :: EdgeMap(:,:)
     REAL(KIND=dp) :: x0,y0,z0,A,S,CX,CY,CZ
     REAL(KIND=dp) :: J11,J12,J13,J21,J22,J23,G11,G12,G21,G22
     LOGICAL :: LongEdge=.FALSE.
!------------------------------------------------------------------------------

     IF(PRESENT(UseLongEdge)) LongEdge = UseLongEdge

     X => Nodes % x
     Y => Nodes % y
     Z => Nodes % z

     Family = Elm % TYPE % ElementCode / 100
     SELECT CASE( Family )

       CASE(1)
         hK = 0.0d0

!------------------------------------------------------------------------------
!       Triangular element
!------------------------------------------------------------------------------
       CASE(3) 
         J11 = X(2) - X(1)
         J12 = Y(2) - Y(1)
         J13 = Z(2) - Z(1)
         J21 = X(3) - X(1)
         J22 = Y(3) - Y(1)
         J23 = Z(3) - Z(1)
         G11 = J11**2  + J12**2  + J13**2
         G12 = J11*J21 + J12*J22 + J13*J23
         G22 = J21**2  + J22**2  + J23**2
         A = SQRT(G11*G22 - G12**2) / 2.0d0

         CX = ( X(1) + X(2) + X(3) ) / 3.0d0
         CY = ( Y(1) + Y(2) + Y(3) ) / 3.0d0
         CZ = ( Z(1) + Z(2) + Z(3) ) / 3.0d0

         s =     (X(1)-CX)**2 + (Y(1)-CY)**2 + (Z(1)-CZ)**2
         s = s + (X(2)-CX)**2 + (Y(2)-CY)**2 + (Z(2)-CZ)**2
         s = s + (X(3)-CX)**2 + (Y(3)-CY)**2 + (Z(3)-CZ)**2

         hK = 16.0d0*A*A / ( 3.0d0 * s )

!------------------------------------------------------------------------------
!      Quadrilateral
!------------------------------------------------------------------------------
       CASE(4)
          CX = (X(2)-X(1))**2 + (Y(2)-Y(1))**2 + (Z(2)-Z(1))**2
          CY = (X(4)-X(1))**2 + (Y(4)-Y(1))**2 + (Z(4)-Z(1))**2
          hk = 2*CX*CY/(CX+CY)

       CASE DEFAULT
         EdgeMap => GetEdgeMap(Family)

         IF(LongEdge) THEN
           hK = -1.0 * HUGE(1.0_dp)
         ELSE
           hK = HUGE(1.0_dp)
         END IF

         DO i=1,SIZE(EdgeMap,1)
           j=EdgeMap(i,1)
           k=EdgeMap(i,2)
           x0 = X(j) - X(k)
           y0 = Y(j) - Y(k)
           z0 = Z(j) - Z(k)
           IF(LongEdge) THEN
             hk = MAX(hK, x0**2 + y0**2 + z0**2)
           ELSE
             hk = MIN(hK, x0**2 + y0**2 + z0**2)
           END IF
         END DO
     END SELECT

     hK = SQRT( hK )
!------------------------------------------------------------------------------
  END FUNCTION ElementDiameter

END MODULE ElementBasis
