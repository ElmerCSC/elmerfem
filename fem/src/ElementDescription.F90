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

MODULE ElementDescription
   USE Integration
   USE GeneralUtils
   USE LinearAlgebra
   USE CoordinateSystems
   ! Use module P element basis functions 
   USE PElementMaps
   USE PElementBase
   ! Vectorized P element basis functions
   USE H1Basis
   USE Lists
   
   IMPLICIT NONE

   INTEGER, PARAMETER,PRIVATE  :: MaxDeg  = 4, MaxDeg3 = MaxDeg**3, &
                           MaxDeg2 = MaxDeg**2

   INTEGER, PARAMETER :: MAX_ELEMENT_NODES = 256

   !
   ! Module global variables
   !
   LOGICAL, PRIVATE :: TypeListInitialized = .FALSE.
   TYPE(ElementType_t), PRIVATE, POINTER :: ElementTypeList
   ! Local workspace for basis function values and mapping
!    REAL(KIND=dp), ALLOCATABLE, PRIVATE :: BasisWrk(:,:), dBasisdxWrk(:,:,:), &
!            LtoGMapsWrk(:,:,:), DetJWrk(:), uWrk(:), vWrk(:), wWrk(:)
!     !$OMP THREADPRIVATE(BasisWrk, dBasisdxWrk, LtoGMapsWrk, DetJWrk, uWrk, vWrk, wWrk)
! !DIR$ ATTRIBUTES ALIGN:64::BasisWrk, dBasisdxWrk
! !DIR$ ATTRIBUTES ALIGN:64::LtoGMapsWrk
! !DIR$ ATTRIBUTES ALIGN:64::DetJWrk
! !DIR$ ATTRIBUTES ALIGN:64::uWrk, vWrk, wWrk

CONTAINS


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
        CASE(3) 
           Element % NumberOfEdges = 3
        CASE(4) 
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
!>    Read the element description input file and add the element types to a
!>    global list. The file is assumed to be found under the name
!>        $ELMER_HOME/lib/elements.def
!>   This is the first routine the user of the element utilities should call
!>   in his/her code.
!------------------------------------------------------------------------------
   SUBROUTINE InitializeElementDescriptions()
!------------------------------------------------------------------------------
!     Local variables
!------------------------------------------------------------------------------
      CHARACTER(LEN=:), ALLOCATABLE :: str
      CHARACTER(LEN=MAX_STRING_LEN) :: tstr,elmer_home

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
      NULLIFY( element % NodeU )
      NULLIFY( element % NodeV )
      NULLIFY( element % NodeW )
      DO k=3,64
        element % NumberOfNodes = k
        element % ElementCode = 100 + k
        CALL AddElementDescription( element,BasisTerms )
      END DO

      ! then the rest of them....
      !--------------------------
#ifdef USE_ISO_C_BINDINGS
      tstr = 'ELMER_LIB'
#else
      tstr = 'ELMER_LIB'//CHAR(0)
#endif
      CALL envir( tstr,elmer_home,k ) 
      
      fexist = .FALSE.
      IF (  k > 0 ) THEN
         WRITE( tstr, '(a,a)' ) elmer_home(1:k),'/elements.def'
	 INQUIRE(FILE=TRIM(tstr), EXIST=fexist)
      END IF
      IF (.NOT. fexist) THEN
#ifdef USE_ISO_C_BINDINGS
        tstr = 'ELMER_HOME'
#else
        tstr = 'ELMER_HOME'//CHAR(0)
#endif
        CALL envir( tstr,elmer_home,k ) 
        IF ( k > 0 ) THEN
           WRITE( tstr, '(a,a)' ) elmer_home(1:k),&
'/share/elmersolver/lib/elements.def'
           INQUIRE(FILE=TRIM(tstr), EXIST=fexist)
        END IF
        IF ((.NOT. fexist) .AND. k > 0) THEN
           WRITE( tstr, '(a,a)' ) elmer_home(1:k),&
                '/elements.def'
           INQUIRE(FILE=TRIM(tstr), EXIST=fexist)
        END IF
     END IF
     IF (.NOT. fexist) THEN
        CALL GetSolverHome(elmer_home, n)
        WRITE(tstr, '(a,a)') elmer_home(1:n), &
                             '/lib/elements.def'
        INQUIRE(FILE=TRIM(tstr), EXIST=fexist)
     END IF
     IF (.NOT. fexist) THEN
        CALL Fatal('InitializeElementDescriptions', &
             'elements.def not found')
     END IF

      OPEN( 1,FILE=TRIM(tstr), STATUS='OLD' )

      ALLOCATE(CHARACTER(MAX_STRING_LEN)::str)
      DO WHILE( ReadAndTrim(1,str) )

        IF ( SEQL(str, 'element') ) THEN

          BasisTerms = 0

          NULLIFY( element % NodeU )
          NULLIFY( element % NodeV )
          NULLIFY( element % NodeW )

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
            IF ( .NOT.ASSOCIATED( element % NodeV ) ) THEN
              ALLOCATE( element % NodeV(element % NumberOfNodes) )
              element % NodeV = 0.0d0
            END IF

            IF ( .NOT.ASSOCIATED( element % NodeW ) ) THEN
              ALLOCATE( element % NodeW(element % NumberOfNodes) )
              element % NodeW = 0.0d0
            END IF

            CALL AddElementDescription( element,BasisTerms )
          ELSE
            IF ( ASSOCIATED( element % NodeU ) ) DEALLOCATE( element % NodeU )
            IF ( ASSOCIATED( element % NodeV ) ) DEALLOCATE( element % NodeV )
            IF ( ASSOCIATED( element % NodeW ) ) DEALLOCATE( element % NodeW )
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
   FUNCTION GetElementType( code,CompStabFlag ) RESULT(element)
!------------------------------------------------------------------------------
      INTEGER :: code
      LOGICAL, OPTIONAL :: CompStabFlag
      TYPE(ElementType_t), POINTER :: element
!------------------------------------------------------------------------------
!     Local variables
!------------------------------------------------------------------------------
      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t), POINTER :: Elm
!------------------------------------------------------------------------------
      element => ElementTypeList

      DO WHILE( ASSOCIATED(element) )
        IF ( code == element % ElementCode ) EXIT
        element => element % NextElementType
      END DO

      IF ( .NOT. ASSOCIATED( element ) ) THEN
        WRITE( message, * ) &
             'Element type code ',code,' not found. Ignoring element.'
        CALL Warn( 'GetElementType', message )
        RETURN
      END IF

      IF ( PRESENT( CompStabFlag ) ) THEN
        IF ( .NOT. CompStabFlag ) RETURN
      END IF

      IF ( Element % StabilizationMK == 0.0d0 ) THEN
        ALLOCATE( Elm )
        Elm % TYPE => element
        Elm % BDOFs  = 0
        Elm % DGDOFs = 0
        NULLIFY( Elm % PDefs )
        NULLIFY( Elm % DGIndexes )
        NULLIFY( Elm % EdgeIndexes )
        NULLIFY( Elm % FaceIndexes )
        NULLIFY( Elm % BubbleIndexes )
        Nodes % x => Element % NodeU
        Nodes % y => Element % NodeV
        Nodes % z => Element % NodeW
        CALL StabParam( Elm, Nodes, Element % NumberOfNodes, &
                 Element % StabilizationMK )

        DEALLOCATE(Elm)
      END IF

   END FUNCTION GetElementType
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute convection diffusion equation stab. parameter  for each and every
!> element of the model by solving the largest eigenvalue of
!
!> Lu = \lambda Gu,
!
!> L = (\nablda^2 u,\nabla^ w), G = (\nabla u,\nabla w)
!------------------------------------------------------------------------------
   SUBROUTINE StabParam(Element,Nodes,n,mK,hK,UseLongEdge)
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t), POINTER :: Element
      INTEGER :: n
      TYPE(Nodes_t) :: Nodes
      REAL(KIND=dp) :: mK
      REAL(KIND=dp), OPTIONAL :: hK
      LOGICAL, OPTIONAL :: UseLongEdge
!------------------------------------------------------------------------------
      INTEGER :: info,p,q,i,j,t,dim
      REAL(KIND=dp) :: EIGR(n),EIGI(n),Beta(n),s,ddp(3),ddq(3),dNodalBasisdx(n,n,3)
      REAL(KIND=dp) :: u,v,w,L(n-1,n-1),G(n-1,n-1),Work(16*n)
      REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3),detJ

      LOGICAL :: stat
      TYPE(GaussIntegrationPoints_t) :: IntegStuff

      IF ( Element % TYPE % BasisFunctionDegree <= 1 ) THEN
         SELECT CASE( Element % TYPE % ElementCode ) 
           CASE( 202, 303, 404, 504, 605, 706  )
              mK = 1.0d0 / 3.0d0
           CASE( 808 )
              mK = 1.0d0 / 6.0d0
         END SELECT
         IF ( PRESENT( hK ) ) hK = ElementDiameter( Element, Nodes, UseLongEdge)
         RETURN
      END IF

      dNodalBasisdx = 0._dp
      DO p=1,n
        u = Element % TYPE % NodeU(p)
        v = Element % TYPE % NodeV(p)
        w = Element % TYPE % NodeW(p)
        stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
        dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
      END DO

      dim = CoordinateSystemDimension()
      IntegStuff = GaussPoints( Element )
      L = 0.0d0
      G = 0.0d0
      DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis, &
                dBasisdx )

        s = detJ * IntegStuff % s(t)

        DO p=2,n
          DO q=2,n
            ddp = 0.0d0
            ddq = 0.0d0
            DO i=1,dim
              G(p-1,q-1) = G(p-1,q-1) + s * dBasisdx(p,i) * dBasisdx(q,i)
              ddp(i) = ddp(i) + SUM( dNodalBasisdx(p,1:n,i) * dBasisdx(1:n,i) )
              ddq(i) = ddq(i) + SUM( dNodalBasisdx(q,1:n,i) * dBasisdx(1:n,i) )
            END DO
            L(p-1,q-1) = L(p-1,q-1) + s * SUM(ddp) * SUM(ddq)
          END DO
        END DO
      END DO

      IF ( ALL(ABS(L) < AEPS) ) THEN
        mK = 1.0d0 / 3.0d0
        IF ( PRESENT(hK) ) THEN
          hK = ElementDiameter( Element,Nodes,UseLongEdge)
        END IF
        RETURN
      END IF


      CALL DSYGV( 1,'N','U',n-1,L,n-1,G,n-1,EIGR,Work,12*n,info )
      mK = EIGR(n-1)

      IF ( mK < 10*AEPS ) THEN
        mK = 1.0d0 / 3.0d0
        IF ( PRESENT(hK) ) THEN
          hK = ElementDiameter( Element,Nodes,UseLongEdge )
        END IF
        RETURN
      END IF

      IF ( PRESENT( hK ) ) THEN
        hK = SQRT( 2.0d0 / (mK * Element % TYPE % StabilizationMK) )
        mK = MIN( 1.0d0 / 3.0d0, Element % TYPE % StabilizationMK )
      ELSE
        SELECT CASE(Element % TYPE % ElementCode / 100)
        CASE(2,4,8) 
          mK = 4 * mK
        END SELECT
        mK = MIN( 1.0d0/3.0d0, 2/mK )
      END IF

!------------------------------------------------------------------------------
   END SUBROUTINE StabParam
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>   Given element structure return value of a quantity x given at element nodes
!>   at local coordinate point u inside the element. Element basis functions are
!>   used to compute the value. This is for 1D elements, and shouldnt propably
!>   be called directly by the user but trough the wrapper routine
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
             s = s + Coeff(i) * u**p(i)
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
          s = s + Coeff(i) * u**p(i)
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
!>   Given element structure return value of a quantity x given at element nodes
!>   at local coordinate point (u,vb) inside the element. Element basis functions
!>   are used to compute the value.This is for 2D elements, and shouldnt propably
!>   be called directly by the user but trough the wrapper routine
!>   InterpolateInElement.
!------------------------------------------------------------------------------
   FUNCTION InterpolateInElement2D( element,x,u,v ) RESULT(y)
!------------------------------------------------------------------------------
     TYPE(Element_t) :: element          !< element structure
     REAL(KIND=dp) :: u                  !< Point at which to evaluate the partial derivative
     REAL(KIND=dp) :: v                  !< Point at which to evaluate the partial derivative
     REAL(KIND=dp), DIMENSION(:) :: x    !< Nodal values of the quantity whose partial derivative we want to know
     REAL(KIND=dp) :: y                  !< value of the quantity y = x(u,v)
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

     elt => element % TYPE
     BasisFunctions => elt % BasisFunctions

     DO n=1,Elt % NumberOfNodes
       p => BasisFunctions(n) % p
       q => BasisFunctions(n) % q
       Coeff => BasisFunctions(n) % Coeff

       s = 0.0d0
       DO i=1,BasisFunctions(n) % n
          s = s + Coeff(i)*u**p(i)*v**q(i)
       END DO
       y(n) = s
     END DO
   END SUBROUTINE NodalBasisFunctions2D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>   Given element structure return value of the first partial derivative with
!>   respect to local coordinate u of i quantity x given at element nodes at local
!>   coordinate point u,v inside the element. Element basis functions are used to
!>   compute the value. 
!------------------------------------------------------------------------------
   FUNCTION FirstDerivativeInU2D( element,x,u,v ) RESULT(y)
!------------------------------------------------------------------------------
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: x(:)
!     INPUT: Nodal values of the quantity whose partial derivative we want to know
!
!    REAL(KIND=dp) :: u,v
!     INPUT: Point at which to evaluate the partial derivative
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: y
!      value of the quantity y = @x(u,v)/@u
!    
!******************************************************************************
   !
   ! Return first partial derivative in u of a quantity x at point u,v
   !
   !
   !

      TYPE(Element_t) :: element

      REAL(KIND=dp) :: u,v
      REAL(KIND=dp), DIMENSION(:) :: x

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

      REAL(KIND=dp) :: y,s,t

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
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: x(:)
!     INPUT: Nodal values of the quantity whose partial derivative we want to know
!
!    REAL(KIND=dp) :: u,v
!     INPUT: Point at which to evaluate the partial derivative
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: y
!      value of the quantity y = @x(u,v)/@v
!    
!------------------------------------------------------------------------------
    !
    ! Return first partial derivative in v of a quantity x at point u,v
    !
    !
    !
      TYPE(Element_t) :: element

      REAL(KIND=dp), DIMENSION(:) :: x
      REAL(KIND=dp) :: u,v

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: y,s,t

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
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: 
!
!    REAL(KIND=dp) :: u,v
!     INPUT: Point at which to evaluate the partial derivative
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: y
!      value of the quantity y = @x(u,v)/@u
!    
!------------------------------------------------------------------------------
   !
   ! Return first partial derivative in u of a quantity x at point u,v
   !
   !
   !

      TYPE(Element_t) :: element
      REAL(KIND=dp) :: u,v,y(:,:)

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

      REAL(KIND=dp) :: s,t

      TYPE(ElementType_t),POINTER :: elt
      REAL(KIND=dp), POINTER :: Coeff(:)
      INTEGER, POINTER :: p(:),q(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)

      INTEGER :: i,n

      elt => element % TYPE
      BasisFunctions => elt % BasisFunctions

      DO n = 1,elt % NumberOfNodes
        p => BasisFunctions(n) % p
        q => BasisFunctions(n) % q
        Coeff => BasisFunctions(n) % Coeff

        s = 0.0d0
        t = 0.0d0
        DO i = 1,BasisFunctions(n) % n
          IF (p(i)>=1) s = s + p(i)*Coeff(i)*u**(p(i)-1)*v**q(i)
          IF (q(i)>=1) t = t + q(i)*Coeff(i)*u**p(i)*v**(q(i)-1)
        END DO
        y(n,1) = s
        y(n,2) = t
      END DO

   END SUBROUTINE NodalFirstDerivatives2D
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!>   Given element structure return value of the second partial derivatives with
!>   respect to local coordinates of a quantity x given at element nodes at local
!>   coordinate point u,v inside the element. Element basis functions are used to
!>   compute the value. 
!------------------------------------------------------------------------------
   FUNCTION SecondDerivatives2D( element,x,u,v ) RESULT(ddx)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: x(:)
!     INPUT: Nodal values of the quantity whose partial derivatives we want to know
!
!    REAL(KIND=dp) :: u,v
!     INPUT: Point at which to evaluate the partial derivative
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: s
!      value of the quantity s = @^2x(u,v)/@v^2
!    
!------------------------------------------------------------------------------

      TYPE(Element_t) :: element

      REAL(KIND=dp), DIMENSION(:) :: x
      REAL(KIND=dp) :: u,v

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      TYPE(ElementType_t),POINTER :: elt
      REAL(KIND=dp), DIMENSION (2,2) :: ddx
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
!>   are used to compute the value. This is for 3D elements, and shouldnt propably
!>   be called directly by the user but trough the wrapper routine
!>   InterpolateInElement.
!------------------------------------------------------------------------------
   FUNCTION InterpolateInElement3D( element,x,u,v,w ) RESULT(y)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: x(:)
!     INPUT: Nodal values of the quantity whose value we want to know
!
!    REAL(KIND=dp) :: u,v,w
!     INPUT: Point at which to evaluate the value
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: y
!      value of the quantity y = x(u,v,w)
!    
!------------------------------------------------------------------------------
   !
   ! Return value of a quantity x at point u,v,w
   !
      TYPE(Element_t) :: element

      REAL(KIND=dp) :: u,v,w
      REAL(KIND=dp), DIMENSION(:) :: x
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: y

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
        y = y + x(1) * ( (1-u) * (1-v) - w + u*v*w * s ) / 4
        y = y + x(2) * ( (1+u) * (1-v) - w - u*v*w * s ) / 4
        y = y + x(3) * ( (1+u) * (1+v) - w + u*v*w * s ) / 4
        y = y + x(4) * ( (1-u) * (1+v) - w - u*v*w * s ) / 4
        y = y + x(5) * w
        RETURN
      ELSE IF ( Elt % ElementCode == 613 ) THEN
        IF ( w == 1 ) w = 1.0d0-1.0d-12
        s = 1.0d0 / (1-w)

        y = 0.0d0
        y = y + x(1)  * (-u-v-1) * ( (1-u) * (1-v) - w + u*v*w * s ) / 4
        y = y + x(2)  * ( u-v-1) * ( (1+u) * (1-v) - w - u*v*w * s ) / 4
        y = y + x(3)  * ( u+v-1) * ( (1+u) * (1+v) - w + u*v*w * s ) / 4
        y = y + x(4)  * (-u+v-1) * ( (1-u) * (1+v) - w - u*v*w * s ) / 4
        y = y + x(5)  * w*(2*w-1)
        y = y + x(6)  * (1+u-w)*(1-u-w)*(1-v-w) * s / 2
        y = y + x(7)  * (1+v-w)*(1-v-w)*(1+u-w) * s / 2
        y = y + x(8)  * (1+u-w)*(1-u-w)*(1+v-w) * s / 2
        y = y + x(9)  * (1+v-w)*(1-v-w)*(1-u-w) * s / 2
        y = y + x(10) * w * (1-u-w) * (1-v-w) * s
        y = y + x(11) * w * (1+u-w) * (1-v-w) * s
        y = y + x(12) * w * (1+u-w) * (1+v-w) * s
        y = y + x(13) * w * (1-u-w) * (1+v-w) * s
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
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: u
!     INPUT: Point at which to evaluate the value
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: y
!      value of the quantity y = x(u)
!    
!------------------------------------------------------------------------------

     TYPE(Element_t) :: element
     REAL(KIND=dp) :: u,v,w,y(:)

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

     elt => element % TYPE
     BasisFunctions => elt % BasisFunctions

     DO n=1,Elt % NumberOfNodes
       p => BasisFunctions(n) % p
       q => BasisFunctions(n) % q
       r => BasisFunctions(n) % r
       Coeff => BasisFunctions(n) % Coeff

       s = 0.0d0
       DO i=1,BasisFunctions(n) % n
          s = s + Coeff(i)*u**p(i)*v**q(i)*w**r(i)
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
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: x(:)
!     INPUT: Nodal values of the quantity whose partial derivative we want to know
!
!    REAL(KIND=dp) :: u,v,w
!     INPUT: Point at which to evaluate the partial derivative
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: y
!      value of the quantity y = @x(u,v,w)/@u
!    
!------------------------------------------------------------------------------
   !
   ! Return first partial derivative in u of a quantity x at point u,v,w
   !

      TYPE(Element_t) :: element

      REAL(KIND=dp) :: u,v,w
      REAL(KIND=dp), DIMENSION(:) :: x

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: y

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
  y = y + x(1) * ( -(1-v) + v*w * s ) / 4
  y = y + x(2) * (  (1-v) - v*w * s ) / 4
  y = y + x(3) * (  (1+v) + v*w * s ) / 4
  y = y + x(4) * ( -(1+v) - v*w * s ) / 4
  RETURN
ELSE IF ( Elt % ElementCode == 613 ) THEN
  IF ( w == 1 ) w = 1.0d0-1.0d-12
  s = 1.0d0 / (1-w)

  y = 0.0d0
  y = y + x(1)  * ( -( (1-u) * (1-v) - w + u*v*w * s ) + &
            (-u-v-1) * ( -(1-v) + v*w * s ) ) / 4

  y = y + x(2)  * (  ( (1+u) * (1-v) - w - u*v*w * s ) + &
            ( u-v-1) * (  (1-v) - v*w * s ) ) / 4

  y = y + x(3)  * (  ( (1+u) * (1+v) - w + u*v*w * s ) + &
            ( u+v-1) * (  (1+v) + v*w * s ) ) / 4

  y = y + x(4)  * ( -( (1-u) * (1+v) - w - u*v*w * s ) + &
            (-u+v-1) * ( -(1+v) - v*w * s ) ) / 4

  y = y + x(5)  * 0.0d0

  y = y + x(6)  * (  (1-u-w)*(1-v-w) - (1+u-w)*(1-v-w) ) * s / 2
  y = y + x(7)  * (  (1+v-w)*(1-v-w) ) * s / 2
  y = y + x(8)  * (  (1-u-w)*(1+v-w) - (1+u-w)*(1+v-w) ) * s / 2
  y = y + x(9)  * ( -(1+v-w)*(1-v-w) ) * s / 2

  y = y - x(10) * w * (1-v-w) * s
  y = y + x(11) * w * (1-v-w) * s
  y = y + x(12) * w * (1+v-w) * s
  y = y - x(13) * w * (1+v-w) * s

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
!
!  DESCRIPTION:
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: x(:)
!     INPUT: Nodal values of the quantity whose partial derivative we want to know
!
!    REAL(KIND=dp) :: u,v,w
!     INPUT: Point at which to evaluate the partial derivative
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: y
!      value of the quantity y = @x(u,v,w)/@v
!    
!------------------------------------------------------------------------------
   !
   ! Return first partial derivative in v of a quantity x at point u,v,w
   !

      TYPE(Element_t) :: element

      REAL(KIND=dp) :: u,v,w
      REAL(KIND=dp), DIMENSION(:) :: x

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: y

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
  y = y + x(1) * ( -(1-u) + u*w * s ) / 4
  y = y + x(2) * ( -(1+u) - u*w * s ) / 4
  y = y + x(3) * (  (1+u) + u*w * s ) / 4
  y = y + x(4) * (  (1-u) - u*w * s ) / 4

  RETURN
ELSE IF ( Elt % ElementCode == 613 ) THEN
  IF ( w == 1 ) w = 1.0d0-1.0d-12
  s = 1.0d0 / (1-w)

  y = 0.0d0
  y = y + x(1)  * ( -( (1-u) * (1-v) - w + u*v*w * s ) +  &
           (-u-v-1) * ( -(1-u) + u*w * s ) ) / 4

  y = y + x(2)  * ( -( (1+u) * (1-v) - w - u*v*w * s ) + &
           ( u-v-1) * ( -(1+u) - u*w * s ) ) / 4

  y = y + x(3)  * (  ( (1+u) * (1+v) - w + u*v*w * s ) + &
           ( u+v-1) * (  (1+u) + u*w * s ) ) / 4

  y = y + x(4)  * (  ( (1-u) * (1+v) - w - u*v*w * s ) + &
           (-u+v-1) * (  (1-u) - u*w * s ) ) / 4

  y = y + x(5)  * 0.0d0

  y = y - x(6)  *  (1+u-w)*(1-u-w) * s / 2
  y = y + x(7)  * ( (1-v-w)*(1+u-w) - (1+v-w)*(1+u-w) ) * s / 2
  y = y + x(8)  *  (1+u-w)*(1-u-w) * s / 2
  y = y + x(9)  * ( (1-v-w)*(1-u-w) - (1+v-w)*(1-u-w) ) * s / 2

  y = y - x(10) *  w * (1-u-w) * s
  y = y - x(11) *  w * (1+u-w) * s
  y = y + x(12) *  w * (1+u-w) * s
  y = y + x(13) *  w * (1-u-w) * s
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
!
!  DESCRIPTION:
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: x(:)
!     INPUT: Nodal values of the quantity whose partial derivative we want to know
!
!    REAL(KIND=dp) :: u,v,w
!     INPUT: Point at which to evaluate the partial derivative
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: y
!      value of the quantity y = @x(u,v,w)/@w
!    
!------------------------------------------------------------------------------
   !
   ! Return first partial derivative in u of a quantity x at point u,v,w
   !
   !

      TYPE(Element_t) :: element

      REAL(KIND=dp) :: u,v,w
      REAL(KIND=dp), DIMENSION(:) :: x

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: y

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
  y = y + x(1) * ( -1 + u*v*s**2 ) / 4
  y = y + x(2) * ( -1 - u*v*s**2 ) / 4
  y = y + x(3) * ( -1 + u*v*s**2 ) / 4
  y = y + x(4) * ( -1 - u*v*s**2 ) / 4
  y = y + x(5)
  RETURN
ELSE IF ( Elt % ElementCode == 613 ) THEN
  IF ( w == 1 ) w = 1.0d0-1.0d-12
  s = 1.0d0 / (1-w)

  y = 0.0d0
  y = y + x(1)  * (-u-v-1) * ( -1 + u*v*s**2 ) / 4
  y = y + x(2)  * ( u-v-1) * ( -1 - u*v*s**2 ) / 4
  y = y + x(3)  * ( u+v-1) * ( -1 + u*v*s**2 ) / 4
  y = y + x(4)  * (-u+v-1) * ( -1 - u*v*s**2 ) / 4

  y = y + x(5)  * (4*w-1)

  y = y + x(6)  * ( ( -(1-u-w)*(1-v-w) - (1+u-w)*(1-v-w) - (1+u-w)*(1-u-w) ) * s + &
                    ( 1+u-w)*(1-u-w)*(1-v-w) * s**2 ) / 2

  y = y + x(7)  * ( ( -(1-v-w)*(1+u-w) - (1+v-w)*(1+u-w) - (1+v-w)*(1-v-w) ) * s + &
                    ( 1+v-w)*(1-v-w)*(1+u-w) * s**2 ) / 2

  y = y + x(8)  * ( ( -(1-u-w)*(1+v-w) - (1+u-w)*(1+v-w) - (1+u-w)*(1-u-w) ) * s + &
                    ( 1+u-w)*(1-u-w)*(1+v-w) * s**2 ) / 2

  y = y + x(9)  * ( ( -(1-v-w)*(1-u-w) - (1+v-w)*(1-u-w) - (1+v-w)*(1-v-w) ) * s + &
                    ( 1+v-w)*(1-v-w)*(1-u-w) * s**2 ) / 2
                    
  y = y + x(10) * ( ( (1-u-w) * (1-v-w) - w * (1-v-w) - w * (1-u-w) ) * s  + &
                   w * (1-u-w) * (1-v-w) * s**2 )

  y = y + x(11) * ( ( (1+u-w) * (1-v-w) - w * (1-v-w) - w * (1+u-w) ) * s  + &
                   w * (1+u-w) * (1-v-w) * s**2 )

  y = y + x(12) * ( ( (1+u-w) * (1+v-w) - w * (1+v-w) - w * (1+u-w) ) * s  + &
                   w * (1+u-w) * (1+v-w) * s**2 )

  y = y + x(13) * ( ( (1-u-w) * (1+v-w) - w * (1+v-w) - w * (1-u-w) ) * s  + &
                   w * (1-u-w) * (1+v-w) * s**2 )
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
   SUBROUTINE NodalFirstDerivatives3D( y,element,u,v,w )
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    REAL(KIND=dp) :: 
!
!    REAL(KIND=dp) :: u,v
!     INPUT: Point at which to evaluate the partial derivative
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: y
!      value of the quantity y = @x(u,v)/@u
!    
!------------------------------------------------------------------------------
   !
   ! Return first partial derivative in u of a quantity x at point u,v
   !

      TYPE(Element_t) :: element
      REAL(KIND=dp) :: u,v,w,y(:,:)

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

      REAL(KIND=dp) :: s,t,z

      TYPE(ElementType_t),POINTER :: elt
      REAL(KIND=dp), POINTER :: Coeff(:)
      INTEGER, POINTER :: p(:),q(:),r(:)
      TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)

      INTEGER :: i,n

      elt => element % TYPE
      BasisFunctions => elt % BasisFunctions

      DO n = 1,elt % NumberOfNodes
        p => BasisFunctions(n) % p
        q => BasisFunctions(n) % q
        r => BasisFunctions(n) % r
        Coeff => BasisFunctions(n) % Coeff

        s = 0.0d0
        t = 0.0d0
        z = 0.0d0
        DO i = 1,BasisFunctions(n) % n
          IF (p(i)>=1) s = s + p(i)*Coeff(i)*u**(p(i)-1)*v**q(i)*w**r(i)
          IF (q(i)>=1) t = t + q(i)*Coeff(i)*u**p(i)*v**(q(i)-1)*w**r(i)
          IF (r(i)>=1) z = z + r(i)*Coeff(i)*u**p(i)*v**q(i)*w**(r(i)-1)
        END DO
        y(n,1) = s
        y(n,2) = t
        y(n,3) = z
      END DO
   END SUBROUTINE NodalFirstDerivatives3D
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of the second partial derivatives with
!>   respect to local coordinates of i quantity x given at element nodes at local
!>   coordinate point u,v inside the element. Element basis functions are used to
!>   compute the value. 
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
!    REAL(KIND=dp) :: u,v
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

      REAL(KIND=dp) :: s
      INTEGER :: i,j,k,l,n,m

!------------------------------------------------------------------------------
      elt => element % TYPE
      k = elt % NumberOfNodes
      BasisFunctions => elt % BasisFunctions

      ddx = 0.0d0

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
   SUBROUTINE NodalBasisFunctions( n, Basis, element, u, v, w)
!------------------------------------------------------------------------------
     INTEGER :: n                 !< The number of (background) element nodes
     REAL(KIND=dp) :: Basis(:)    !< The values of reference element basis
     TYPE(Element_t) :: element   !< The element structure
     REAL(KIND=dp) :: u,v,w       !< The coordinates of the reference element point
!------------------------------------------------------------------------------
     INTEGER   :: i, q, dim
     REAL(KIND=dp) :: NodalBasis(n)

     dim = Element % TYPE % DIMENSION

     IF ( isActivePElement(Element) ) THEN
       SELECT CASE(dim)
       CASE(1)
         CALL NodalBasisFunctions1D( Basis, element, u )
       CASE(2)
         IF (isPTriangle(Element)) THEN
           DO q=1,n
             Basis(q) = TriangleNodalPBasis(q, u, v)
           END DO
         ELSE IF (isPQuad(Element)) THEN
           DO q=1,n
             Basis(q) = QuadNodalPBasis(q, u, v)
           END DO
         END IF
       CASE(3)
         IF (isPTetra( Element )) THEN
           DO q=1,n
             Basis(q) = TetraNodalPBasis(q, u, v, w)
           END DO
         ELSE IF (isPWedge( Element )) THEN
           DO q=1,n
             Basis(q) = WedgeNodalPBasis(q, u, v, w)
           END DO
         ELSE IF (isPPyramid( Element )) THEN
           DO q=1,n
             Basis(q) = PyramidNodalPBasis(q, u, v, w)
           END DO
         ELSE IF (isPBrick( Element )) THEN
           DO q=1,n
             Basis(q) = BrickNodalPBasis(q, u, v, w)
           END DO
         END IF
       END SELECT
     ELSE
       SELECT CASE( dim )
       CASE(1)
         CALL NodalBasisFunctions1D( Basis, element, u )
       CASE(2)
         CALL NodalBasisFunctions2D( Basis, element, u,v )
       CASE(3)
         IF ( Element % TYPE % ElementCode/100==6 ) THEN
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
   SUBROUTINE NodalFirstDerivatives( n, dLBasisdx, element, u, v, w)
!------------------------------------------------------------------------------
     INTEGER :: n                    !< The number of (background) element nodes
     REAL(KIND=dp) :: dLBasisdx(:,:) !< The gradient of reference element basis functions
     TYPE(Element_t) :: element      !< The element structure
     REAL(KIND=dp) :: u,v,w          !< The coordinates of the reference element point
!------------------------------------------------------------------------------
     INTEGER   :: i, q, dim
     REAL(KIND=dp) :: NodalBasis(n)
!------------------------------------------------------------------------------
     dim = Element % TYPE % DIMENSION

     IF ( IsActivePElement(Element) ) THEN
       SELECT CASE(dim)
       CASE(1)
         CALL NodalFirstDerivatives1D( dLBasisdx, element, u )
       CASE(2)
         IF (isPTriangle(Element)) THEN
           DO q=1,n
             dLBasisdx(q,1:2) = dTriangleNodalPBasis(q, u, v)
           END DO
         ELSE IF (isPQuad(Element)) THEN
           DO q=1,n
             dLBasisdx(q,1:2) = dQuadNodalPBasis(q, u, v)
           END DO
         END IF
       CASE(3)
         IF (isPTetra( Element )) THEN
           DO q=1,n
             dLBasisdx(q,1:3) = dTetraNodalPBasis(q, u, v, w)
           END DO
         ELSE IF (isPWedge( Element )) THEN
           DO q=1,n
             dLBasisdx(q,1:3) = dWedgeNodalPBasis(q, u, v, w)
           END DO
         ELSE IF (isPPyramid( Element )) THEN
           DO q=1,n
             dLBasisdx(q,1:3) = dPyramidNodalPBasis(q, u, v, w)
           END DO
         ELSE IF (isPBrick( Element )) THEN
           DO q=1,n
             dLBasisdx(q,1:3) = dBrickNodalPBasis(q, u, v, w)
           END DO
         END IF
       END SELECT
     ELSE
       SELECT CASE(dim)
       CASE(1)
         CALL NodalFirstDerivatives1D( dLBasisdx, element, u )
       CASE(2)
         CALL NodalFirstDerivatives2D( dLBasisdx, element, u,v )
       CASE(3)
         IF ( Element % TYPE % ElementCode / 100 == 6 ) THEN
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
   SUBROUTINE ElementBasisDegree( Element, BasisDegree )
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Element_t), TARGET :: Element             !< Element structure
     INTEGER :: BasisDegree(:)!< Degree of each basis function in Basis(:) vector. 
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: t,s
     LOGICAL :: invert, degrees
     INTEGER :: i, j, k, l, q, p, f, n, nb, dim, cdim, locali, localj,  &
          tmp(4), direction(4)

     TYPE(Element_t) :: Bubble
     TYPE(Element_t), POINTER :: Edge, Face
!------------------------------------------------------------------------------

     n    = Element % TYPE % NumberOfNodes
     dim  = Element % TYPE % DIMENSION
     cdim = CoordinateSystemDimension()

     BasisDegree = 0
     BasisDegree(1:n) = Element % Type % BasisFunctionDegree

     IF ( isActivePElement(element) ) THEN

       ! Check for need of P basis degrees and set degree of
       ! linear basis if vector asked:
       ! ---------------------------------------------------
       BasisDegree(1:n) = 1
       q = n

!------------------------------------------------------------------------------
     SELECT CASE( Element % TYPE % ElementCode ) 
!------------------------------------------------------------------------------

     ! P element code for line element:
     ! --------------------------------
     CASE(202)
        ! Bubbles of line element
        IF (Element % BDOFs > 0) THEN
           ! For each bubble in line element get value of basis function
           DO i=1, Element % BDOFs
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
        IF ( Element % BDOFs > 0 ) THEN
           ! Get element p
           p = Element % PDefs % P

           nb = MAX( GetBubbleDOFs( Element, p ), Element % BDOFs )
           p = CEILING( ( 3.0d0+SQRT(1.0d0+8.0d0*nb) ) / 2.0d0 )
           
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
        IF ( Element % BDOFs > 0 ) THEN
          ! Get element P
           p = Element % PDefs % P

           nb = MAX( GetBubbleDOFs( Element, p ), Element % BDOFs )
           p = CEILING( ( 5.0d0+SQRT(1.0d0+8.0d0*nb) ) / 2.0d0 )
          
           DO i=2,(p-2)
              DO j=2,(p-i)
                 IF ( q >= SIZE(BasisDegree) ) CYCLE
                 q = q + 1
                 BasisDegree(q) = i+j
              END DO
           END DO
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
        IF ( Element % BDOFs > 0 ) THEN
           p = Element % PDefs % P

           nb = MAX( GetBubbleDOFs(Element, p), Element % BDOFs )
           p=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+2)

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
                 DO i=2,p-2
                    DO j=2,p-i
                       IF ( q >= SIZE(BasisDegree) ) CYCLE
                       q = q + 1
                       BasisDegree(q) = i+j
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
        IF (Element % BDOFs >= 0) THEN 
           ! Get element p
           p = Element % PDefs % p
           nb = MAX( GetBubbleDOFs(Element, p), Element % BDOFs )
           p=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+2)

           ! Calculate value of bubble functions from indexes
           ! i,j,k=0,..,p-4 i+j+k=0,..,p-4
           DO i=0,p-4
              DO j=0,p-i-4
                 DO k=0,p-i-j-4
                    IF ( q >= SIZE(BasisDegree)) CYCLE
                    q = q + 1
                    BasisDegree(q) = 4+i+j+k
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

                 ! Use basis compatible with pyramid if necessary
                 ! @todo Correct this!
                 IF (Edge % PDefs % pyramidQuadEdge) THEN
                    CALL Fatal('ElementInfo','Pyramid compatible wedge edge basis NIY!')
                 END IF
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
                 DO i=2,p-2
                    DO j=2,p-i
                       IF ( q >= SIZE(BasisDegree) ) CYCLE
                       q = q + 1
                       BasisDegree(q) = i+j
                    END DO
                 END DO
              END SELECT
                           
           END DO
        END IF

        ! Bubbles of P Wedge
        IF ( Element % BDOFs > 0 ) THEN
           ! Get p from element
           p = Element % PDefs % P
           nb = MAX( GetBubbleDOFs( Element, p ), Element % BDOFs )
           p=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+3)
           
           ! For each bubble calculate value of basis function and its derivative
           ! for index pairs i,j=0,..,p-5 k=2,..,p-3 i+j+k=2,..,p-3
           DO i=0,p-5
              DO j=0,p-5-i
                 DO k=2,p-3-i-j
                    IF ( q >= SIZE(BasisDegree) ) CYCLE
                    q = q + 1
                    BasisDegree(q) = 3+i+j+k
                 END DO
              END DO
           END DO
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
              DO i=2,p-2
                 DO j=2,p-i
                    IF ( q >= SIZE(BasisDegree) ) CYCLE
                    q = q + 1
                    BasisDegree(q) = i+j
                 END DO
              END DO
           END DO
        END IF

        ! Bubbles of p brick
        IF ( Element % BDOFs > 0 ) THEN
           ! Get p from bubble DOFs 
           p = Element % PDefs % P
           nb = MAX( GetBubbleDOFs(Element, p), Element % BDOFs )
           p=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+4)

           ! For each bubble calculate value of basis function and its derivative
           ! for index pairs i,j,k=2,..,p-4, i+j+k=6,..,p
           DO i=2,p-4
              DO j=2,p-i-2
                 DO k=2,p-i-j
                    IF ( q >= SIZE(BasisDegree) ) CYCLE
                    q = q + 1
                    BasisDegree(q) = i+j+k
                 END DO
              END DO
           END DO
        END IF

     END SELECT
     END IF ! P element flag check
!------------------------------------------------------------------------------
   END SUBROUTINE ElementBasisDegree
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Return the referencial description b(f(p)) of the basis function b(x),
!>  with f mapping points p on a reference element to points x on a physical
!>  element. The referencial description of the spatial gradient field grad b 
!>  and, if requested, the second spatial derivatives may also be returned.
!>  Also return the square root of the determinant of the metric tensor
!>  (=sqrt(det(J^TJ))) related to the mapping f.
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ElementInfo( Element, Nodes, u, v, w, detJ, &
     Basis, dBasisdx, ddBasisddx, SecondDerivatives, Bubbles, BasisDegree, EdgeBasis, RotBasis ) RESULT(stat)
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Element_t), TARGET :: Element             !< Element structure
     TYPE(Nodes_t)   :: Nodes                       !< Element nodal coordinates.
     REAL(KIND=dp) :: u                             !< 1st local coordinate at which to calculate the basis function.
     REAL(KIND=dp) :: v                             !< 2nd local coordinate.
     REAL(KIND=dp) :: w                             !< 3rd local coordinate.
     REAL(KIND=dp) :: detJ                          !< Square root of determinant of element coordinate system metric
     REAL(KIND=dp) :: Basis(:)                      !< Basis function values at p=(u,v,w)
     REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)       !< Global first derivatives of basis functions at (u,v,w)
     REAL(KIND=dp), OPTIONAL :: ddBasisddx(:,:,:)   !< Global second derivatives of basis functions at (u,v,w) if requested
     INTEGER, OPTIONAL :: BasisDegree(:)            !< Degree of each basis function in Basis(:) vector. 
	                                                !! May be used with P element basis functions
     LOGICAL, OPTIONAL :: SecondDerivatives         !< Are the second derivatives needed? (still present for historical reasons)
     LOGICAL, OPTIONAL :: Bubbles                   !< Are the bubbles to be avaluated.
     REAL(KIND=dp), OPTIONAL :: EdgeBasis(:,:)      !< If present, the values of H(curl)-conforming basis functions B(f(p))
     REAL(KIND=dp), OPTIONAL :: RotBasis(:,:)       !< The referencial description of the spatial curl of B
     LOGICAL :: Stat                                !< If .FALSE. element is degenerate.
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: BubbleValue, dBubbledx(3), t, s, LtoGMap(3,3)
     LOGICAL :: invert, degrees
     INTEGER :: i, j, k, l, q, p, f, n, nb, dim, cdim, locali, localj,  &
          tmp(4), direction(4)
     REAL(KIND=dp) :: LinBasis(8), dLinBasisdx(8,3), ElmMetric(3,3)

     REAL(KIND=dp) :: NodalBasis(Element % TYPE % NumberOfNodes), &
             dLBasisdx(MAX(SIZE(Nodes % x),SIZE(Basis)),3)

     TYPE(Element_t) :: Bubble
     TYPE(Element_t), POINTER :: Edge, Face
!------------------------------------------------------------------------------
     IF(PRESENT(EdgeBasis)) THEN
       stat = EdgeElementInfo(Element,Nodes,u,v,w,detF=Detj,Basis=Basis, &
            EdgeBasis=EdgeBasis,RotBasis=RotBasis,dBasisdx=dBasisdx,ApplyPiolaTransform=.TRUE.)
       RETURN
     END IF

     stat = .TRUE.
     n    = Element % TYPE % NumberOfNodes
     dim  = Element % TYPE % DIMENSION
     cdim = CoordinateSystemDimension()

     IF ( Element % TYPE % ElementCode == 101 ) THEN
        detJ = 1.0d0
        Basis(1) = 1.0d0
        IF ( PRESENT(dBasisdx) ) dBasisdx(1,:) = 0.0d0
        RETURN
     END IF

     Basis = 0.0d0
     CALL NodalBasisFunctions(n, Basis, element, u, v, w)

     dLbasisdx = 0.0d0
     CALL NodalFirstDerivatives(n, dLBasisdx, element, u, v, w)

     q = n

     ! P ELEMENT CODE:
     ! ---------------
     IF ( isActivePElement(element) ) THEN

      ! Check for need of P basis degrees and set degree of
      ! linear basis if vector asked:
      ! ---------------------------------------------------
      degrees = .FALSE.
      IF ( PRESENT(BasisDegree)) THEN 
        degrees = .TRUE.
        BasisDegree = 0
        BasisDegree(1:n) = 1
      END IF

!------------------------------------------------------------------------------
     SELECT CASE( Element % TYPE % ElementCode ) 
!------------------------------------------------------------------------------

     ! P element code for line element:
     ! --------------------------------
     CASE(202)
        ! Bubbles of line element
        IF (Element % BDOFs > 0) THEN
           ! For boundary element integration check direction
           invert = .FALSE.
           IF ( Element % PDefs % isEdge .AND. &
                Element % NodeIndexes(1)>Element % NodeIndexes(2) ) invert = .TRUE.

           ! For each bubble in line element get value of basis function
           DO i=1, Element % BDOFs
              IF (q >= SIZE(Basis)) CYCLE
              q = q + 1
              
              Basis(q) = LineBubblePBasis(i+1,u,invert)
              dLBasisdx(q,1) = dLineBubblePBasis(i+1,u,invert)
              
              ! Polynomial degree of basis function to vector
              IF (degrees) BasisDegree(q) = 1+i
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

              ! Get local number of edge start and endpoint nodes
              tmp(1:2) = getTriangleEdgeMap(i)
              locali = tmp(1)
              localj = tmp(2)

              ! Invert edge for parity if needed
              invert = .FALSE.
              IF ( Element % NodeIndexes(locali)>Element % NodeIndexes(localj) ) invert=.TRUE.

              ! For each dof in edge get value of p basis function 
              DO k=1,Edge % BDOFs
                 IF (q >= SIZE(Basis)) CYCLE
                 q = q + 1
                 
                 ! Value of basis functions for edge=i and i=k+1 by parity
                 Basis(q) = TriangleEdgePBasis(i, k+1, u, v, invert)
                 ! Value of derivative of basis function
                 dLBasisdx(q,1:2) = dTriangleEdgePBasis(i, k+1, u, v, invert)
                 
                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
              END DO
           END DO 
        END IF

        ! Bubbles of p triangle      
        IF ( Element % BDOFs > 0 ) THEN
           ! Get element p
           p = Element % PDefs % P

           nb = MAX( GetBubbleDOFs( Element, p ), Element % BDOFs )
           p = CEILING( ( 3.0d0+SQRT(1.0d0+8.0d0*nb) ) / 2.0d0 )
           
           ! For boundary element direction needs to be calculated
           IF (Element % PDefs % isEdge) THEN
              direction = 0
              ! Get direction of this face (mask for face = boundary element nodes)
              direction(1:3) = getTriangleFaceDirection(Element, [ 1,2,3 ])
           END IF

           DO i = 0,p-3
              DO j = 0,p-i-3
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! Get bubble basis functions and their derivatives
                 ! 3d Boundary element has a direction
                 IF (Element % PDefs % isEdge) THEN
                    Basis(q) = TriangleEBubblePBasis(i,j,u,v,direction) 
                    dLBasisdx(q,1:2) = dTriangleEBubblePBasis(i,j,u,v,direction)
                 ELSE
                 ! 2d element bubbles have no direction
                    Basis(q) = TriangleBubblePBasis(i,j,u,v) 
                    dLBasisdx(q,1:2) = dTriangleBubblePBasis(i,j,u,v)
                 END IF
                 
                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 3+i+j
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

              ! Choose correct parity by global edge dofs
              tmp(1:2) = getQuadEdgeMap(i)
              locali = tmp(1)
              localj = tmp(2)
              
              ! Invert parity if needed
              invert = .FALSE.
              IF (Element % NodeIndexes(locali) > Element % NodeIndexes(localj)) invert = .TRUE. 

              ! For each DOF in edge calculate value of p basis function
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! For pyramid square face edges use different basis
                 IF (Edge % PDefs % pyramidQuadEdge) THEN
                    Basis(q) = QuadPyraEdgePBasis(i,k+1,u,v,invert)
                    dLBasisdx(q,1:2) = dQuadPyraEdgePBasis(i,k+1,u,v,invert)
                 ! Normal case, use basis of quadrilateral
                 ELSE
                    ! Get values of basis functions for edge=i and i=k+1 by parity
                    Basis(q) = QuadEdgePBasis(i,k+1,u,v,invert)
                    ! Get value of derivatives of basis functions
                    dLBasisdx(q,1:2) = dQuadEdgePBasis(i,k+1,u,v,invert)
                 END IF
                 
                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
              END DO              
           END DO         
        END IF

        ! Bubbles of p quadrilateral
        IF ( Element % BDOFs > 0 ) THEN
          ! Get element P
           p = Element % PDefs % P

           nb = MAX( GetBubbleDOFs( Element, p ), Element % BDOFs )
           p = CEILING( ( 5.0d0+SQRT(1.0d0+8.0d0*nb) ) / 2.0d0 )

           ! For boundary element direction needs to be calculated
           IF (Element % PDefs % isEdge) THEN
              direction = 0
              direction = getSquareFaceDirection(Element, [ 1,2,3,4 ])
           END IF
          
           ! For each bubble calculate value of p basis function
           ! and their derivatives for index pairs i,j>=2, i+j=4,...,p
           DO i=2,(p-2)
              DO j=2,(p-i)
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1
                 
                 ! Get values of bubble functions
                 ! 3D boundary elements have a direction
                 IF (Element % PDefs % isEdge) THEN
                    Basis(q) = QuadBubblePBasis(i,j,u,v,direction)
                    dLBasisdx(q,1:2) = dQuadBubblePBasis(i,j,u,v,direction)
                 ELSE
                 ! 2d element bubbles have no direction
                    Basis(q) = QuadBubblePBasis(i,j,u,v)
                    dLBasisdx(q,1:2) = dQuadBubblePBasis(i,j,u,v)
                 END IF

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = i+j
              END DO
           END DO
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
                 IF (q >= SIZE(Basis)) CYCLE
                 q = q + 1

                 Basis(q) = TetraEdgePBasis(i,k+1,u,v,w, Element % PDefs % TetraType)
                 dLBasisdx(q,1:3) = dTetraEdgePBasis(i,k+1,u,v,w, Element % PDefs % TetraType)

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
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
                    IF (q >= SIZE(Basis)) CYCLE
                    q = q + 1 
                    
                    Basis(q) = TetraFacePBasis(F,i,j,u,v,w, Element % PDefs % TetraType)
                    dLBasisdx(q,1:3) = dTetraFacePBasis(F,i,j,u,v,w, Element % PDefs % TetraType)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = 3+i+j
                 END DO
              END DO
           END DO
        END IF

        ! Bubbles of p tetrahedron
        IF ( Element % BDOFs > 0 ) THEN
           p = Element % PDefs % P

           nb = MAX( GetBubbleDOFs(Element, p), Element % BDOFs )
           p=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+2)

           ! For each DOF in bubbles calculate value of bubble functions
           ! and their derivatives for index pairs
           ! i,j,k=0,..,p-4 i+j+k=0,..,p-4
           DO i=0,p-4
              DO j=0,p-i-4
                 DO k=0,p-i-j-4
                    IF (q >= SIZE(Basis)) CYCLE
                    q = q + 1

                    Basis(q) = TetraBubblePBasis(i,j,k,u,v,w)
                    dLBasisdx(q,1:3) = dTetraBubblePBasis(i,j,k,u,v,w)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = 4+i+j+k
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
              
              ! Get local indexes of current edge
              tmp(1:2) = getPyramidEdgeMap(i)
              locali = tmp(1)
              localj = tmp(2)

              ! Determine edge direction
              invert = .FALSE.
              
              ! Invert edge if local first node has greater global index than second one
              IF ( Element % NodeIndexes(locali) > Element % NodeIndexes(localj) ) invert = .TRUE.

              ! For each DOF in edge calculate values of edge functions
              ! and their derivatives for edge=i and i=k+1
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! Get values of edge basis functions and their derivatives
                 Basis(q) = PyramidEdgePBasis(i,k+1,u,v,w,invert)
                 dLBasisdx(q,1:3) = dPyramidEdgePBasis(i,k+1,u,v,w,invert)

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
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
                 direction = 0
                 ! Get global direction vector for enforcing parity
                 tmp(1:4) = getPyramidFaceMap(F)
                 direction(1:4) = getSquareFaceDirection( Element, tmp(1:4) )
                 
                 ! For each face calculate values of functions from index
                 ! pairs i,j=2,..,p-2 i+j=4,..,p
                 DO i=2,p-2
                    DO j=2,p-i
                       IF ( q >= SIZE(Basis) ) CYCLE
                       q = q + 1
                       
                       Basis(q) = PyramidFacePBasis(F,i,j,u,v,w,direction)
                       dLBasisdx(q,:) = dPyramidFacePBasis(F,i,j,u,v,w,direction)
                       
                       ! Polynomial degree of basis function to vector
                       IF (degrees) BasisDegree(q) = i+j
                    END DO
                 END DO

              CASE (2,3,4,5)
                 direction = 0
                 ! Get global direction vector for enforcing parity
                 tmp(1:4) = getPyramidFaceMap(F) 
                 direction(1:3) = getTriangleFaceDirection( Element, tmp(1:3) )
                 
                 ! For each face calculate values of functions from index
                 ! pairs i,j=0,..,p-3 i+j=0,..,p-3
                 DO i=0,p-3
                    DO j=0,p-i-3
                       IF ( q >= SIZE(Basis) ) CYCLE
                       q = q + 1

                       Basis(q) = PyramidFacePBasis(F,i,j,u,v,w,direction)
                       dLBasisdx(q,:) = dPyramidFacePBasis(F,i,j,u,v,w,direction)

                       ! Polynomial degree of basis function to vector
                       IF (degrees) BasisDegree(q) = 3+i+j
                    END DO
                 END DO
              END SELECT    
           END DO
        END IF

        ! Bubbles of P Pyramid
        IF (Element % BDOFs >= 0) THEN 
           ! Get element p
           p = Element % PDefs % p
           nb = MAX( GetBubbleDOFs(Element, p), Element % BDOFs )
           p=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+2)

           ! Calculate value of bubble functions from indexes
           ! i,j,k=0,..,p-4 i+j+k=0,..,p-4
           DO i=0,p-4
              DO j=0,p-i-4
                 DO k=0,p-i-j-4
                    IF ( q >= SIZE(Basis)) CYCLE
                    q = q + 1

                    Basis(q) = PyramidBubblePBasis(i,j,k,u,v,w)
                    dLBasisdx(q,:) = dPyramidBubblePBasis(i,j,k,u,v,w)
                    
                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = 4+i+j+k
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
              
              ! Get local indexes of current edge
              tmp(1:2) = getWedgeEdgeMap(i)
              locali = tmp(1)
              localj = tmp(2)

              ! Determine edge direction
              invert = .FALSE.
              ! Invert edge if local first node has greater global index than second one
              IF ( Element % NodeIndexes(locali) > Element % NodeIndexes(localj) ) invert = .TRUE.
       
              ! For each DOF in edge calculate values of edge functions
              ! and their derivatives for edge=i and i=k+1
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! Use basis compatible with pyramid if necessary
                 ! @todo Correct this!
                 IF (Edge % PDefs % pyramidQuadEdge) THEN
                    CALL Fatal('ElementInfo','Pyramid compatible wedge edge basis NIY!')
                 END IF

                 ! Get values of edge basis functions and their derivatives
                 Basis(q) = WedgeEdgePBasis(i,k+1,u,v,w,invert)
                 dLBasisdx(q,1:3) = dWedgeEdgePBasis(i,k+1,u,v,w,invert)

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
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
                 direction = 0
                 ! Get global direction vector for enforcing parity
                 tmp(1:4) = getWedgeFaceMap(F) 
                 direction(1:3) = getTriangleFaceDirection( Element, tmp(1:3) )
                 
                 ! For each face calculate values of functions from index
                 ! pairs i,j=0,..,p-3 i+j=0,..,p-3
                 DO i=0,p-3
                    DO j=0,p-i-3
                       IF ( q >= SIZE(Basis) ) CYCLE
                       q = q + 1

                       Basis(q) = WedgeFacePBasis(F,i,j,u,v,w,direction)
                       dLBasisdx(q,:) = dWedgeFacePBasis(F,i,j,u,v,w,direction)

                       ! Polynomial degree of basis function to vector
                       IF (degrees) BasisDegree(q) = 3+i+j
                    END DO
                 END DO
              CASE (3,4,5)
                 direction = 0
                 ! Get global direction vector for enforcing parity
                 invert = .FALSE.
                 tmp(1:4) = getWedgeFaceMap(F)
                 direction(1:4) = getSquareFaceDirection( Element, tmp(1:4) )
                 
                 ! First and second node must form a face in upper or lower triangle
                 IF (.NOT. wedgeOrdering(direction)) THEN
                    invert = .TRUE.
                    tmp(1) = direction(2)
                    direction(2) = direction(4)
                    direction(4) = tmp(1)
                 END IF

                 ! For each face calculate values of functions from index
                 ! pairs i,j=2,..,p-2 i+j=4,..,p
                 DO i=2,p-2
                    DO j=2,p-i
                       IF ( q >= SIZE(Basis) ) CYCLE
                       q = q + 1

                       IF (.NOT. invert) THEN
                          Basis(q) = WedgeFacePBasis(F,i,j,u,v,w,direction)
                          dLBasisdx(q,:) = dWedgeFacePBasis(F,i,j,u,v,w,direction)
                       ELSE
                          Basis(q) = WedgeFacePBasis(F,j,i,u,v,w,direction)
                          dLBasisdx(q,:) = dWedgeFacePBasis(F,j,i,u,v,w,direction)
                       END IF

                       ! Polynomial degree of basis function to vector
                       IF (degrees) BasisDegree(q) = i+j
                    END DO
                 END DO
              END SELECT
                           
           END DO
        END IF

        ! Bubbles of P Wedge
        IF ( Element % BDOFs > 0 ) THEN
           ! Get p from element
           p = Element % PDefs % P
           nb = MAX( GetBubbleDOFs( Element, p ), Element % BDOFs )
           p=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+3)
           
           ! For each bubble calculate value of basis function and its derivative
           ! for index pairs i,j=0,..,p-5 k=2,..,p-3 i+j+k=2,..,p-3
           DO i=0,p-5
              DO j=0,p-5-i
                 DO k=2,p-3-i-j
                    IF ( q >= SIZE(Basis) ) CYCLE
                    q = q + 1

                    Basis(q) = WedgeBubblePBasis(i,j,k,u,v,w)
                    dLBasisdx(q,:) = dWedgeBubblePBasis(i,j,k,u,v,w)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = 3+i+j+k
                 END DO
              END DO
           END DO
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
              
              ! Get local indexes of current edge
              tmp(1:2) = getBrickEdgeMap(i)
              locali = tmp(1)
              localj = tmp(2)
              
              ! Determine edge direction
              invert = .FALSE.
              
              ! Invert edge if local first node has greater global index than second one
              IF ( Element % NodeIndexes(locali) > Element % NodeIndexes(localj) ) invert = .TRUE.
              
              ! For each DOF in edge calculate values of edge functions
              ! and their derivatives for edge=i and i=k+1
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! For edges connected to pyramid square face, use different basis
                 IF (Edge % PDefs % pyramidQuadEdge) THEN
                    ! Get values of edge basis functions and their derivatives
                    Basis(q) = BrickPyraEdgePBasis(i,k+1,u,v,w,invert)
                    dLBasisdx(q,1:3) = dBrickPyraEdgePBasis(i,k+1,u,v,w,invert)
                 ! Normal case. Use standard brick edge functions
                 ELSE
                    ! Get values of edge basis functions and their derivatives
                    Basis(q) = BrickEdgePBasis(i,k+1,u,v,w,invert)
                    dLBasisdx(q,1:3) = dBrickEdgePBasis(i,k+1,u,v,w,invert)
                 END IF

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
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
              
              ! Generate direction vector for this face
              tmp(1:4) = getBrickFaceMap(F)
              direction(1:4) = getSquareFaceDirection(Element, tmp)
              
              ! For each face calculate values of functions from index
              ! pairs i,j=2,..,p-2 i+j=4,..,p
              DO i=2,p-2
                 DO j=2,p-i
                    IF ( q >= SIZE(Basis) ) CYCLE
                    q = q + 1
                    Basis(q) = BrickFacePBasis(F,i,j,u,v,w,direction)
                    dLBasisdx(q,:) = dBrickFacePBasis(F,i,j,u,v,w,direction)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = i+j
                 END DO
              END DO
           END DO
        END IF

        ! Bubbles of p brick
        IF ( Element % BDOFs > 0 ) THEN
           ! Get p from bubble DOFs 
           p = Element % PDefs % P
           nb = MAX( GetBubbleDOFs(Element, p), Element % BDOFs )
           p=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+4)

           
           ! For each bubble calculate value of basis function and its derivative
           ! for index pairs i,j,k=2,..,p-4, i+j+k=6,..,p
           DO i=2,p-4
              DO j=2,p-i-2
                 DO k=2,p-i-j
                    IF ( q >= SIZE(Basis) ) CYCLE
                    q = q + 1
                    Basis(q) = BrickBubblePBasis(i,j,k,u,v,w)
                    dLBasisdx(q,:) = dBrickBubblePBasis(i,j,k,u,v,w)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = i+j+k
                 END DO
              END DO
           END DO
        END IF

     END SELECT
     END IF ! P element flag check
!------------------------------------------------------------------------------

     ! Element (contravariant) metric and square root of determinant
     !--------------------------------------------------------------
     IF ( .NOT. ElementMetric( q, Element, Nodes, &
           ElmMetric, detJ, dLBasisdx, LtoGMap ) ) THEN
        stat = .FALSE.
        RETURN
     END IF

     ! Get global first derivatives:
     !------------------------------
     IF ( PRESENT(dBasisdx) ) THEN
       dBasisdx = 0.0d0
       DO i=1,q
         DO j=1,cdim
            DO k=1,dim
              dBasisdx(i,j) = dBasisdx(i,j) + dLBasisdx(i,k)*LtoGMap(j,k)
            END DO
         END DO
       END DO
     END IF

     ! Get matrix of second derivatives, if needed:
     !---------------------------------------------
     IF ( PRESENT(ddBasisddx) .AND. PRESENT(SecondDerivatives) ) THEN
       IF ( SecondDerivatives ) THEN
         NodalBasis = 0.0d0
         ddBasisddx(1:n,:,:) = 0.0d0
         DO q=1,n
           NodalBasis(q) = 1.0d0
           CALL GlobalSecondDerivatives(Element,Nodes,NodalBasis, &
               ddBasisddx(q,:,:),u,v,w,ElmMetric,dLBasisdx )
           NodalBasis(q) = 0.0d0
         END DO
       END IF
     END IF

!------------------------------------------------------------------------------
!    Generate bubble basis functions, if requested. Bubble basis is as follows:
!    B_i (=(N_(i+n)) = B * N_i, where N_i:s are the nodal basis functions of
!    the element, and B the basic bubble, i.e. the product of nodal basis
!    functions of the corresponding linear element for triangles and tetras,
!    and product of two diagonally opposed nodal basisfunctions of the
!    correspoding (bi-,tri-)linear element for 1d-elements, quads and hexas.
!------------------------------------------------------------------------------
     IF ( PRESENT( Bubbles ) ) THEN
       Bubble % BDOFs = 0
       NULLIFY( Bubble % PDefs )
       NULLIFY( Bubble % EdgeIndexes )
       NULLIFY( Bubble % FaceIndexes )
       NULLIFY( Bubble % BubbleIndexes )

       IF ( Bubbles .AND. SIZE(Basis) >= 2*n ) THEN

         SELECT CASE(Element % TYPE % ElementCode / 100)
           CASE(2)

              IF ( Element % TYPE % ElementCode == 202 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(202)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                          LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(2)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(1,j) * LinBasis(2)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(2,j) * LinBasis(1)
                END DO
              END DO

           CASE(3)

              IF ( Element % TYPE % ElementCode == 303 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(303)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                            LinBasis, dLinBasisdx )
              END IF
  
              BubbleValue = LinBasis(1) * LinBasis(2) * LinBasis(3)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(1,j) * LinBasis(2) * LinBasis(3)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(2,j) * LinBasis(1) * LinBasis(3)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(3,j) * LinBasis(1) * LinBasis(2)
                END DO
              END DO

           CASE(4)

              IF ( Element % TYPE % ElementCode == 404 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(404)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                             LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(3)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                         dLinBasisdx(1,j) * LinBasis(3)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                         dLinBasisdx(3,j) * LinBasis(1)
                END DO
              END DO

           CASE(5)

              IF ( Element % TYPE % ElementCode == 504 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(504)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                            LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(2) * LinBasis(3) * LinBasis(4)
              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(1,j) * &
                                    LinBasis(2) * LinBasis(3) * LinBasis(4)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(2,j) * &
                                    LinBasis(1) * LinBasis(3) * LinBasis(4)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(3,j) * &
                                    LinBasis(1) * LinBasis(2) * LinBasis(4)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(4,j) * &
                                    LinBasis(1) * LinBasis(2) * LinBasis(3)
                END DO
              END DO

           CASE(8)

              IF ( Element % TYPE % ElementCode == 808 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(808)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                  LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(7)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                        dLinBasisdx(1,j) * LinBasis(7)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                        dLinBasisdx(7,j) * LinBasis(1)
                END DO
              END DO

         CASE DEFAULT
 
              WRITE( Message, '(a,i4,a)' ) 'Bubbles for element: ', &
               Element % TYPE % ElementCode, ' are not implemented.'
              CALL Error( 'ElementInfo', Message )
              CALL Fatal( 'ElementInfo', 'Please use p-element basis instead.' )

         END SELECT
       END IF
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ElementInfo
!------------------------------------------------------------------------------
   
   ! SUBROUTINE ElementInfoVec_InitWork(m, n)
   !   IMPLICIT NONE

   !   INTEGER, INTENT(IN) :: m, n
   !   INTEGER :: allocstat

   !   allocstat = 0
   !   IF (.NOT. ALLOCATED(BasisWrk)) THEN
   !     ALLOCATE(BasisWrk(m,n), &
   !             dBasisdxWrk(m,n,3), &
   !             LtoGMapsWrk(m,3,3), &
   !             DetJWrk(m), &
   !             uWrk(m), vWrk(m), wWrk(m), STAT=allocstat)
   !   ELSE IF (SIZE(BasisWrk,1) /= m .OR. SIZE(BasisWrk,2) /= n) THEN
   !     DEALLOCATE(BasisWrk, dBasisdxWrk, LtoGMapsWrk, DetJWrk, uWrk, vWrk, wWrk)
   !     ALLOCATE(BasisWrk(m,n), &
   !             dBasisdxWrk(m,n,3), &
   !             LtoGMapsWrk(m,3,3), &
   !             DetJWrk(m), &
   !             uWrk(m), vWrk(m), wWrk(m), STAT=allocstat)
   !   END IF

   !   ! Check memory allocation status
   !   IF (allocstat /= 0) THEN
   !     CALL Error('ElementInfo_InitWork','Storage allocation for local element basis failed')
   !   END IF
   ! END SUBROUTINE ElementInfoVec_InitWork

   ! SUBROUTINE ElementInfoVec_FreeWork()
   !   IMPLICIT NONE

   !   IF (ALLOCATED(BasisWrk)) THEN
   !     DEALLOCATE(BasisWrk, dBasisdxWrk, LtoGMapsWrk, DetJWrk, uWrk, vWrk, wWrk)
   !   END IF
   ! END SUBROUTINE ElementInfoVec_FreeWork

! ElementInfoVec currently uses only P element definitions for basis
! functions, even for purely nodal elements. Support for standard nodal elements
! will be implemented in the future. 
!------------------------------------------------------------------------------
   FUNCTION ElementInfoVec( Element, Nodes, nc, u, v, w, detJ, nbmax, Basis, dBasisdx ) RESULT(retval)
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Element_t), TARGET :: Element    !< Element structure
     TYPE(Nodes_t)   :: Nodes              !< Element nodal coordinates.
     INTEGER, INTENT(IN) :: nc             !< Number of local coordinates to compute values of the basis function
     REAL(KIND=dp), POINTER CONTIG :: u(:)  !< 1st local coordinates at which to calculate the basis function.
     REAL(KIND=dp), POINTER CONTIG :: v(:)  !< 2nd local coordinates.
     REAL(KIND=dp), POINTER CONTIG :: w(:)  !< 3rd local coordinates.
     REAL(KIND=dp) CONTIG, INTENT(OUT) :: detJ(:) !< Square roots of determinants of element coordinate system metric at coordinates
     INTEGER, INTENT(IN) :: nbmax          !< Maximum number of basis functions to compute
     REAL(KIND=dp) CONTIG :: Basis(:,:)    !< Basis function values at (u,v,w)
     REAL(KIND=dp) CONTIG, OPTIONAL :: dBasisdx(:,:,:)    !< Global first derivatives of basis functions at (u,v,w)
     LOGICAL :: retval                             !< If .FALSE. element is degenerate. or if local storage allocation fails

     ! Internal work arrays (always needed)
     REAL(KIND=dp) :: uWrk(VECTOR_BLOCK_LENGTH), vWrk(VECTOR_BLOCK_LENGTH), wWrk(VECTOR_BLOCK_LENGTH)
     REAL(KIND=dp) :: BasisWrk(VECTOR_BLOCK_LENGTH,nbmax)
     REAL(KIND=dp) :: dBasisdxWrk(VECTOR_BLOCK_LENGTH,nbmax,3)
     REAL(KIND=dp) :: DetJWrk(VECTOR_BLOCK_LENGTH)
     REAL(KIND=dp) :: LtoGMapsWrk(VECTOR_BLOCK_LENGTH,3,3)
     
     INTEGER :: i
!DIR$ ATTRIBUTES ALIGN:64::uWrk, vWrk, wWrk, BasisWrk, dBasisdxWrk, DetJWrk, LtoGMapsWrk
     
     !------------------------------------------------------------------------------
     ! Special case, Element: POINT
     IF (Element % TYPE % ElementCODE == 101) THEN
       DetJ(1:nc) = REAL(1, dp)
       Basis(1:nc,1) = REAL(1, dp)
       IF (PRESENT(dBasisdx)) THEN
         DO i=1,nc
           dBasisdx(i,1,1) = REAL(0, dp)
         END DO
       END IF
       retval = .TRUE.
       RETURN
     END IF
     
     ! Set up workspace arrays 
     ! CALL ElementInfoVec_InitWork(VECTOR_BLOCK_LENGTH, nbmax)
     IF ( nbmax < Element % TYPE % NumberOfNodes ) THEN
       CALL Fatal('ElementInfoVec','Not enough storage to compute local element basis')
     END IF

     retval =  ElementInfoVec_ComputePElementBasis(Element,Nodes,nc,u,v,w,detJ,nbmax,Basis,&
           uWrk,vWrk,wWrk,BasisWrk,dBasisdxWrk,DetJWrk,LtoGmapsWrk,dBasisdx)
   END FUNCTION ElementInfoVec
     
   FUNCTION ElementInfoVec_ComputePElementBasis(Element, Nodes, nc, u, v, w, DetJ, nbmax, Basis, &
                                                uWrk, vWrk, wWrk, BasisWrk, dBasisdxWrk, &
                                                DetJWrk, LtoGmapsWrk, dBasisdx) RESULT(retval)
     IMPLICIT NONE
     TYPE(Element_t), TARGET :: Element    !< Element structure
     TYPE(Nodes_t)   :: Nodes              !< Element nodal coordinates.
     INTEGER, INTENT(IN) :: nc             !< Number of local coordinates to compute values of the basis function
     REAL(KIND=dp), POINTER CONTIG :: u(:)  !< 1st local coordinates at which to calculate the basis function.
     REAL(KIND=dp), POINTER CONTIG :: v(:)  !< 2nd local coordinates.
     REAL(KIND=dp), POINTER CONTIG :: w(:)  !< 3rd local coordinates.
     REAL(KIND=dp) CONTIG, INTENT(OUT) :: detJ(:) !< Square roots of determinants of element coordinate system metric at coordinates
     INTEGER, INTENT(IN) :: nbmax          !< Maximum number of basis functions to compute
     REAL(KIND=dp) CONTIG :: Basis(:,:)    !< Basis function values at (u,v,w)
     ! Internal work arrays
     REAL(KIND=dp) :: uWrk(VECTOR_BLOCK_LENGTH), vWrk(VECTOR_BLOCK_LENGTH), wWrk(VECTOR_BLOCK_LENGTH)
     REAL(KIND=dp) :: BasisWrk(VECTOR_BLOCK_LENGTH,nbmax)
     REAL(KIND=dp) :: dBasisdxWrk(VECTOR_BLOCK_LENGTH,nbmax,3)
     REAL(KIND=dp) :: DetJWrk(VECTOR_BLOCK_LENGTH)
     REAL(KIND=dp) :: LtoGMapsWrk(VECTOR_BLOCK_LENGTH,3,3)
     REAL(KIND=dp) CONTIG, OPTIONAL :: dBasisdx(:,:,:)    !< Global first derivatives of basis functions at (u,v,w)
     LOGICAL :: retval                             !< If .FALSE. element is degenerate. or if local storage allocation fails


     !------------------------------------------------------------------------------
     !    Local variables
     !------------------------------------------------------------------------------
     INTEGER :: EdgeDegree(H1Basis_MaxPElementEdges), &
           FaceDegree(H1Basis_MaxPElementFaces), &
           EdgeDirection(H1Basis_MaxPElementEdgeNodes,H1Basis_MaxPElementEdges), &
           FaceDirection(H1Basis_MaxPElementFaceNodes,H1Basis_MaxPElementFaces)

     INTEGER :: cdim, dim, i, j, k, l, ll, lln, ncl, ip, n, p, &
           nbp, nbdxp, allocstat, ncpad, EdgeMaxDegree, FaceMaxDegree

     LOGICAL :: invertBubble, elem
!DIR$ ATTRIBUTES ALIGN:64::EdgeDegree, FaceDegree
!DIR$ ATTRIBUTES ALIGN:64::EdgeDirection, FaceDirection
!DIR$ ASSUME_ALIGNED uWrk:64, vWrk:64, wWrk:64, BasisWrk:64, dBasisdxWrk:64, DetJWrk:64, LtoGMapsWrk:64

     retval = .TRUE.
     n    = Element % TYPE % NumberOfNodes
     dim  = Element % TYPE % DIMENSION
     cdim = CoordinateSystemDimension()


     ! Block the computation for large values of input points
     DO ll=1,nc,VECTOR_BLOCK_LENGTH
       lln = MIN(ll+VECTOR_BLOCK_LENGTH-1,nc)
       ncl = lln-ll+1

       ! Set number of computed basis functions
       nbp = 0
       nbdxp = 0

       ! Block copy input
       uWrk(1:ncl) = u(ll:lln)
       IF (cdim > 1) THEN
         vWrk(1:ncl) = v(ll:lln)
       END IF
       IF (cdim > 2) THEN
         wWrk(1:ncl) = w(ll:lln)
       END IF

       ! Compute local p element basis
       SELECT CASE (Element % Type % ElementCode)
         ! Element: LINE
       CASE (202)
         ! Compute nodal basis
         CALL H1Basis_LineNodal(ncl, uWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dLineNodal(ncl, uWrk, nbmax, dBasisdxWrk, nbdxp)

         ! Element bubble functions
         IF (Element % BDOFS > 0) THEN 
           ! For first round of blocked loop, compute edge direction
           IF (ll==1) THEN
             ! Compute P from bubble dofs
             P = Element % BDOFS + 1

             IF (Element % PDefs % isEdge .AND. &
                   Element % NodeIndexes(1)> Element % NodeIndexes(2)) THEN
               invertBubble = .TRUE.
             ELSE
               invertBubble = .FALSE.
             END IF
           END IF

           CALL H1Basis_LineBubbleP(ncl, uWrk, P, nbmax, BasisWrk, nbp, invertBubble)
           CALL H1Basis_dLineBubbleP(ncl, uWrk, P, nbmax, dBasisdxWrk, nbdxp, invertBubble)
         END IF

         ! Element: TRIANGLE
       CASE (303)
         ! Compute nodal basis
         CALL H1Basis_TriangleNodalP(ncl, uWrk, vWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dTriangleNodalP(ncl, uWrk, vWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             CALL GetElementMeshEdgeInfo(CurrentModel % Solver % Mesh, &
                   Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree > 1) THEN
             CALL H1Basis_TriangleEdgeP(ncl, uWrk, vWrk, EdgeDegree, nbmax, BasisWrk, &
                   nbp, EdgeDirection)
             CALL H1Basis_dTriangleEdgeP(ncl, uWrk, vWrk, EdgeDegree, nbmax, dBasisdxWrk, &
                   nbdxp, EdgeDirection)
           END IF
         END IF

         ! Element bubble functions
         IF (Element % BDOFS > 0) THEN 
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             ! Compute P from bubble dofs
             P = CEILING( ( 3.0d0+SQRT(1.0d0+8.0d0*(Element % BDOFS)) ) / 2.0d0 )

             IF (Element % PDefs % isEdge) THEN
               ! Get 2D face direction
               CALL H1Basis_GetFaceDirection(Element % Type % ElementCode, &
                     1, &
                     Element % NodeIndexes, &
                     FaceDirection)
             END IF
           END IF
           IF (Element % PDefs % isEdge) THEN
             CALL H1Basis_TriangleBubbleP(ncl, uWrk, vWrk, P, nbmax, BasisWrk, nbp, &
                   FaceDirection(1:3,1))
             CALL H1Basis_dTriangleBubbleP(ncl, uWrk, vWrk, P, nbmax, dBasisdxWrk, nbdxp, &
                   FaceDirection(1:3,1))
           ELSE
             CALL H1Basis_TriangleBubbleP(ncl, uWrk, vWrk, P, nbmax, BasisWrk, nbp)
             CALL H1Basis_dTriangleBubbleP(ncl, uWrk, vWrk, P, nbmax, dBasisdxWrk, nbdxp)
           END IF
         END IF

         ! QUADRILATERAL
       CASE (404)
         ! Compute nodal basis
         CALL H1Basis_QuadNodal(ncl, uWrk, vWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dQuadNodal(ncl, uWrk, vWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             CALL GetElementMeshEdgeInfo(CurrentModel % Solver % Mesh, &
                   Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree > 1) THEN
             CALL H1Basis_QuadEdgeP(ncl, uWrk, vWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                   EdgeDirection)
             CALL H1Basis_dQuadEdgeP(ncl, uWrk, vWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                   EdgeDirection)
           END IF
         END IF

         ! Element bubble functions
         IF (Element % BDOFS > 0) THEN 
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             ! Compute P from bubble dofs
             P = CEILING( ( 5.0d0+SQRT(1.0d0+8.0d0*(Element % BDOFS)) ) / 2.0d0 )

             IF (Element % PDefs % isEdge) THEN
               ! Get 2D face direction
               CALL H1Basis_GetFaceDirection(Element % Type % ElementCode, &
                     1, &
                     Element % NodeIndexes, &
                     FaceDirection)
             END IF
           END IF

           IF (Element % PDefs % isEdge) THEN
             CALL H1Basis_QuadBubbleP(ncl, uWrk, vWrk, P, nbmax, BasisWrk, nbp, &
                   FaceDirection(1:4,1))
             CALL H1Basis_dQuadBubbleP(ncl, uWrk, vWrk, P, nbmax, dBasisdxWrk, nbdxp, &
                   FaceDirection(1:4,1))
           ELSE
             CALL H1Basis_QuadBubbleP(ncl, uWrk, vWrk, P, nbmax, BasisWrk, nbp)
             CALL H1Basis_dQuadBubbleP(ncl, uWrk, vWrk, P, nbmax, dBasisdxWrk, nbdxp)
           END IF
         END IF

         ! TETRAHEDRON
       CASE (504)
         ! Compute nodal basis
         CALL H1Basis_TetraNodalP(ncl, uWrk, vWrk, wWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dTetraNodalP(ncl, uWrk, vWrk, wWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             ! Get polynomial degree of each edge
             EdgeMaxDegree = 0
             IF (CurrentModel % Solver % Mesh % MinEdgeDOFs == &
                   CurrentModel % Solver % Mesh % MaxEdgeDOFs) THEN
               EdgeMaxDegree = Element % BDOFs+1
               EdgeDegree(1:Element % Type % NumberOfFaces) = EdgeMaxDegree
             ELSE
               DO i=1,6
                 EdgeDegree(i) = CurrentModel % Solver % &
                       Mesh % Edges( Element % EdgeIndexes(i) ) % BDOFs + 1
                 EdgeMaxDegree = MAX(EdgeDegree(i),EdgeMaxDegree)
               END DO
             END IF

             ! Tetrahedral directions are enforced by tetra element types
             IF (EdgeMaxDegree > 1) THEN
               CALL H1Basis_GetTetraEdgeDirection(Element % PDefs % TetraType, EdgeDirection)
             END IF
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree > 1) THEN
             CALL H1Basis_TetraEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                   EdgeDirection)
             CALL H1Basis_dTetraEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                   EdgeDirection)
           END IF
         END IF

         IF (ASSOCIATED( Element % FaceIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! face directions
           IF (ll==1) THEN
             ! Get polynomial degree of each face
             FaceMaxDegree = 0
             IF (CurrentModel % Solver % Mesh % MinFaceDOFs == &
                   CurrentModel % Solver % Mesh % MaxFaceDOFs) THEN
               FaceMaxDegree = CurrentModel % Solver % Mesh % Faces( Element % FaceIndexes(1) ) % PDefs % P
               FaceDegree(1:Element % Type % NumberOfFaces) = FaceMaxDegree
             ELSE
               DO i=1,4
                 IF (CurrentModel % Solver % Mesh % &
                       Faces( Element % FaceIndexes(i) ) % BDOFs /= 0) THEN
                   FaceDegree(i) = CurrentModel % Solver % Mesh % &
                         Faces( Element % FaceIndexes(i) ) % PDefs % P
                   FaceMaxDegree = MAX(FaceDegree(i), FaceMaxDegree)
                 ELSE
                   FaceDegree(i) = 0
                 END IF
               END DO
             END IF

             ! Tetrahedral directions are enforced by tetra element types
             IF (FaceMaxDegree > 1) THEN
               CALL H1Basis_GetTetraFaceDirection(Element % PDefs % TetraType, FaceDirection)
             END IF
           END IF

           ! Compute basis function values
           IF (FaceMaxDegree > 1) THEN
             CALL H1Basis_TetraFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, BasisWrk, nbp, &
                   FaceDirection)
             CALL H1Basis_dTetraFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, dBasisdxWrk, nbdxp, &
                   FaceDirection)
           END IF
         END IF

         ! Element bubble functions
         IF (Element % BDOFS > 0) THEN 
           ! Compute P based on bubble dofs
           P=CEILING(1/3d0*(81*(Element % BDOFS) + &
                 3*SQRT(-3d0+729*(Element % BDOFS)**2))**(1/3d0) + &
                 1d0/(81*(Element % BDOFS)+ &
                 3*SQRT(-3d0+729*(Element % BDOFS)**2))**(1/3d0)+2)

           CALL H1Basis_TetraBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, BasisWrk, nbp)
           CALL H1Basis_dTetraBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, dBasisdxWrk, nbdxp)
         END IF

         ! WEDGE
       CASE (706)
         ! Compute nodal basis
         CALL H1Basis_WedgeNodalP(ncl, uWrk, vWrk, wWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dWedgeNodalP(ncl, uWrk, vWrk, wWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             CALL GetElementMeshEdgeInfo(CurrentModel % Solver % Mesh, &
                   Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree > 1) THEN
             CALL H1Basis_WedgeEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                   EdgeDirection)
             CALL H1Basis_dWedgeEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                   EdgeDirection)
           END IF
         END IF

         IF (ASSOCIATED( Element % FaceIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! face directions
           IF (ll==1) THEN
             CALL GetElementMeshFaceInfo(CurrentModel % Solver % Mesh, &
                   Element, FaceDegree, FaceDirection, FaceMaxDegree)
           END IF

           ! Compute basis function values
           IF (FaceMaxDegree > 1) THEN
             CALL H1Basis_WedgeFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, BasisWrk, nbp, &
                   FaceDirection)
             CALL H1Basis_dWedgeFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, dBasisdxWrk, nbdxp, &
                   FaceDirection)
           END IF
         END IF

         ! Element bubble functions
         IF (Element % BDOFS > 0) THEN 
           ! Compute P from bubble dofs
           P=CEILING(1/3d0*(81*(Element % BDOFS) + &
                 3*SQRT(-3d0+729*(Element % BDOFS)**2))**(1/3d0) + &
                 1d0/(81*(Element % BDOFS)+ &
                 3*SQRT(-3d0+729*(Element % BDOFS)**2))**(1/3d0)+3)

           CALL H1Basis_WedgeBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, BasisWrk, nbp)
           CALL H1Basis_dWedgeBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, dBasisdxWrk, nbdxp)
         END IF

         ! HEXAHEDRON
       CASE (808)
         ! Compute local basis
         CALL H1Basis_BrickNodal(ncl, uWrk, vWrk, wWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dBrickNodal(ncl, uWrk, vWrk, wWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             CALL GetElementMeshEdgeInfo(CurrentModel % Solver % Mesh, &
                   Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree > 1) THEN
             CALL H1Basis_BrickEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                   EdgeDirection)
             CALL H1Basis_dBrickEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                   EdgeDirection)
           END IF
         END IF

         IF (ASSOCIATED( Element % FaceIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! face directions
           IF (ll==1) THEN
             CALL GetElementMeshFaceInfo(CurrentModel % Solver % Mesh, &
                   Element, FaceDegree, FaceDirection, FaceMaxDegree)
           END IF

           ! Compute basis function values
           IF (FaceMaxDegree > 1) THEN
             CALL H1Basis_BrickFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, BasisWrk, nbp, &
                   FaceDirection)
             CALL H1Basis_dBrickFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, dBasisdxWrk, nbdxp, &
                   FaceDirection)
           END IF
         END IF

         ! Element bubble functions
         IF (Element % BDOFS > 0) THEN 
           ! Compute P from bubble dofs
           P=CEILING(1/3d0*(81*Element % BDOFS + &
                 3*SQRT(-3d0+729*Element % BDOFS**2))**(1/3d0) + &
                 1d0/(81*Element % BDOFS+3*SQRT(-3d0+729*Element % BDOFS**2))**(1/3d0)+4)

           CALL H1Basis_BrickBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, BasisWrk, nbp)
           CALL H1Basis_dBrickBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, dBasisdxWrk, nbdxp)
         END IF

       CASE DEFAULT
         WRITE( Message, '(a,i4,a)' ) 'Vectorized basis for element: ', &
               Element % TYPE % ElementCode, ' not implemented.'
         CALL Error( 'ElementInfoVec', Message )
         CALL Fatal( 'ElementInfoVec', 'ElementInfoVec is still does not include pyramids.' )
       END SELECT

       ! Copy basis function values to global array
       DO j=1,nbp
         DO i=1,ncl
           Basis(i+ll-1,j)=BasisWrk(i,j)
         END DO
       END DO

       !--------------------------------------------------------------
       ! Element (contravariant) metric and square root of determinant
       !--------------------------------------------------------------
       elem = ElementMetricVec( Element, Nodes, ncl, nbp, DetJWrk, &
             nbmax, dBasisdxWrk, LtoGMapsWrk )
       IF (.NOT. elem) THEN
         retval = .FALSE.
         RETURN
       END IF

       !_ELMER_OMP_SIMD
       DO i=1,ncl
         DetJ(i+ll-1)=DetJWrk(i)
       END DO

       ! Get global basis functions
       !--------------------------------------------------------------
       ! First derivatives
       IF (PRESENT(dBasisdx)) THEN
!DIR$ FORCEINLINE
         CALL ElementInfoVec_ElementBasisToGlobal(ncl, nbp, nbmax, dBasisdxWrk, dim, cdim, LtoGMapsWrk, ll, dBasisdx)
       END IF
     END DO ! Block over Gauss points
  CONTAINS
   
     SUBROUTINE GetElementMeshEdgeInfo(Mesh, Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
       IMPLICIT NONE
       
       TYPE(Mesh_t), INTENT(IN) :: Mesh
       TYPE(Element_t), INTENT(IN) :: Element
       INTEGER, INTENT(OUT) :: EdgeDegree(H1Basis_MaxPElementEdges), &
               EdgeDirection(H1Basis_MaxPElementEdgeNodes,H1Basis_MaxPElementEdges)
       INTEGER, INTENT(OUT) :: EdgeMaxDegree
       INTEGER :: i

       EdgeMaxDegree = 0
       IF (Mesh % MinEdgeDOFs == Mesh % MaxEdgeDOFs) THEN
          EdgeDegree(1:Element % Type % NumberOfEdges) = Mesh % MaxEdgeDOFs + 1
          EdgeMaxDegree = Mesh % MaxEdgeDOFs + 1
       ELSE
       ! Get polynomial degree of each edge separately
!DIR$ LOOP COUNT MAX=12
          DO i=1,Element % Type % NumberOfEdges
             EdgeDegree(i) = Mesh % Edges( Element % EdgeIndexes(i) ) % BDOFs + 1
             EdgeMaxDegree = MAX(EdgeDegree(i), EdgeMaxDegree)
          END DO
       END IF

       ! Get edge directions if needed
       IF (EdgeMaxDegree > 1) THEN
         CALL H1Basis_GetEdgeDirection(Element % Type % ElementCode, &
                                       Element % Type % NumberOfEdges, &
                                       Element % NodeIndexes, &
                                       EdgeDirection)
       END IF
     END SUBROUTINE GetElementMeshEdgeInfo
     
     SUBROUTINE GetElementMeshFaceInfo(Mesh, Element, FaceDegree, FaceDirection, FaceMaxDegree)
       IMPLICIT NONE
       
       TYPE(Mesh_t), INTENT(IN) :: Mesh
       TYPE(Element_t), INTENT(IN) :: Element
       INTEGER, INTENT(OUT) :: FaceDegree(H1Basis_MaxPElementFaces), &
               FaceDirection(H1Basis_MaxPElementFaceNodes,H1Basis_MaxPElementFaces)
       INTEGER, INTENT(OUT) :: FaceMaxDegree
       INTEGER :: i

       ! Get polynomial degree of each face
       FaceMaxDegree = 0
       IF (Mesh % MinFaceDOFs == Mesh % MaxFaceDOFs) THEN
          FaceMaxDegree = Mesh % Faces( Element % FaceIndexes(1) ) % PDefs % P
          FaceDegree(1:Element % Type % NumberOfFaces) = FaceMaxDegree
       ELSE
!DIR$ LOOP COUNT MAX=6
          DO i=1,Element % Type % NumberOfFaces
             IF (Mesh % Faces( Element % FaceIndexes(i) ) % BDOFs /= 0) THEN
                FaceDegree(i) = Mesh % Faces( Element % FaceIndexes(i) ) % PDefs % P
                FaceMaxDegree = MAX(FaceDegree(i), FaceMaxDegree)
             ELSE
                FaceDegree(i) = 0
             END IF
          END DO
       END IF

       ! Get face directions
       IF (FaceMaxDegree > 1) THEN
         CALL H1Basis_GetFaceDirection(Element % Type % ElementCode, &
                                       Element % Type % NumberOfFaces, &
                                       Element % NodeIndexes, &
                                       FaceDirection)
       END IF
     END SUBROUTINE GetElementMeshFaceInfo     
!------------------------------------------------------------------------------
   END FUNCTION ElementInfoVec_ComputePElementBasis
!------------------------------------------------------------------------------
   
   SUBROUTINE ElementInfoVec_ElementBasisToGlobal(npts, nbasis, nbmax, dLBasisdx, dim, cdim, LtoGMap, offset, dBasisdx)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: npts
     INTEGER, INTENT(IN) :: nbasis
     INTEGER, INTENT(IN) :: nbmax
     REAL(KIND=dp), INTENT(IN) :: dLBasisdx(VECTOR_BLOCK_LENGTH,nbmax,3)
     INTEGER, INTENT(IN) :: dim
     INTEGER, INTENT(IN) :: cdim
     REAL(KIND=dp), INTENT(IN) :: LtoGMap(VECTOR_BLOCK_LENGTH,3,3)
     INTEGER, INTENT(IN) :: offset
     REAL(KIND=dp) CONTIG :: dBasisdx(:,:,:)

     INTEGER :: i, j, l
!DIR$ ASSUME_ALIGNED dLBasisdx:64, LtoGMap:64

     ! Map local basis function to global
     SELECT CASE (dim)
     CASE(1)
       !DIR$ LOOP COUNT MAX=3
       DO j=1,cdim
         DO i=1,nbasis
           !_ELMER_OMP_SIMD
           DO l=1,npts
             dBasisdx(l+offset-1,i,j) = dLBasisdx(l,i,1)*LtoGMap(l,j,1)
           END DO
         END DO
       END DO
     CASE(2)
       !DIR$ LOOP COUNT MAX=3
       DO j=1,cdim
         DO i=1,nbasis
           !_ELMER_OMP_SIMD
           DO l=1,npts
             ! Map local basis function to global
             dBasisdx(l+offset-1,i,j) = dLBasisdx(l,i,1)*LtoGMap(l,j,1)+ &
                   dLBasisdx(l,i,2)*LtoGMap(l,j,2)
           END DO
         END DO
       END DO
     CASE(3)
       !DIR$ LOOP COUNT MAX=3
       DO j=1,cdim
         DO i=1,nbasis
           !_ELMER_OMP_SIMD
           DO l=1,npts
             ! Map local basis function to global
             dBasisdx(l+offset-1,i,j) = dLBasisdx(l,i,1)*LtoGMap(l,j,1)+ &
                   dLBasisdx(l,i,2)*LtoGMap(l,j,2)+ &
                   dLBasisdx(l,i,3)*LtoGMap(l,j,3)
           END DO
         END DO
       END DO
     END SELECT

   END SUBROUTINE ElementInfoVec_ElementBasisToGlobal

   
!------------------------------------------------------------------------------
!>  Returns just the size of the element at its center.
!>  providing a more economical way than calling ElementInfo. 
!------------------------------------------------------------------------------
   FUNCTION ElementSize( Element, Nodes ) RESULT ( detJ )

     TYPE(Element_t) :: Element
     TYPE(Nodes_t) :: Nodes
     REAL(KIND=dp) :: detJ

     REAL(KIND=dp) :: u,v,w
     REAL(KIND=dp), ALLOCATABLE :: Basis(:)
     INTEGER :: n,family
     LOGICAL :: Stat


     family = Element % TYPE % ElementCode / 100
     n = Element % TYPE % NumberOfNodes
     ALLOCATE( Basis(n) )

     SELECT CASE ( family )
       
       CASE ( 1 )
         DetJ = 1.0_dp
         RETURN

       CASE ( 2 )
         u = 0.0_dp
         v = 0.0_dp

       CASE ( 3 )
         u = 0.5_dp
         v = 0.5_dp
         
       CASE ( 4 )
         u = 0.0_dp
         v = 0.0_dp

       CASE ( 5 )
         u = 0.5_dp
         v = 0.5_dp
         w = 0.5_dp

       CASE ( 8 ) 
         u = 0.0_dp
         v = 0.0_dp
         w = 0.0_dp
         
       CASE DEFAULT
         CALL Fatal('ElementSize','Not implemented for elementtype')

       END SELECT

       Stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )

     END FUNCTION ElementSize
!------------------------------------------------------------------------------


!----------------------------------------------------------------------------------
!>  Return H(div)-conforming face element basis function values and their divergence 
!>  with respect to the reference element coordinates at a given point on the
!>  reference element. Here the basis for a real element K is constructed by  
!>  transforming the basis functions defined on the reference element k via the 
!>  Piola transformation. The data for performing the Piola transformation is also returned.
!>  Note that the reference element is chosen as in the p-approximation so that
!>  the reference element edges/faces have the same length/area. This choice simplifies 
!>  the associated assembly procedure.
!---------------------------------------------------------------------------------
     RECURSIVE FUNCTION FaceElementInfo( Element, Nodes, u, v, w, F, detF, &
          Basis, FBasis, DivFBasis, BDM, Dual, BasisDegree ) RESULT(stat)
!------------------------------------------------------------------------------
       IMPLICIT NONE

       TYPE(Element_t), TARGET :: Element     !< Element structure
       TYPE(Nodes_t) :: Nodes                 !< Data corresponding to the classic element nodes
       REAL(KIND=dp) :: u                     !< 1st reference element coordinate at which the basis functions are evaluated
       REAL(KIND=dp) :: v                     !< 2nd reference element coordinate
       REAL(KIND=dp) :: w                     !< 3rd reference element coordinate
       REAL(KIND=dp) :: F(3,3)                !< The gradient F=Grad f, with f the element map f:k->K
       REAL(KIND=dp) :: detF                  !< The determinant of the gradient matrix F
       REAL(KIND=dp) :: Basis(:)              !< Standard nodal basis functions evaluated at (u,v,w)
       REAL(KIND=dp) :: FBasis(:,:)           !< Face element basis functions spanning the reference element space   
       REAL(KIND=dp) :: DivFBasis(:)          !< The divergence of basis functions with respect to the local coordinates
       LOGICAL, OPTIONAL :: BDM               !< If .TRUE., a basis for BDM space is constructed
       LOGICAL, OPTIONAL :: Dual              !< If .TRUE., create an alternate dual basis
       INTEGER, OPTIONAL :: BasisDegree(:)    !< This a dummy parameter at the moment
       LOGICAL :: Stat                        !< Should be .FALSE. for a degenerate element but this is not yet checked
!-----------------------------------------------------------------------------------------------------------------
!      Local variables
!------------------------------------------------------------------------------------------------------------
       INTEGER :: n, dim, q, i, j, k, ni, nj, nk, A, B, C, D, I1, I2
       INTEGER :: FDofMap(4,3), DofsPerFace, FaceIndeces(4)
       REAL(KIND=dp) :: dLbasisdx(MAX(SIZE(Nodes % x),SIZE(Basis)),3), t1(3), t2(3), m(3), e(3), S, D1, D2
       REAL(KIND=dp) :: BDMBasis(12,3), BDMDivBasis(12), WorkBasis(2,3), WorkDivBasis(2)     
       INTEGER, POINTER :: EdgeMap(:,:), FaceMap(:,:), Ind(:)
       INTEGER, TARGET :: TetraFaceMap(4,3)
       INTEGER :: SquareFaceMap(4)
       LOGICAL :: RevertSign(4), RevertSign2(4), CheckSignReversions, CreateBDMBasis, Parallel
       LOGICAL :: CreateDualBasis
       TYPE(Mesh_t), POINTER :: Mesh
!-----------------------------------------------------------------------------------------------------
       Mesh => CurrentModel % Solver % Mesh
       Parallel = ASSOCIATED(Mesh % ParallelInfo % Interface)

       TetraFaceMap(1,:) = (/ 2, 1, 3 /)
       TetraFaceMap(2,:) = (/ 1, 2, 4 /)
       TetraFaceMap(3,:) = (/ 2, 3, 4 /) 
       TetraFaceMap(4,:) = (/ 3, 1, 4 /)
       !--------------------------------------------------------------
       ! Check whether BDM or dual basis functions should be created 
       !--------------------------------------------------------------
       CreateBDMBasis = .FALSE.
       IF ( PRESENT(BDM) ) CreateBDMBasis = BDM
       CreateDualBasis = .FALSE.
       IF ( PRESENT(Dual) ) CreateDualBasis = Dual
       !-----------------------------------------------------------------------------------------------------
       stat = .TRUE.
       Basis = 0.0d0
       FBasis = 0.0d0
       DivFBasis = 0.0d0
       F = 0.0d0

       dLbasisdx = 0.0d0      
       n = Element % TYPE % NumberOfNodes
       dim = Element % TYPE % DIMENSION

       IF ( Element % TYPE % ElementCode == 101 ) THEN
          detF = 1.0d0
          Basis(1) = 1.0d0
          RETURN
       END IF

       !-----------------------------------------------------------------------
       ! The standard nodal basis functions on the reference element and
       ! their derivatives with respect to the local coordinates. These define 
       ! the mapping of the reference element to an actual element on the 
       ! background mesh but are not the basis functions for face element approximation.
       ! Remark: Using reference elements having the faces of the same area
       ! simplifies the implementation of element assembly procedures.
       !-----------------------------------------------------------------------
       SELECT CASE(Element % TYPE % ElementCode / 100)
       CASE(3)
          DO q=1,n
             Basis(q) = TriangleNodalPBasis(q, u, v)
             dLBasisdx(q,1:2) = dTriangleNodalPBasis(q, u, v) 
          END DO
       CASE(4)
          DO q=1,n
             Basis(q) = QuadNodalPBasis(q, u, v)
             dLBasisdx(q,1:2) = dQuadNodalPBasis(q, u, v) 
          END DO
       CASE(5)
          DO q=1,n
             Basis(q) = TetraNodalPBasis(q, u, v, w)
             dLBasisdx(q,1:3) = dTetraNodalPBasis(q, u, v, w)
          END DO
       CASE DEFAULT
          CALL Fatal('ElementDescription::FaceElementInfo','Unsupported element type')
       END SELECT          

       !-----------------------------------------------------------------------
       ! Get data for performing the Piola transformation...
       !-----------------------------------------------------------------------
       stat = PiolaTransformationData(n, Element, Nodes, F, detF, dLBasisdx) 
       !------------------------------------------------------------------------
       ! ... in order to define the basis for the element space X(K) via 
       ! applying the Piola transformation as
       !    X(K) = { B | B = 1/(det F) F b(f^{-1}(x)) }
       ! with b giving the face element basis function on the reference element k,
       ! f mapping k to the actual element K, i.e. K = f(k) and F = Grad f. This 
       ! function returns the local basis functions b and their divergence (with respect
       ! to local coordinates) evaluated at the integration point. The effect of 
       ! the Piola transformation need to be considered when integrating, so we 
       ! shall return also the values of F and det F.
       !
       ! The construction of face element bases could be done in an alternate way for 
       ! triangles and tetrahedra, while the chosen approach has the benefit that
       ! it generalizes to other cases. For example general quadrilaterals may now 
       ! be handled in the same way.
       !---------------------------------------------------------------------------

       SELECT CASE(Element % TYPE % ElementCode / 100)
       CASE(3)
          !----------------------------------------------------------------
          ! Note that the global orientation of face normal is taken to be
          ! n = t x e_z where the tangent vector t is aligned with
          ! the element edge and points towards the node that has
          ! a larger global index.
          !---------------------------------------------------------------
          EdgeMap => LGetEdgeMap(3)
          !EdgeMap => GetEdgeMap(GetElementFamily(Element))

          IF (CreateBDMBasis) THEN
             !-------------------------------------------------
             ! The following is for the BDM space of degree k=1.
             ! First two basis functions defined on face 12.
             !-------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
                ! The sign and order of basis functions are reverted as
                ! compared with the other possibility
                FBasis(1,1) = -sqrt(3.0d0)/6.0d0 * (sqrt(3.0d0) + u - v)             
                FBasis(1,2) = -sqrt(3.0d0)/6.0d0 * (-sqrt(3.0d0) - 3.0d0 * u + v)
                DivFBasis(1) = -sqrt(3.0d0)/3.0d0 

                FBasis(2,1) = -sqrt(3.0d0)/6.0d0 * (-sqrt(3.0d0) + u + v)             
                FBasis(2,2) = -sqrt(3.0d0)/6.0d0 * (-sqrt(3.0d0) + 3.0d0 * u + v)
                DivFBasis(2) = -sqrt(3.0d0)/3.0d0
             ELSE
                FBasis(1,1) = sqrt(3.0d0)/6.0d0 * (-sqrt(3.0d0) + u + v)             
                FBasis(1,2) = sqrt(3.0d0)/6.0d0 * (-sqrt(3.0d0) + 3.0d0 * u + v)
                DivFBasis(1) = sqrt(3.0d0)/3.0d0

                FBasis(2,1) = sqrt(3.0d0)/6.0d0 * (sqrt(3.0d0) + u - v)             
                FBasis(2,2) = sqrt(3.0d0)/6.0d0 * (-sqrt(3.0d0) - 3.0d0 * u + v)
                DivFBasis(2) = sqrt(3.0d0)/3.0d0
             END IF

             ! Two basis functions defined on face 23
             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
                FBasis(3,1) = -1.0d0/6.0d0 * (-3.0d0+sqrt(3.0d0)+(-3.0d0+sqrt(3.0d0))*u + 2.0d0*sqrt(3.0d0)*v)
                FBasis(3,2) = -1.0d0/6.0d0 * ( 3.0d0+sqrt(3.0d0) ) * v
                DivFBasis(3) = -sqrt(3.0d0)/3.0d0

                FBasis(4,1) = -1.0d0/(3.0d0+sqrt(3.0d0)) * (2.0d0+sqrt(3.0d0)+(2.0d0+sqrt(3.0d0))*u-(1.0d0+sqrt(3.0d0))*v)
                FBasis(4,2) = -1.0d0/6.0d0 * ( -3.0d0+sqrt(3.0d0) ) * v
                DivFBasis(4) = -sqrt(3.0d0)/3.0d0
             ELSE
                FBasis(3,1) = 1.0d0/(3.0d0+sqrt(3.0d0)) * (2.0d0+sqrt(3.0d0)+(2.0d0+sqrt(3.0d0))*u-(1.0d0+sqrt(3.0d0))*v)
                FBasis(3,2) = 1.0d0/6.0d0 * ( -3.0d0+sqrt(3.0d0) ) * v
                DivFBasis(3) = sqrt(3.0d0)/3.0d0

                FBasis(4,1) = 1.0d0/6.0d0 * (-3.0d0+sqrt(3.0d0)+(-3.0d0+sqrt(3.0d0))*u + 2.0d0*sqrt(3.0d0)*v)
                FBasis(4,2) = 1.0d0/6.0d0 * ( 3.0d0+sqrt(3.0d0) ) * v
                DivFBasis(4) = sqrt(3.0d0)/3.0d0
             END IF

             ! Two basis functions defined on face 31
             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
                FBasis(5,1) = -1.0d0/6.0d0 * (-3.0d0-sqrt(3.0d0)+(3.0d0+sqrt(3.0d0))*u + 2.0d0*sqrt(3.0d0)*v)
                FBasis(5,2) = -1.0d0/6.0d0 * ( -3.0d0+sqrt(3.0d0) ) * v
                DivFBasis(5) = -sqrt(3.0d0)/3.0d0

                FBasis(6,1) = -1.0d0/( 3.0d0+sqrt(3.0d0) ) * ( 1.0d0 - u - v - sqrt(3.0d0)*v )
                FBasis(6,2) = -( 3.0d0+2.0d0*sqrt(3.0d0) ) * v /(3.0d0*(1.0d0+sqrt(3.0d0)))
                DivFBasis(6) = -sqrt(3.0d0)/3.0d0
             ELSE
                FBasis(5,1) = 1.0d0/( 3.0d0+sqrt(3.0d0) ) * ( 1.0d0 - u - v - sqrt(3.0d0)*v ) 
                FBasis(5,2) = ( 3.0d0+2.0d0*sqrt(3.0d0) ) * v /(3.0d0*(1.0d0+sqrt(3.0d0)))
                DivFBasis(5) = sqrt(3.0d0)/3.0d0

                FBasis(6,1) = 1.0d0/6.0d0 * (-3.0d0-sqrt(3.0d0)+(3.0d0+sqrt(3.0d0))*u + 2.0d0*sqrt(3.0d0)*v)
                FBasis(6,2) = 1.0d0/6.0d0 * ( -3.0d0+sqrt(3.0d0) ) * v
                DivFBasis(6) = sqrt(3.0d0)/3.0d0
             END IF

          ELSE

             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(1,1) = SQRT(3.0d0)/6.0d0 * u
             FBasis(1,2) = -0.5d0 + SQRT(3.0d0)/6.0d0 * v
             DivFBasis(1) =  SQRT(3.0d0)/3.0d0
             IF (nj<ni) THEN
                FBasis(1,:) = -FBasis(1,:)
                DivFBasis(1) = -DivFBasis(1)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(2,1) = SQRT(3.0d0)/6.0d0 * (1.0d0 + u)
             FBasis(2,2) = SQRT(3.0d0)/6.0d0 * v
             DivFBasis(2) =  SQRT(3.0d0)/3.0d0        
             IF (nj<ni) THEN
                FBasis(2,:) = -FBasis(2,:)
                DivFBasis(2) = -DivFBasis(2)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(3,1) = SQRT(3.0d0)/6.0d0 * (-1.0d0 + u)
             FBasis(3,2) = SQRT(3.0d0)/6.0d0 * v
             DivFBasis(3) =  SQRT(3.0d0)/3.0d0          
             IF (nj<ni) THEN
                FBasis(3,:) = -FBasis(3,:)
                DivFBasis(3) = -DivFBasis(3)
             END IF

          END IF
          
       CASE(4)
          !--------------------------------------------------------------------
          ! Quadrilateral Arnold-Boffi-Falk (ABF) element basis of degree k=0
          !--------------------------------------------------------------------
          EdgeMap => LGetEdgeMap(4)
          SquareFaceMap(:) = (/ 1,2,3,4 /)          
          Ind => Element % Nodeindexes

          IF (.NOT. CreateDualBasis) THEN
             !-------------------------------------------------
             ! Four basis functions defined on the edges
             !-------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(1,1) = 0.0d0
             FBasis(1,2) = -((-1.0d0 + v)*v)/4.0d0
             DivFBasis(1) = (1.0d0 - 2*v)/4.0d0
             IF (nj<ni) THEN
                FBasis(1,:) = -FBasis(1,:)
                DivFBasis(1) = -DivFBasis(1)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(2,1) = (u*(1.0d0 + u))/4.0d0
             FBasis(2,2) = 0.0d0
             DivFBasis(2) = (1 + 2.0d0*u)/4.0d0
             IF (nj<ni) THEN
                FBasis(2,:) = -FBasis(2,:)
                DivFBasis(2) = -DivFBasis(2)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(3,1) = 0.0d0
             FBasis(3,2) = (v*(1.0d0 + v))/4.0d0
             DivFBasis(3) = (1.0d0 + 2.0d0*v)/4.0d0
             IF (nj<ni) THEN
                FBasis(3,:) = -FBasis(3,:)
                DivFBasis(3) = -DivFBasis(3)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(4,1) = -((-1.0d0 + u)*u)/4.0d0
             FBasis(4,2) = 0.0d0
             DivFBasis(4) = (1.0d0 - 2.0d0*u)/4.0d0
             IF (nj<ni) THEN
                FBasis(4,:) = -FBasis(4,:)
                DivFBasis(4) = -DivFBasis(4)
             END IF

             !--------------------------------------------------------------------
             ! Additional two basis functions associated with the element interior
             !-------------------------------------------------------------------
             WorkBasis(1,:) = 0.0d0
             WorkBasis(2,:) = 0.0d0
             WorkDivBasis(:) = 0.0d0

             WorkBasis(1,1) = 0.0d0
             WorkBasis(1,2) = (-1.0d0 + v**2)/2.0d0
             WorkDivBasis(1) = v

             WorkBasis(2,1) = (1.0d0 - u**2)/2.0d0
             WorkBasis(2,2) = 0.0d0
             WorkDivBasis(2) = -u

             DO j=1,4
                FaceIndeces(j) = Ind(SquareFaceMap(j))
             END DO
             IF (Parallel) THEN
                DO j=1,4
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             FBasis(5,:) = D1 * WorkBasis(I1,:)
             DivFBasis(5) = D1 * WorkDivBasis(I1)
             FBasis(6,:) = D2 * WorkBasis(I2,:)
             DivFBasis(6) = D2 * WorkDivBasis(I2)   
          ELSE
             !---------------------------------------------------------------------------
             ! Create alternate basis functions for the ABF space so that these basis
             ! functions are dual to the standard basis functions when the mesh is regular.
             ! First four basis functions which are dual to the standard edge basis 
             ! functions:
             !----------------------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(1,1) = 0.0d0
             FBasis(1,2) = (-3.0d0*(-1.0d0 - 2.0d0*v + 5.0d0*v**2))/4.0d0
             DivFBasis(1) = (-3.0d0*(-1.0d0 + 5.0d0*v))/2.0d0
             IF (nj<ni) THEN
                FBasis(1,:) = -FBasis(1,:)
                DivFBasis(1) = -DivFBasis(1)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(2,1) = (3.0d0*(-1.0d0 + 2.0d0*u + 5.0d0*u**2))/4.0d0
             FBasis(2,2) = 0.0d0
             DivFBasis(2) = (3.0d0*(1.0d0 + 5.0d0*u))/2.0d0
             IF (nj<ni) THEN
                FBasis(2,:) = -FBasis(2,:)
                DivFBasis(2) = -DivFBasis(2)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(3,1) = 0.0d0
             FBasis(3,2) = (3.0d0*(-1.0d0 + 2.0d0*v + 5.0d0*v**2))/4.0d0
             DivFBasis(3) = (3.0d0*(1.0d0 + 5.0d0*v))/2.0d0
             IF (nj<ni) THEN
                FBasis(3,:) = -FBasis(3,:)
                DivFBasis(3) = -DivFBasis(3)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             FBasis(4,1) = (-3.0d0*(-1.0d0 - 2.0d0*u + 5.0d0*u**2))/4.0d0
             FBasis(4,2) = 0.0d0
             DivFBasis(4) = (-3.0d0*(-1.0d0 + 5.0d0*u))/2.0d0
             IF (nj<ni) THEN
                FBasis(4,:) = -FBasis(4,:)
                DivFBasis(4) = -DivFBasis(4)
             END IF

             !-------------------------------------------------------------------------
             ! Additional two dual basis functions associated with the element interior
             !-------------------------------------------------------------------------
             WorkBasis(1,:) = 0.0d0
             WorkBasis(2,:) = 0.0d0
             WorkDivBasis(:) = 0.0d0

             WorkBasis(1,1) = 0.0d0
             WorkBasis(1,2) = (3.0d0*(-3.0d0 + 5.0d0*v**2))/8.0d0
             WorkDivBasis(1) = 15.0d0*v/4.0d0

             WorkBasis(2,1) = (3.0d0*(3.0d0 - 5.0d0*u**2))/8.0d0
             WorkBasis(2,2) = 0.0d0
             WorkDivBasis(2) = -15.0d0*u/4.0d0

             DO j=1,4
                FaceIndeces(j) = Ind(SquareFaceMap(j))
             END DO
             IF (Parallel) THEN
                DO j=1,4
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             FBasis(5,:) = D1 * WorkBasis(I1,:)
             DivFBasis(5) = D1 * WorkDivBasis(I1)
             FBasis(6,:) = D2 * WorkBasis(I2,:)
             DivFBasis(6) = D2 * WorkDivBasis(I2)
          END IF

       CASE(5)
          !-----------------------------------------
          ! This branch is for handling tetrahedra
          !-----------------------------------------
          FaceMap => TetraFaceMap
          Ind => Element % Nodeindexes

          DO q=1,4
             !-----------------------------------------------------------------------------------
             ! Check first whether a sign reversion will be needed as face dofs have orientation.
             ! If the sign is not reverted, a positive value of the degree of freedom produces
             ! positive outward flux from the element through the face handled.
             !-----------------------------------------------------------------------------------
             RevertSign(q) = .FALSE.

             DO j=1,3
                FaceIndeces(j) = Ind(FaceMap(q,j))
             END DO
             IF (Parallel) THEN
                DO j=1,3
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                END DO
             END IF
             
             IF ( (FaceIndeces(1) < FaceIndeces(2)) .AND. (FaceIndeces(1) < FaceIndeces(3)) ) THEN
                IF (FaceIndeces(3) < FaceIndeces(2)) THEN
                   RevertSign(q) = .TRUE.
                END IF
             ELSE IF ( ( FaceIndeces(2) < FaceIndeces(1) ) .AND. ( FaceIndeces(2) < FaceIndeces(3) ) ) THEN
                IF ( FaceIndeces(1) < FaceIndeces(3) ) THEN
                   RevertSign(q) = .TRUE.
                END IF
             ELSE  
                IF ( FaceIndeces(2) < FaceIndeces(1) ) THEN
                   RevertSign(q) = .TRUE.
                END IF
             END IF

          END DO

          !----------------------------------------------------------------------
          ! Another way for finding sign reversions. This code is retained here,
          ! although it was used for verification purposes...
          !----------------------------------------------------------------------
          CheckSignReversions = .FALSE.
          IF (CheckSignReversions) THEN
             DO q=1,4
                RevertSign2(q) = .FALSE.
                i = FaceMap(q,1)
                j = FaceMap(q,2)
                k = FaceMap(q,3)

                IF ( ( Ind(i) < Ind(j) ) .AND. ( Ind(i) < Ind(k) ) ) THEN
                   A = i
                   IF (Ind(j) < Ind(k)) THEN
                      B = j
                      C = k
                   ELSE
                      B = k
                      C = j
                   END IF
                ELSE IF ( ( Ind(j) < Ind(i) ) .AND. ( Ind(j) < Ind(k) ) ) THEN
                   A = j
                   IF (Ind(i) < Ind(k)) THEN
                      B = i
                      C = k
                   ELSE
                      B = k
                      C = i
                   END IF
                ELSE
                   A = k
                   IF (Ind(i) < Ind(j)) THEN
                      B = i
                      C = j
                   ELSE
                      B = j
                      C = i
                   END IF
                END IF

                t1(1) = Nodes % x(B) - Nodes % x(A)
                t1(2) = Nodes % y(B) - Nodes % y(A)              
                t1(3) = Nodes % z(B) - Nodes % z(A)

                t2(1) = Nodes % x(C) - Nodes % x(A)
                t2(2) = Nodes % y(C) - Nodes % y(A)              
                t2(3) = Nodes % z(C) - Nodes % z(A)

                m(1:3) = CrossProduct(t1,t2)

                SELECT CASE(q)
                CASE(1)
                   D = 4
                CASE(2)
                   D = 3 
                CASE(3)
                   D = 1
                CASE(4)
                   D = 2                   
                END SELECT

                e(1) = Nodes % x(D) - Nodes % x(A)
                e(2) = Nodes % y(D) - Nodes % y(A)                
                e(3) = Nodes % z(D) - Nodes % z(A)  

                IF ( SUM(m(1:3) * e(1:3)) > 0.0d0 ) RevertSign2(q) = .TRUE.

             END DO

             IF ( ANY(RevertSign(1:4) .NEQV. RevertSign2(1:4)) ) THEN
                PRINT *, 'CONFLICTING SIGN REVERSIONS SUGGESTED'
                PRINT *, RevertSign(1:4)
                PRINT *, RevertSign2(1:4)
                STOP
             END IF

          END IF

          IF (CreateBDMBasis) THEN
             DofsPerFace = 3 ! This choice is used for the BDM space of degree k=1
             !----------------------------------------------------------------------------
             ! Create a table of BDM basis functions in the default order
             !----------------------------------------------------------------------------
             ! Face {213}:
             BDMBasis(1,1) = (3*Sqrt(6.0d0) + 2*Sqrt(6.0d0)*u - 3*Sqrt(2.0d0)*v - 3*w)/12.0
             BDMBasis(1,2) = (-2*Sqrt(2.0d0) - 3*Sqrt(2.0d0)*u + Sqrt(3.0d0)*w)/12.0
             BDMBasis(1,3) = (-8 - 12*u + 4*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w)/12.0

             BDMBasis(2,1) = (2*Sqrt(6.0d0)*u + 3*(-Sqrt(6.0d0) + Sqrt(2.0d0)*v + w))/12.0
             BDMBasis(2,2) = (-2*Sqrt(2.0d0) + 3*Sqrt(2.0d0)*u + Sqrt(3.0d0)*w)/12.0
             BDMBasis(2,3) = u + (-8 + 4*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w)/12.0

             BDMBasis(3,1) = -u/(2.0*Sqrt(6.0d0))
             BDMBasis(3,2) = (Sqrt(2.0d0) + 3*Sqrt(6.0d0)*v - 2*Sqrt(3.0d0)*w)/12.0
             BDMBasis(3,3) = (4 - 8*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w)/12.0

             ! Face {124}:
             BDMBasis(4,1) = (2*Sqrt(6.0d0)*u + 3*(-Sqrt(6.0d0) + Sqrt(2.0d0)*v + w))/12.0
             BDMBasis(4,2) = (-6*Sqrt(2.0d0) + 9*Sqrt(2.0d0)*u + 2*Sqrt(6.0d0)*v + 3*Sqrt(3.0d0)*w)/12.0
             BDMBasis(4,3) = -w/(2.0*Sqrt(6.0d0))
             BDMBasis(5,1) = (3*Sqrt(6.0d0) + 2*Sqrt(6.0d0)*u - 3*Sqrt(2.0d0)*v - 3*w)/12.0
             BDMBasis(5,2) = (-6*Sqrt(2.0d0) - 9*Sqrt(2.0d0)*u + 2*Sqrt(6.0d0)*v + 3*Sqrt(3.0d0)*w)/12.0
             BDMBasis(5,3) = -w/(2.0*Sqrt(6.0d0))
             BDMBasis(6,1) = -u/(2.0*Sqrt(6.0d0))
             BDMBasis(6,2) = (3*Sqrt(2.0d0) - Sqrt(6.0d0)*v - 6*Sqrt(3.0d0)*w)/12.0
             BDMBasis(6,3) = (5*w)/(2.0*Sqrt(6.0d0))

             ! Face {234}:
             BDMBasis(7,1) = (5*Sqrt(6.0d0) + 5*Sqrt(6.0d0)*u - 6*Sqrt(2.0d0)*v - 6*w)/12.0
             BDMBasis(7,2) = -v/(2.0*Sqrt(6.0d0))
             BDMBasis(7,3) = -w/(2.0*Sqrt(6.0d0))
             BDMBasis(8,1) = (-Sqrt(6.0d0) - Sqrt(6.0d0)*u + 6*Sqrt(2.0d0)*v - 3*w)/12.0
             BDMBasis(8,2) = (5*Sqrt(6.0)*v - 3*Sqrt(3.0d0)*w)/12.0
             BDMBasis(8,3) = -w/(2.0*Sqrt(6.0d0))
             BDMBasis(9,1) = (-Sqrt(6.0d0) - Sqrt(6.0d0)*u + 9*w)/12.0
             BDMBasis(9,2) = (-(Sqrt(6.0d0)*v) + 3*Sqrt(3.0d0)*w)/12.0
             BDMBasis(9,3) = (5*w)/(2.0*Sqrt(6.0d0))

             ! Face {314}:             
             BDMBasis(10,1) = (Sqrt(6.0d0) - Sqrt(6.0d0)*u - 6*Sqrt(2.0d0)*v + 3*w)/12.0
             BDMBasis(10,2) = (5*Sqrt(6.0d0)*v - 3*Sqrt(3.0d0)*w)/12.0
             BDMBasis(10,3) = -w/(2.0*Sqrt(6.0d0))
             BDMBasis(11,1) = (-5*Sqrt(6.0d0) + 5*Sqrt(6.0d0)*u + 6*Sqrt(2.0d0)*v + 6*w)/12.0
             BDMBasis(11,2) = -v/(2.0*Sqrt(6.0d0))
             BDMBasis(11,3) = -w/(2.0*Sqrt(6.0d0))
             BDMBasis(12,1) = (Sqrt(6.0d0) - Sqrt(6.0d0)*u - 9*w)/12.0
             BDMBasis(12,2) = (-(Sqrt(6.0d0)*v) + 3*Sqrt(3.0d0)*w)/12.0
             BDMBasis(12,3) = (5*w)/(2.0*Sqrt(6.0d0))

             !----------------------------------------------------------------------
             ! Find out how face basis functions must be ordered so that the global
             ! indexing convention is respected. 
             !-----------------------------------------------------------------------
             FDofMap = 0
             DO q=1,4
                !i = FaceMap(q,1)
                !j = FaceMap(q,2)
                !k = FaceMap(q,3)

                DO j=1,3
                   FaceIndeces(j) = Ind(FaceMap(q,j))
                END DO
                IF (Parallel) THEN
                   DO j=1,3
                      FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                   END DO
                END IF
                
                IF ( ( FaceIndeces(1) < FaceIndeces(2) ) .AND. ( FaceIndeces(1) < FaceIndeces(3) ) ) THEN
                   ! A = i
                   FDofMap(q,1) = 1
                   IF (FaceIndeces(2) < FaceIndeces(3)) THEN
                      !B = j
                      !C = k
                      FDofMap(q,2) = 2
                      FDofMap(q,3) = 3                      
                   ELSE
                      !B = k
                      !C = j
                      FDofMap(q,2) = 3
                      FDofMap(q,3) = 2
                   END IF
                ELSE IF ( ( FaceIndeces(2) < FaceIndeces(1) ) .AND. ( FaceIndeces(2) < FaceIndeces(3) ) ) THEN
                   !A = j
                   FDofMap(q,1) = 2
                   IF (FaceIndeces(1) < FaceIndeces(3)) THEN
                      !B = i
                      !C = k
                      FDofMap(q,2) = 1
                      FDofMap(q,3) = 3
                   ELSE
                      !B = k
                      !C = i
                      FDofMap(q,2) = 3
                      FDofMap(q,3) = 1
                   END IF
                ELSE
                   !A = k
                   FDofMap(q,1) = 3
                   IF (FaceIndeces(1) < FaceIndeces(2)) THEN
                      !B = i
                      !C = j
                      FDofMap(q,2) = 1
                      FDofMap(q,3) = 2 
                   ELSE
                      !B = j
                      !C = i
                      FDofMap(q,2) = 2
                      FDofMap(q,3) = 1 
                   END IF
                END IF
             END DO

             !-----------------------------------------------------
             ! Now do the actual reordering and sign reversion
             !-----------------------------------------------------
             DO q=1,4
                IF (RevertSign(q)) THEN
                   S = -1.0d0
                ELSE
                   S = 1.0d0
                END IF

                DO j=1,DofsPerFace
                   k = FDofMap(q,j)
                   i = (q-1)*DofsPerFace + j
                   FBasis(i,:) = S * BDMBasis((q-1)*DofsPerFace+k,:)
                   DivFBasis(i) = S * sqrt(3.0d0)/(2.0d0*sqrt(2.0d0))
                END DO
             END DO

          ELSE
             !-------------------------------------------------------------------------
             ! The basis functions that define RT space on reference element
             !-----------------------------------------------------------------------
             FBasis(1,1) = SQRT(2.0d0)/4.0d0 * u
             FBasis(1,2) = -SQRT(6.0d0)/12.0d0 + SQRT(2.0d0)/4.0d0 * v
             FBasis(1,3) = -1.0d0/SQRT(3.0d0) + SQRT(2.0d0)/4.0d0 * w
             DivFBasis(1) = 3.0d0*SQRT(2.0d0)/4.0d0
             IF ( RevertSign(1) ) THEN
                FBasis(1,:) = -FBasis(1,:)
                DivFBasis(1) = -DivFBasis(1)
             END IF

             FBasis(2,1) = SQRT(2.0d0)/4.0d0 * u
             FBasis(2,2) = -SQRT(6.0d0)/4.0d0 + SQRT(2.0d0)/4.0d0 * v
             FBasis(2,3) = SQRT(2.0d0)/4.0d0 * w
             DivFBasis(2) = 3.0d0*SQRT(2.0d0)/4.0d0
             IF ( RevertSign(2) ) THEN
                FBasis(2,:) = -FBasis(2,:)
                DivFBasis(2) = -DivFBasis(2)
             END IF

             FBasis(3,1) = SQRT(2.0d0)/4.0d0 + SQRT(2.0d0)/4.0d0 * u
             FBasis(3,2) = SQRT(2.0d0)/4.0d0 * v
             FBasis(3,3) = SQRT(2.0d0)/4.0d0 * w
             DivFBasis(3) = 3.0d0*SQRT(2.0d0)/4.0d0
             IF ( RevertSign(3) ) THEN
                FBasis(3,:) = -FBasis(3,:)
                DivFBasis(3) = -DivFBasis(3)
             END IF

             FBasis(4,1) = -SQRT(2.0d0)/4.0d0 + SQRT(2.0d0)/4.0d0 * u
             FBasis(4,2) = SQRT(2.0d0)/4.0d0 * v
             FBasis(4,3) = SQRT(2.0d0)/4.0d0 * w
             DivFBasis(4) = 3.0d0*SQRT(2.0d0)/4.0d0
             IF ( RevertSign(4) ) THEN
                FBasis(4,:) = -FBasis(4,:)
                DivFBasis(4) = -DivFBasis(4)
             END IF
          END IF
       CASE DEFAULT
          CALL Fatal('ElementDescription::FaceElementInfo','Unsupported element type')
       END SELECT
    
!-----------------------------------------------------------------------------
     END FUNCTION FaceElementInfo
!------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------
!> This function returns data for performing the Piola transformation 
!------------------------------------------------------------------------------------------------
     FUNCTION PiolaTransformationData(nn,Element,Nodes,F,DetF,dLBasisdx) RESULT(Success)
!-------------------------------------------------------------------------------------------------
       INTEGER :: nn                   !< The number of classic nodes used in the element mapping
       TYPE(Element_t) :: Element      !< Element structure
       TYPE(Nodes_t) :: Nodes          !< Data corresponding to the classic element nodes
       REAL(KIND=dp) :: F(:,:)         !< The gradient of the element mapping
       REAL(KIND=dp) :: DetF           !< The determinant of the gradient matrix (or the Jacobian matrix)
       REAL(KIND=dp) :: dLBasisdx(:,:) !< Derivatives of nodal basis functions with respect to local coordinates
       LOGICAL :: Success              !< Could and should return .FALSE. if the element is degenerate
!-----------------------------------------------------------------------------------------------------
!      Local variables
!-------------------------------------------------------------------------------------------------
       REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
       INTEGER :: cdim,dim,n,i
!-------------------------------------------------------------------------------------------------
       x => Nodes % x
       y => Nodes % y
       z => Nodes % z     

       ! cdim = CoordinateSystemDimension()
       n = MIN( SIZE(x), nn )
       dim  = Element % TYPE % DIMENSION

       !------------------------------------------------------------------------------
       ! The gradient of the element mapping K = f(k), with k the reference element
       !------------------------------------------------------------------------------
       F = 0.0d0
       DO i=1,dim
          F(1,i) = SUM( x(1:n) * dLBasisdx(1:n,i) )
          F(2,i) = SUM( y(1:n) * dLBasisdx(1:n,i) )
          !IF (dim == 3) &
          ! In addition to the case dim = 3, the following entries may be useful  
          ! with dim=2 when natural BCs in 3-D are handled. 
          F(3,i) = SUM( z(1:n) * dLBasisdx(1:n,i) )
       END DO

       SELECT CASE( dim )    
       CASE (2)
          DetF = F(1,1)*F(2,2) - F(1,2)*F(2,1)
       CASE(3)
          DetF = F(1,1) * ( F(2,2)*F(3,3) - F(2,3)*F(3,2) ) + &
               F(1,2) * ( F(2,3)*F(3,1) - F(2,1)*F(3,3) ) + &
               F(1,3) * ( F(2,1)*F(3,2) - F(2,2)*F(3,1) )
       END SELECT

       success = .TRUE.
!------------------------------------------------
     END FUNCTION PiolaTransformationData
!------------------------------------------------


!------------------------------------------------------------------------------
!> Perform the cross product of two vectors
!------------------------------------------------------------------------------
     FUNCTION CrossProduct( v1, v2 ) RESULT( v3 )
!------------------------------------------------------------------------------
       IMPLICIT NONE
       REAL(KIND=dp) :: v1(3), v2(3), v3(3)
       v3(1) =  v1(2)*v2(3) - v1(3)*v2(2)
       v3(2) = -v1(1)*v2(3) + v1(3)*v2(1)
       v3(3) =  v1(1)*v2(2) - v1(2)*v2(1)
!------------------------------------------------------------------------------
     END FUNCTION CrossProduct
!------------------------------------------------------------------------------


!----------------------------------------------------------------------------------
!>  Return H(curl)-conforming edge element basis function values and their Curl  
!>  with respect to the reference element coordinates at a given point on the
!>  reference element. Here the basis for a real element K is constructed by  
!>  transforming the basis functions defined on the reference element k via a version
!>  of the Piola transformation designed for functions in H(curl). This construction
!>  differs from the approach taken in the alternate subroutine GetEdgeBasis, which
!>  does not make reference to the Piola transformation and hence may have limitations
!>  in its extendability. The data for performing the Piola transformation is also returned.
!>  Note that the reference element is chosen as in the p-approximation so that
!>  the reference element edges/faces have the same lenghth/area. This choice simplifies 
!>  the associated assembly procedure.
!>     With giving the optional argument ApplyPiolaTransform = .TRUE., this function
!>  also performs the Piola transform, so that the basis functions and their spatial
!>  curl as defined on the physical element are returned.
!>     In the lowest-order case this function returns the basis functions belonging
!>  to the optimal family which is not subject to degradation of convergence on
!>  meshes consting of non-affine physical elements. The second-order elements
!>  are members of the Nedelec's first family and are constructed in a hierarchic
!>  fashion (the lowest-order basis functions give a partial construction of
!>  the second-order basis).
!---------------------------------------------------------------------------------
     FUNCTION EdgeElementInfo( Element, Nodes, u, v, w, F, G, detF, &
          Basis, EdgeBasis, RotBasis, dBasisdx, SecondFamily, BasisDegree, &
          ApplyPiolaTransform, ReadyEdgeBasis, ReadyRotBasis, &
          TangentialTrMapping) RESULT(stat)
!------------------------------------------------------------------------------
       IMPLICIT NONE

       TYPE(Element_t), TARGET :: Element        !< Element structure
       TYPE(Nodes_t) :: Nodes                    !< Data corresponding to the classic element nodes
       REAL(KIND=dp) :: u                        !< 1st reference element coordinate at which the basis functions are evaluated
       REAL(KIND=dp) :: v                        !< 2nd local coordinate
       REAL(KIND=dp) :: w                        !< 3rd local coordinate
       REAL(KIND=dp), OPTIONAL :: F(3,3)         !< The gradient F=Grad f, with f the element map f:k->K
       REAL(KIND=dp), OPTIONAL :: G(3,3)         !< The transpose of the inverse of the gradient F
       REAL(KIND=dp) :: detF                     !< The determinant of the gradient matrix F
       REAL(KIND=dp) :: Basis(:)                 !< H1-conforming basis functions evaluated at (u,v,w)
       REAL(KIND=dp) :: EdgeBasis(:,:)           !< The basis functions b spanning the reference element space
       REAL(KIND=dp), OPTIONAL :: RotBasis(:,:)  !< The Curl of the edge basis functions with respect to the local coordinates
       REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)  !< The first derivatives of the H1-conforming basis functions at (u,v,w)
       LOGICAL, OPTIONAL :: SecondFamily         !< If .TRUE., a Nedelec basis of the second kind is returned (only simplicial elements)
       INTEGER, OPTIONAL :: BasisDegree          !< The approximation degree 2 is also supported
       LOGICAL, OPTIONAL :: ApplyPiolaTransform  !< If  .TRUE., perform the Piola transform so that, instead of b
                                                 !< and Curl b, return  B(f(p)) and (curl B)(f(p)) with B(x) the basis 
                                                 !< functions on the physical element and curl the spatial curl operator.
                                                 !< In this case the absolute value of detF is returned.
       REAL(KIND=dp), OPTIONAL :: ReadyEdgeBasis(:,:) !< A pretabulated edge basis function can be given
       REAL(KIND=dp), OPTIONAL :: ReadyRotBasis(:,:)  !< The preretabulated Curl of the edge basis function
       LOGICAL, OPTIONAL :: TangentialTrMapping  !< To return b x n, with n=(0,0,1) the normal to the 2D reference element.
                                                 !< The Piola transform is then the usual div-conforming version.    
       LOGICAL :: Stat                           !< .FALSE. for a degenerate element
!-----------------------------------------------------------------------------------------------------------------
!      Local variables
!------------------------------------------------------------------------------------------------------------
       TYPE(Mesh_t), POINTER :: Mesh
       INTEGER :: n, dim, cdim, q, i, j, k, l, ni, nj, A, I1, I2, FaceIndeces(4)
       REAL(KIND=dp) :: dLbasisdx(MAX(SIZE(Nodes % x),SIZE(Basis)),3), WorkBasis(4,3), WorkCurlBasis(4,3)
       REAL(KIND=dp) :: D1, D2, B(3), curlB(3), GT(3,3), LG(3,3), LF(3,3)
       REAL(KIND=dp) :: ElmMetric(3,3), detJ, CurlBasis(54,3)
       REAL(KIND=dp) :: t(3), s(3), v1, v2, v3, h1, h2, h3, dh1, dh2, dh3, grad(2)
       REAL(KIND=dp) :: LBasis(Element % TYPE % NumberOfNodes), Beta(4), EdgeSign(16)
       LOGICAL :: Create2ndKindBasis, PerformPiolaTransform, UsePretabulatedBasis, Parallel
       LOGICAL :: SecondOrder, ApplyTraceMapping
       INTEGER, POINTER :: EdgeMap(:,:), Ind(:)
       INTEGER :: TriangleFaceMap(3), SquareFaceMap(4), BrickFaceMap(6,4), PrismSquareFaceMap(3,4), DOFs
!----------------------------------------------------------------------------------------------------------
       Mesh => CurrentModel % Solver % Mesh
       !Parallel = ParEnv % PEs>1
       Parallel = ASSOCIATED(Mesh % ParallelInfo % Interface)

       stat = .TRUE.
       Basis = 0.0d0
       EdgeBasis = 0.0d0
       WorkBasis = 0.0d0
       CurlBasis = 0.0d0
       LG = 0.0d0
       !--------------------------------------------------------------------------------------------
       ! Check whether ready edge basis fuction values are available to reduce computation.
       ! If they are available, this function is used primarily to obtain the Piola transformation.
       !--------------------------------------------------------------------------------------------
       UsePretabulatedBasis = .FALSE.
       IF ( PRESENT(ReadyEdgeBasis) .AND. PRESENT(ReadyRotBasis) ) UsePretabulatedBasis = .TRUE.
       !------------------------------------------------------------------------------------------
       ! Check whether the Nedelec basis functions of the second kind or higher order basis
       ! functions should be created and whether the Piola transform is already applied within 
       ! this function.
       !------------------------------------------------------------------------------------------
       Create2ndKindBasis = .FALSE.
       IF ( PRESENT(SecondFamily) ) Create2ndKindBasis = SecondFamily
       SecondOrder = .FALSE.
       IF ( PRESENT(BasisDegree) ) THEN
         SecondOrder = BasisDegree > 1
       END IF
       PerformPiolaTransform = .FALSE.
       IF ( PRESENT(ApplyPiolaTransform) ) PerformPiolaTransform = ApplyPiolaTransform

       ApplyTraceMapping = .FALSE.
       IF ( PRESENT(TangentialTrMapping) ) ApplyTraceMapping = TangentialTrMapping
       !-------------------------------------------------------------------------------------------
       dLbasisdx = 0.0d0      
       n = Element % TYPE % NumberOfNodes
       dim = Element % TYPE % DIMENSION
       cdim = CoordinateSystemDimension()

       IF ( Element % TYPE % ElementCode == 101 ) THEN
         detF = 1.0d0
         Basis(1) = 1.0d0
         IF ( PRESENT(dBasisdx) ) dBasisdx(1,:) = 0.0d0
         RETURN
       END IF

       IF (cdim == 2 .AND. dim==1) THEN
         CALL Warn('EdgeElementInfo', 'Traces of 2-D edge elements have not been implemented yet')
         RETURN
       END IF

       !-----------------------------------------------------------------------
       ! The standard nodal basis functions on the reference element and
       ! their derivatives with respect to the local coordinates. These define 
       ! the mapping of the reference element to an actual element on the background 
       ! mesh but are not the basis functions for the edge element approximation.
       ! Remark: Using reference elements having the edges of the same length
       ! simplifies the implementation of element assembly procedures.
       !-----------------------------------------------------------------------
       SELECT CASE(Element % TYPE % ElementCode / 100)
       CASE(3)
         IF (SecondOrder) THEN
           ! DOFs is the number of H(curl)-conforming basis functions: 
           DOFs = 8
           IF (n == 6) THEN
             ! Here the element of the background mesh is of type 306.
             ! The Lagrange interpolation basis on the p-approximation reference element:
             Basis(1) = (3.0d0*u**2 + v*(-Sqrt(3.0d0) + v) + u*(-3.0d0 + 2.0d0*Sqrt(3.0d0)*v))/6.0d0
             dLBasisdx(1,1) = -0.5d0 + u + v/Sqrt(3.0d0)
             dLBasisdx(1,2) = (-Sqrt(3.0d0) + 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v)/6.0d0
             Basis(2) = (3.0d0*u**2 + v*(-Sqrt(3.0d0) + v) + u*(3.0d0 - 2.0d0*Sqrt(3.0d0)*v))/6.0d0
             dLBasisdx(2,1) = 0.5d0 + u - v/Sqrt(3.d0)
             dLBasisdx(2,2) = (-Sqrt(3.0d0) - 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v)/6.0d0
             Basis(3) = (v*(-Sqrt(3.0d0) + 2.0d0*v))/3.0d0
             dLBasisdx(3,1) = 0.0d0
             dLBasisdx(3,2) =  -(1.0d0/Sqrt(3.0d0)) + (4.0d0*v)/3.0d0
             Basis(4) = (3.0d0 - 3.0d0*u**2 - 2.0d0*Sqrt(3.0d0)*v + v**2)/3.0d0
             dLBasisdx(4,1) = -2.0d0*u
             dLBasisdx(4,2) = (-2.0d0*(Sqrt(3.0d0) - v))/3.0d0
             Basis(5) = (2.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - v)*v)/3.0d0
             dLBasisdx(5,1) =  (2.0d0*v)/Sqrt(3.0d0)
             dLBasisdx(5,2) = (2.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 2.0d0*v))/3.0d0
             Basis(6) = (-2.0d0*v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v))/3.0d0           
             dLBasisdx(6,1) = (-2.0d0*v)/Sqrt(3.0d0)
             dLBasisdx(6,2) = (-2.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 2.0d0*v))/3.0d0
           ELSE
             ! Here the element of the background mesh is of type 303:
             DO q=1,3
               Basis(q) = TriangleNodalPBasis(q, u, v)
               dLBasisdx(q,1:2) = dTriangleNodalPBasis(q, u, v) 
             END DO
           END IF
         ELSE
           DO q=1,n
             Basis(q) = TriangleNodalPBasis(q, u, v)
             dLBasisdx(q,1:2) = dTriangleNodalPBasis(q, u, v) 
           END DO
           IF (Create2ndKindBasis) THEN
             DOFs = 6
           ELSE
             DOFs = 3
           END IF
         END IF
       CASE(4)
         IF (SecondOrder) THEN
           ! The second-order quad from the Nedelec's first family: affine physical elements may be needed
           DOFs = 12
         ELSE
           ! The lowest-order quad from the optimal family (ABF_0)
           DOFs = 6
         END IF
         IF (n>4) THEN
           ! Here the background mesh is supposed to be of type 408/409
           CALL NodalBasisFunctions2D(Basis, Element, u, v)
           CALL NodalFirstDerivatives(n, dLBasisdx, Element, u, v, w)
         ELSE
           ! Here the background mesh is of type 404           
           DO q=1,4
             Basis(q) = QuadNodalPBasis(q, u, v)
             dLBasisdx(q,1:2) = dQuadNodalPBasis(q, u, v) 
           END DO
         END IF
       CASE(5)
         IF (SecondOrder) THEN
           DOFs = 20
           IF (n == 10) THEN
             ! Here the element of the background mesh is of type 510.
             ! The Lagrange interpolation basis on the p-approximation reference element:
             Basis(1) = (6.0d0*u**2 - 2.0d0*Sqrt(3.0d0)*v + 2.0d0*v**2 - Sqrt(6.0d0)*w + 2.0d0*Sqrt(2.0d0)*v*w + &
                 w**2 + 2.0d0*u*(-3.0d0 + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/12.0d0
             dLBasisdx(1,1) = -0.5d0 + u + v/Sqrt(3.0d0) + w/Sqrt(6.0d0)
             dLBasisdx(1,2) = (-Sqrt(3.0d0) + 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v + Sqrt(2.0d0)*w)/6.0d0
             dLBasisdx(1,3) = (-Sqrt(6.0d0) + 2.0d0*Sqrt(6.0d0)*u + 2.0d0*Sqrt(2.0d0)*v + 2.0d0*w)/12.0d0
             Basis(2) = (6.0d0*u**2 - 2.0d0*Sqrt(3.0d0)*v + 2.0d0*v**2 - Sqrt(6.0d0)*w + 2.0d0*Sqrt(2.0d0)*v*w + &
                 w**2 - 2.0d0*u*(-3.0d0 + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/12.0d0
             dLBasisdx(2,1) = 0.5d0 + u - v/Sqrt(3.0d0) - w/Sqrt(6.0d0)
             dLBasisdx(2,2) = (-Sqrt(3.0d0) - 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v + Sqrt(2.0d0)*w)/6.0d0
             dLBasisdx(2,3) = (-Sqrt(6.0d0) - 2.0d0*Sqrt(6.0d0)*u + 2.0d0*Sqrt(2.0d0)*v + 2.0d0*w)/12.0d0
             Basis(3) =  (8.0d0*v**2 + w*(Sqrt(6.0d0) + w) - 4.0d0*v*(Sqrt(3.0d0) + Sqrt(2.0d0)*w))/12.0d0
             dLBasisdx(3,1) = 0.0d0
             dLBasisdx(3,2) = (-Sqrt(3.0d0) + 4.0d0*v - Sqrt(2.0d0)*w)/3.0d0
             dLBasisdx(3,3) = (Sqrt(6.0d0) - 4.0d0*Sqrt(2.0d0)*v + 2.0d0*w)/12.0d0
             Basis(4) = (w*(-Sqrt(6.0d0) + 3.0d0*w))/4.0d0
             dLBasisdx(4,1) = 0.0d0
             dLBasisdx(4,2) = 0.0d0
             dLBasisdx(4,3) = (-Sqrt(6.0d0) + 6.0d0*w)/4.0d0
             Basis(5) =  (6.0d0 - 6.0d0*u**2 - 4.0d0*Sqrt(3.0d0)*v + 2.0d0*v**2 - 2.0d0*Sqrt(6.0d0)*w + &
                 2.0d0*Sqrt(2.0d0)*v*w + w**2)/6.0d0
             dLBasisdx(5,1) = -2.0d0*u
             dLBasisdx(5,2) = (-2.0d0*Sqrt(3.0d0) + 2.0d0*v + Sqrt(2.0d0)*w)/3.0d0
             dLBasisdx(5,3) = (-Sqrt(6.0d0) + Sqrt(2.0d0)*v + w)/3.0d0
             Basis(6) =  (-4.0d0*v**2 + w*(-Sqrt(6.0d0) - Sqrt(6.0d0)*u + w) + v*(4.0d0*Sqrt(3.0d0) + &
                 4.0d0*Sqrt(3.0d0)*u - Sqrt(2.0d0)*w))/6.0d0
             dLBasisdx(6,1) = (2.0d0*v)/Sqrt(3.0d0) - w/Sqrt(6.0d0)
             dLBasisdx(6,2) = (4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u - 8.0d0*v - Sqrt(2.0d0)*w)/6.0d0
             dLBasisdx(6,3) = (-Sqrt(6.0d0) - Sqrt(6.0d0)*u - Sqrt(2.0d0)*v + 2.0d0*w)/6.0d0
             Basis(7) =  (-4.0d0*v**2 + w*(-Sqrt(6.0d0) + Sqrt(6.0d0)*u + w) - &
                 v*(-4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u + Sqrt(2.0d0)*w))/6.0d0
             dLBasisdx(7,1) = (-2.0d0*v)/Sqrt(3.0d0) + w/Sqrt(6.0d0)
             dLBasisdx(7,2) = (4.0d0*Sqrt(3.0d0) - 4.0d0*Sqrt(3.0d0)*u - 8.0d0*v - Sqrt(2.0d0)*w)/6.0d0
             dLBasisdx(7,3) = (-Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v + 2.0d0*w)/6.0d0
             Basis(8) = -(w*(-Sqrt(6.0d0) + Sqrt(6.0d0)*u + Sqrt(2.0d0)*v + w))/2.0d0
             dLBasisdx(8,1) = -(Sqrt(1.5d0)*w)
             dLBasisdx(8,2) = -(w/Sqrt(2.0d0))
             dLBasisdx(8,3) = (Sqrt(6.0d0) - Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - 2.0d0*w)/2.0d0
             Basis(9) = ((Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - w)*w)/2.0d0
             dLBasisdx(9,1) = Sqrt(1.5d0)*w
             dLBasisdx(9,2) = -(w/Sqrt(2.0d0))
             dLBasisdx(9,3) = (Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - 2.0d0*w)/2.0d0
             Basis(10) = Sqrt(2.0d0)*v*w - w**2/2.0d0
             dLBasisdx(10,1) = 0.0d0
             dLBasisdx(10,2) = Sqrt(2.0d0)*w
             dLBasisdx(10,3) = Sqrt(2.0d0)*v - w
           ELSE
             ! Here the element of the background mesh is of type 504: 
             DO q=1,4
               Basis(q) = TetraNodalPBasis(q, u, v, w)
               dLBasisdx(q,1:3) = dTetraNodalPBasis(q, u, v, w)
             END DO
           END IF
         ELSE
           DO q=1,n
             Basis(q) = TetraNodalPBasis(q, u, v, w)
             dLBasisdx(q,1:3) = dTetraNodalPBasis(q, u, v, w)
           END DO
           IF (Create2ndKindBasis) THEN
             DOFs = 12
           ELSE
             DOFs = 6
           END IF
         END IF
       CASE(6)
         IF (SecondOrder) THEN
           ! The second-order pyramid from the Nedelec's first family
           DOFs = 31
         ELSE
           ! The lowest-order pyramid from the optimal family
           DOFs = 10
         END IF

         IF (n==13) THEN
           ! Here the background mesh is supposed to be of type 613. The difference between the standard
           ! reference element and the p-reference element can be taken into account by a simple scaling:
           CALL NodalBasisFunctions3D(Basis, Element, u, v, sqrt(2.0d0)*w)
           CALL NodalFirstDerivatives(n, dLBasisdx, Element, u, v, sqrt(2.0d0)*w)
           dLBasisdx(1:n,3) = sqrt(2.0d0) * dLBasisdx(1:n,3)
         ELSE
           ! Background mesh elements of the type 605:
           DO q=1,n
             Basis(q) = PyramidNodalPBasis(q, u, v, w)
             dLBasisdx(q,1:3) = dPyramidNodalPBasis(q, u, v, w)
           END DO
         END IF

       CASE(7)
         IF (SecondOrder) THEN
           ! The second-order prism from the Nedelec's first family: affine physical elements may be needed
           DOFs = 36
         ELSE
           ! The lowest-order prism from the optimal family
           DOFs = 15
         END IF

         IF (n==15) THEN
           ! Here the background mesh is of type 715.
           ! The Lagrange interpolation basis on the p-approximation reference element:

           h1 = -0.5d0*w + 0.5d0*w**2
           h2 = 0.5d0*w + 0.5d0*w**2
           h3 = 1.0d0 - w**2
           dh1 = -0.5d0 + w
           dh2 = 0.5d0 + w
           dh3 = -2.0d0 * w
           
           WorkBasis(1,1) = (3.0d0*u**2 + v*(-Sqrt(3.0d0) + v) + u*(-3.0d0 + 2.0d0*Sqrt(3.0d0)*v))/6
           grad(1) = -0.5d0 + u + v/Sqrt(3.0d0)
           grad(2) = (-Sqrt(3.0d0) + 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v)/6.0d0
           Basis(1) = WorkBasis(1,1) * h1
           dLBasisdx(1,1:2) = grad(1:2) * h1
           dLBasisdx(1,3) = WorkBasis(1,1) * dh1
           Basis(4) = WorkBasis(1,1) * h2
           dLBasisdx(4,1:2) = grad(1:2) * h2
           dLBasisdx(4,3) = WorkBasis(1,1) * dh2
           Basis(13) = WorkBasis(1,1) * h3
           dLBasisdx(13,1:2) = grad(1:2) * h3
           dLBasisdx(13,3) = WorkBasis(1,1) * dh3

           WorkBasis(1,1) = (3.0d0*u**2 + v*(-Sqrt(3.0d0) + v) + u*(3.0d0 - 2.0d0*Sqrt(3.0d0)*v))/6.0d0
           grad(1) = 0.5d0 + u - v/Sqrt(3.d0)
           grad(2) = (-Sqrt(3.0d0) - 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v)/6.0d0
           Basis(2) = WorkBasis(1,1) * h1
           dLBasisdx(2,1:2) = grad(1:2) * h1
           dLBasisdx(2,3) = WorkBasis(1,1) * dh1
           Basis(5) = WorkBasis(1,1) * h2
           dLBasisdx(5,1:2) = grad(1:2) * h2
           dLBasisdx(5,3) = WorkBasis(1,1) * dh2
           Basis(14) = WorkBasis(1,1) * h3
           dLBasisdx(14,1:2) = grad(1:2) * h3
           dLBasisdx(14,3) = WorkBasis(1,1) * dh3

           WorkBasis(1,1) = (v*(-Sqrt(3.0d0) + 2.0d0*v))/3.0d0
           grad(1) = 0.0d0
           grad(2) = -(1.0d0/Sqrt(3.0d0)) + (4.0d0*v)/3.0d0
           Basis(3) = WorkBasis(1,1) * h1
           dLBasisdx(3,1:2) = grad(1:2) * h1
           dLBasisdx(3,3) = WorkBasis(1,1) * dh1
           Basis(6) = WorkBasis(1,1) * h2
           dLBasisdx(6,1:2) = grad(1:2) * h2
           dLBasisdx(6,3) = WorkBasis(1,1) * dh2
           Basis(15) = WorkBasis(1,1) * h3
           dLBasisdx(15,1:2) = grad(1:2) * h3
           dLBasisdx(15,3) = WorkBasis(1,1) * dh3

           h1 = 0.5d0 * (1.0d0 - w)
           dh1 = -0.5d0
           h2 = 0.5d0 * (1.0d0 + w)
           dh2 = 0.5d0

           WorkBasis(1,1) = (3.0d0 - 3.0d0*u**2 - 2.0d0*Sqrt(3.0d0)*v + v**2)/3.0d0
           grad(1) = -2.0d0*u
           grad(2) = (-2.0d0*(Sqrt(3.0d0) - v))/3.0d0
           Basis(7) = WorkBasis(1,1) * h1
           dLBasisdx(7,1:2) = grad(1:2) * h1
           dLBasisdx(7,3) = WorkBasis(1,1) * dh1
           Basis(10) = WorkBasis(1,1) * h2
           dLBasisdx(10,1:2) = grad(1:2) * h2
           dLBasisdx(10,3) = WorkBasis(1,1) * dh2

           WorkBasis(1,1) = (2.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - v)*v)/3.0d0
           grad(1) = (2.0d0*v)/Sqrt(3.0d0)
           grad(2) = (2.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 2.0d0*v))/3.0d0
           Basis(8) = WorkBasis(1,1) * h1
           dLBasisdx(8,1:2) = grad(1:2) * h1
           dLBasisdx(8,3) = WorkBasis(1,1) * dh1
           Basis(11) = WorkBasis(1,1) * h2
           dLBasisdx(11,1:2) = grad(1:2) * h2
           dLBasisdx(11,3) = WorkBasis(1,1) * dh2

           WorkBasis(1,1) = (-2.0d0*v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v))/3.0d0
           grad(1) = (-2.0d0*v)/Sqrt(3.0d0)
           grad(2) = (-2.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 2.0d0*v))/3.0d0
           Basis(9) = WorkBasis(1,1) * h1
           dLBasisdx(9,1:2) = grad(1:2) * h1
           dLBasisdx(9,3) = WorkBasis(1,1) * dh1
           Basis(12) = WorkBasis(1,1) * h2
           dLBasisdx(12,1:2) = grad(1:2) * h2
           dLBasisdx(12,3) = WorkBasis(1,1) * dh2
         ELSE
           ! Here the background mesh is of type 706
           DO q=1,n
             Basis(q) = WedgeNodalPBasis(q, u, v, w)
             dLBasisdx(q,1:3) = dWedgeNodalPBasis(q, u, v, w)
           END DO
         END IF
       CASE(8)
         IF (SecondOrder) THEN
           ! The second-order brick from the Nedelec's first family: affine physical elements may be needed
           DOFs = 54
         ELSE
           ! The lowest-order brick from the optimal family
           DOFs = 27
         END IF
         IF (n>8) THEN
           ! Here the background mesh is supposed to be of type 820/827
           CALL NodalBasisFunctions3D(Basis, Element, u, v, w)
           CALL NodalFirstDerivatives(n, dLBasisdx, Element, u, v, w) 
         ELSE
           ! Here the background mesh is of type 808
           DO q=1,n
             Basis(q) = BrickNodalPBasis(q, u, v, w)
             dLBasisdx(q,1:3) = dBrickNodalPBasis(q, u, v, w)
           END DO
         END IF
       CASE DEFAULT
         CALL Fatal('ElementDescription::EdgeElementInfo','Unsupported element type')
       END SELECT

       !-----------------------------------------------------------------------
       ! Get data for performing the Piola transformation...
       !-----------------------------------------------------------------------
       stat = PiolaTransformationData(n, Element, Nodes, LF, detF, dLBasisdx) 
       !------------------------------------------------------------------------
       ! ... in order to define the basis for the element space X(K) via 
       ! applying a version of the Piola transformation as
       !    X(K) = { B | B = F^{-T}(f^{-1}(x)) b(f^{-1}(x)) }
       ! with b giving the edge basis function on the reference element k,
       ! f mapping k to the actual element K, i.e. K = f(k) and F = Grad f. This 
       ! function returns the local basis functions b and their Curl (with respect
       ! to local coordinates) evaluated at the integration point. The effect of 
       ! the Piola transformation need to be considered when integrating, so we 
       ! shall return also the values of F, G=F^{-T} and det F.
       !
       ! The construction of edge element bases could be done in an alternate way for 
       ! triangles and tetrahedra, while the chosen approach has the benefit that
       ! it generalizes to other cases. For example general quadrilaterals may now 
       ! be handled in the same way.
       !---------------------------------------------------------------------------
       IF (cdim == dim) THEN
          SELECT CASE(Element % TYPE % ElementCode / 100)
          CASE(3,4)
             LG(1,1) = 1.0d0/detF * LF(2,2)
             LG(1,2) = -1.0d0/detF * LF(1,2)
             LG(2,1) = -1.0d0/detF * LF(2,1)
             LG(2,2) = 1.0d0/detF * LF(1,1)
          CASE(5,6,7,8)
             CALL InvertMatrix3x3(LF,LG,detF)       
          CASE DEFAULT
             CALL Fatal('ElementDescription::EdgeElementInfo','Unsupported element type')
          END SELECT
          LG(1:dim,1:dim) = TRANSPOSE( LG(1:dim,1:dim) )
       END IF

       IF (UsePretabulatedBasis) THEN
         DO i=1,DOFs
           EdgeBasis(i,1:3) = ReadyEdgeBasis(i,1:3)
           CurlBasis(i,1:3) = ReadyRotBasis(i,1:3)
         END DO
       ELSE
         SELECT CASE(Element % TYPE % ElementCode / 100)
         CASE(3)
           !--------------------------------------------------------------
           ! This branch is for handling triangles. Note that
           ! the global orientation of the edge tanget t is defined such that
           ! t points towards the node that has a larger global index.
           !--------------------------------------------------------------
           EdgeMap => LGetEdgeMap(3)
           !EdgeMap => GetEdgeMap(GetElementFamily(Element))

           IF (Create2ndKindBasis) THEN
             !-------------------------------------------------
             ! Two basis functions defined on the edge 12.
             !-------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)             
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
               ! The sign and order of basis functions are reverted as
               ! compared with the other possibility
               EdgeBasis(1,1) = -(3.0d0 + 3.0d0*Sqrt(3.0d0)*u - Sqrt(3.0d0)*v)/6.0d0
               EdgeBasis(1,2) = -(3.0d0 + Sqrt(3.0d0)*u - Sqrt(3.0d0)*v)/6.0d0
               CurlBasis(1,3) = -1.0d0/Sqrt(3.0d0)

               EdgeBasis(2,1) = -(3.0d0 + Sqrt(3.0d0) - 3.0d0*(1.0d0 + Sqrt(3.0d0))*u - &
                   (1.0d0 + Sqrt(3.0d0))*v)/(2.0d0*(3.0d0 + Sqrt(3.0d0)))
               EdgeBasis(2,2) = -(-3.0d0 - Sqrt(3.0d0) + u + Sqrt(3.0d0)*u + v + Sqrt(3.0d0)*v)/ &
                   (2.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(2,3) = -1.0d0/Sqrt(3.0d0)
             ELSE
               EdgeBasis(1,1) = (3.0d0 + Sqrt(3.0d0) - 3.0d0*(1.0d0 + Sqrt(3.0d0))*u - &
                   (1.0d0 + Sqrt(3.0d0))*v)/(2.0d0*(3.0d0 + Sqrt(3.0d0)))
               EdgeBasis(1,2) = (-3.0d0 - Sqrt(3.0d0) + u + Sqrt(3.0d0)*u + v + Sqrt(3.0d0)*v)/ &
                   (2.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(1,3) = 1.0d0/Sqrt(3.0d0)

               EdgeBasis(2,1) = (3.0d0 + 3.0d0*Sqrt(3.0d0)*u - Sqrt(3.0d0)*v)/6.0d0
               EdgeBasis(2,2) = (3.0d0 + Sqrt(3.0d0)*u - Sqrt(3.0d0)*v)/6.0d0
               CurlBasis(2,3) = 1.0d0/Sqrt(3.0d0)                 
             END IF

             !-------------------------------------------------
             ! Two basis functions defined on the edge 23.
             !-------------------------------------------------
             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
               ! The sign and order of basis functions are reverted as
               ! compared with the other possibility
               EdgeBasis(3,1) = ((3.0d0 + Sqrt(3.0d0))*v)/6.0d0
               EdgeBasis(3,2) = -(-3.0d0 + Sqrt(3.0d0) + (-3.0d0 + Sqrt(3.0d0))*u + 2.0d0*Sqrt(3.0d0)*v)/6.0d0
               CurlBasis(3,3) = -1.0d0/Sqrt(3.0d0)

               EdgeBasis(4,1) = ((-3.0d0 + Sqrt(3.0d0))*v)/6.0d0
               EdgeBasis(4,2) = -(2.0d0 + Sqrt(3.0d0) + (2.0d0 + Sqrt(3.0d0))*u - &
                   (1.0d0 + Sqrt(3.0d0))*v)/(3.0d0 + Sqrt(3.0d0))
               CurlBasis(4,3) = -1.0d0/Sqrt(3.0d0)
             ELSE
               EdgeBasis(3,1) = -((-3.0d0 + Sqrt(3.0d0))*v)/6.0d0
               EdgeBasis(3,2) = (2.0d0 + Sqrt(3.0d0) + (2.0d0 + Sqrt(3.0d0))*u - &
                   (1.0d0 + Sqrt(3.0d0))*v)/(3.0d0 + Sqrt(3.0d0))
               CurlBasis(3,3) = 1.0d0/Sqrt(3.0d0)

               EdgeBasis(4,1) = -((3.0d0 + Sqrt(3.0d0))*v)/6.0d0
               EdgeBasis(4,2) = (-3.0d0 + Sqrt(3.0d0) + (-3.0d0 + Sqrt(3.0d0))*u + 2.0d0*Sqrt(3.0d0)*v)/6.0d0
               CurlBasis(4,3) = 1.0d0/Sqrt(3.0d0)                 
             END IF

             !-------------------------------------------------
             ! Two basis functions defined on the edge 31.
             !-------------------------------------------------
             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
               ! The sign and order of basis functions are reverted as
               ! compared with the other possibility
               EdgeBasis(5,1) = ((-3.0d0 + Sqrt(3.0d0))*v)/6.0d0
               EdgeBasis(5,2) = -(-3.0d0 - Sqrt(3.0d0) + (3.0d0 + Sqrt(3.0d0))*u + 2.0d0*Sqrt(3.0d0)*v)/6.0d0
               CurlBasis(5,3) = -1.0d0/Sqrt(3.0d0)

               EdgeBasis(6,1) = ((3.0d0 + 2.0d0*Sqrt(3.0d0))*v)/(3.0d0*(1.0d0 + Sqrt(3.0d0)))
               EdgeBasis(6,2) = ((-1.0d0 + u + v + Sqrt(3.0d0)*v)/(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(6,3) = -1.0d0/Sqrt(3.0d0)
             ELSE
               EdgeBasis(5,1) = -((3.0d0 + 2.0d0*Sqrt(3.0d0))*v)/(3.0d0*(1.0d0 + Sqrt(3.0d0)))
               EdgeBasis(5,2) = -((-1.0d0 + u + v + Sqrt(3.0d0)*v)/(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(5,3) = 1.0d0/Sqrt(3.0d0)

               EdgeBasis(6,1) = -((-3.0d0 + Sqrt(3.0d0))*v)/6.0d0
               EdgeBasis(6,2) = (-3.0d0 - Sqrt(3.0d0) + (3.0d0 + Sqrt(3.0d0))*u + 2.0d0*Sqrt(3.0d0)*v)/6.0d0
               CurlBasis(6,3) = 1.0d0/Sqrt(3.0d0)                 
             END IF

           ELSE
             
             !------------------------------------------------------------
             ! The optimal/Nedelec basis functions of the first kind. We employ
             ! a hierarchic basis, so the lowest-order basis functions are
             ! also utilized in the construction of the second-order basis. 
             ! First the edge 12 ...
             !------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(1,1) = (3.0d0 - Sqrt(3.0d0)*v)/6.0d0
             EdgeBasis(1,2) = u/(2.0d0*Sqrt(3.0d0))
             CurlBasis(1,3) = 1.0d0/Sqrt(3.0d0)
             IF (nj<ni) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,3) = -CurlBasis(1,3)
             END IF
             IF (SecondOrder) THEN
               EdgeBasis(2,1) = -(u*(-3.0d0 + Sqrt(3.0d0)*v))/2.0d0
               EdgeBasis(2,2) = (Sqrt(3.0d0)*u**2)/2.0d0
               CurlBasis(2,3) = (3.0d0*Sqrt(3.0d0)*u)/2.0d0                     
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the edge 23:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 3
               EdgeBasis(4,1) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*v)/4.0d0
               EdgeBasis(4,2) = (Sqrt(3.0d0)*(1.0d0 + u)*(-1.0d0 - u + Sqrt(3.0d0)*v))/4.0d0
               CurlBasis(4,3) = (-3.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v))/4.0d0
             ELSE
               k = 2
             END IF
             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(k,1) = -v/(2.0d0*Sqrt(3.0d0))
             EdgeBasis(k,2) = (1 + u)/(2.0d0*Sqrt(3.0d0))
             CurlBasis(k,3) =  1.0d0/Sqrt(3.0d0)
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,3) = -CurlBasis(k,3)
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the edge 31:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 5
               EdgeBasis(6,1) = (v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0
               EdgeBasis(6,2) = -(Sqrt(3.0d0)*(-1.0d0 + u)*(-1.0d0 + u + Sqrt(3.0d0)*v))/4.0d0
               CurlBasis(6,3) = (-3.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0                     
             ELSE
               k = 3
             END IF
             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(k,1) = -v/(2.0d0*Sqrt(3.0d0))
             EdgeBasis(k,2) = (-1 + u)/(2.0d0*Sqrt(3.0d0))
             CurlBasis(k,3) = 1.0d0/Sqrt(3.0d0)
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,3) = -CurlBasis(k,3)
             END IF

             IF (SecondOrder) THEN
               !-------------------------------------------------
               ! Two basis functions defined on the face 123:
               !-------------------------------------------------
               TriangleFaceMap(:) = (/ 1,2,3 /)          
               Ind => Element % Nodeindexes

               DO j=1,3
                 FaceIndeces(j) = Ind(TriangleFaceMap(j))
               END DO
               IF (Parallel) THEN
                 DO j=1,3
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                 END DO
               END IF
               CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

               WorkBasis(1,1) = ((Sqrt(3.0d0) - v)*v)/6.0d0
               WorkBasis(1,2) = (u*v)/6.0d0
               WorkCurlBasis(1,3) = (-Sqrt(3.0d0) + 3.0d0*v)/6.0d0
               WorkBasis(2,1) = (v*(1.0d0 + u - v/Sqrt(3.0d0)))/(4.0d0*Sqrt(3.0d0))
               WorkBasis(2,2) = ((-1.0d0 + u)*(-3.0d0 - 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
               WorkCurlBasis(2,3) =(-Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u + 3.0d0*v)/12.0d0
               WorkBasis(3,1) = (v*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
               WorkBasis(3,2) = -((1.0d0 + u)*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
               WorkCurlBasis(3,3) = (Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u - 3.0d0*v)/12.0d0

               EdgeBasis(7,:) = D1 * WorkBasis(I1,:)
               CurlBasis(7,3) = D1 * WorkCurlBasis(I1,3)
               EdgeBasis(8,:) = D2 * WorkBasis(I2,:)
               CurlBasis(8,3) = D2 * WorkCurlBasis(I2,3)  

             END IF
           END IF

         CASE(4)
           !--------------------------------------------------------------
           ! This branch is for handling quadrilaterals
           !--------------------------------------------------------------
           EdgeMap => LGetEdgeMap(4)
           IF (SecondOrder) THEN
             !---------------------------------------------------------------
             ! The second-order element from the Nedelec's first family with
             ! a hierarchic basis. This element may not be optimally accurate
             ! if the physical element is not affine.
             ! First, the eight basis functions associated with the edges:
             !--------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(1,1) = 0.1D1 / 0.4D1 - v / 0.4D1
             CurlBasis(1,3) = 0.1D1 / 0.4D1
             IF (nj<ni) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,3) = -CurlBasis(1,3)
             END IF
             EdgeBasis(2,1) = 0.3D1 * u * (0.1D1 / 0.4D1 - v / 0.4D1)
             CurlBasis(2,3) = 0.3D1 / 0.4D1 * u

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(3,2) = 0.1D1 / 0.4D1 + u / 0.4D1 
             CurlBasis(3,3) = 0.1D1 / 0.4D1
             IF (nj<ni) THEN
               EdgeBasis(3,:) = -EdgeBasis(3,:)
               CurlBasis(3,3) = -CurlBasis(3,3)
             END IF
             EdgeBasis(4,2) = 0.3D1 * v * (0.1D1 / 0.4D1 + u / 0.4D1)
             CurlBasis(4,3) = 0.3D1 / 0.4D1 * v

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(5,1) = -0.1D1 / 0.4D1 - v / 0.4D1
             CurlBasis(5,3) = 0.1D1 / 0.4D1
             IF (nj<ni) THEN
               EdgeBasis(5,:) = -EdgeBasis(5,:)
               CurlBasis(5,3) = -CurlBasis(5,3)
             END IF
             EdgeBasis(6,1) = -0.3D1 * u * (-0.1D1 / 0.4D1 - v / 0.4D1)
             CurlBasis(6,3) = -0.3D1 / 0.4D1 * u

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(7,2) = -0.1D1 / 0.4D1 + u / 0.4D1
             CurlBasis(7,3) = 0.1D1 / 0.4D1
             IF (nj<ni) THEN
               EdgeBasis(7,:) = -EdgeBasis(7,:)
               CurlBasis(7,3) = -CurlBasis(7,3)
             END IF
             EdgeBasis(8,2) = -0.3D1 * v * (-0.1D1 / 0.4D1 + u / 0.4D1)
             CurlBasis(8,3) = -0.3D1 / 0.4D1 * v

             !--------------------------------------------------------------------
             ! Additional four basis functions associated with the element interior
             !-------------------------------------------------------------------
             SquareFaceMap(:) = (/ 1,2,3,4 /)          
             Ind => Element % Nodeindexes

             WorkBasis = 0.0d0
             WorkCurlBasis = 0.0d0

             WorkBasis(1,1) = 0.2D1 * (0.1D1 / 0.2D1 - v / 0.2D1) * (0.1D1 / 0.2D1 + v / 0.2D1)
             WorkCurlBasis(1,3) = v
             WorkBasis(2,1) = 0.12D2 * u * (0.1D1 / 0.2D1 - v / 0.2D1) * (0.1D1 / 0.2D1 + v / 0.2D1)
             WorkCurlBasis(2,3) = 0.6D1 * u * (0.1D1 / 0.2D1 + v / 0.2D1) - &
                 0.6D1 * u * (0.1D1 / 0.2D1 - v / 0.2D1)

             WorkBasis(3,2) = 0.2D1 * (0.1D1 / 0.2D1 - u / 0.2D1) * (0.1D1 / 0.2D1 + u / 0.2D1)
             WorkCurlBasis(3,3) = -u
             WorkBasis(4,2) = 0.12D2 * v * (0.1D1 / 0.2D1 - u / 0.2D1) * (0.1D1 / 0.2D1 + u / 0.2D1)
             WorkCurlBasis(4,3) = -0.6D1 * v * (0.1D1 / 0.2D1 + u / 0.2D1) + &
                 0.6D1 * v * (0.1D1 / 0.2D1 - u / 0.2D1)

             DO j=1,4
               FaceIndeces(j) = Ind(SquareFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(9,:) = D1 * WorkBasis(2*(I1-1)+1,:)
             CurlBasis(9,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
             EdgeBasis(10,:) = WorkBasis(2*(I1-1)+2,:)
             CurlBasis(10,:) = WorkCurlBasis(2*(I1-1)+2,:)
             EdgeBasis(11,:) = D2 * WorkBasis(2*(I2-1)+1,:)
             CurlBasis(11,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
             EdgeBasis(12,:) = WorkBasis(2*(I2-1)+2,:)
             CurlBasis(12,:) = WorkCurlBasis(2*(I2-1)+2,:)

           ELSE
             !------------------------------------------------------
             ! The Arnold-Boffi-Falk element of degree k=0 which is
             ! a member of the optimal edge element family. 
             ! First, four basis functions defined on the edges
             !-------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(1,1) = ((-1.0d0 + v)*v)/4.0d0
             EdgeBasis(1,2) = 0.0d0
             CurlBasis(1,3) = (1.0d0 - 2*v)/4.0d0
             IF (nj<ni) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,3) = -CurlBasis(1,3)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(2,1) = 0.0d0
             EdgeBasis(2,2) = (u*(1.0d0 + u))/4.0d0
             CurlBasis(2,3) = (1.0d0 + 2*u)/4.0d0
             IF (nj<ni) THEN
               EdgeBasis(2,:) = -EdgeBasis(2,:)
               CurlBasis(2,3) = -CurlBasis(2,3)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(3,1) = -(v*(1.0d0 + v))/4.0d0
             EdgeBasis(3,2) = 0.0d0
             CurlBasis(3,3) = (1.0d0 + 2*v)/4.0d0
             IF (nj<ni) THEN
               EdgeBasis(3,:) = -EdgeBasis(3,:)
               CurlBasis(3,3) = -CurlBasis(3,3)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(4,1) = 0.0d0
             EdgeBasis(4,2) = -((-1 + u)*u)/4.0d0
             CurlBasis(4,3) = (1.0d0 - 2*u)/4.0d0
             IF (nj<ni) THEN
               EdgeBasis(4,:) = -EdgeBasis(4,:)
               CurlBasis(4,3) = -CurlBasis(4,3)
             END IF

             !--------------------------------------------------------------------
             ! Additional two basis functions associated with the element interior
             !-------------------------------------------------------------------
             SquareFaceMap(:) = (/ 1,2,3,4 /)          
             Ind => Element % Nodeindexes

             WorkBasis(1,:) = 0.0d0
             WorkBasis(2,:) = 0.0d0
             WorkCurlBasis(1,:) = 0.0d0
             WorkCurlBasis(2,:) = 0.0d0         

             WorkBasis(1,1) = (1.0d0 - v**2)/2.0d0
             WorkBasis(1,2) = 0.0d0
             WorkCurlBasis(1,3) = v

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = (1.0d0 - u**2)/2.0d0
             WorkCurlBasis(2,3) = -u

             DO j=1,4
               FaceIndeces(j) = Ind(SquareFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(5,:) = D1 * WorkBasis(I1,:)
             CurlBasis(5,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(6,:) = D2 * WorkBasis(I2,:)
             CurlBasis(6,:) = D2 * WorkCurlBasis(I2,:)         
           END IF

         CASE(5)
           !--------------------------------------------------------------
           ! This branch is for handling tetrahedra
           !--------------------------------------------------------------
           EdgeMap => LGetEdgeMap(5)

           IF (Create2ndKindBasis) THEN
             !-------------------------------------------------
             ! Two basis functions defined on the edge 12.
             !-------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
               ! The sign and order of basis functions are reverted as
               ! compared with the other possibility
               EdgeBasis(1,1) = -(6.0d0 + 6.0d0*Sqrt(3.0d0)*u - 2.0d0*Sqrt(3.0d0)*v - Sqrt(6.0d0)*w)/12.0d0
               EdgeBasis(1,2) = -(6.0d0 + 2.0d0*Sqrt(3.0d0)*u - 2.0d0*Sqrt(3.0d0)*v - Sqrt(6.0d0)*w)/12.0d0
               EdgeBasis(1,3) = -(3.0d0*Sqrt(2.0d0) + Sqrt(6.0d0)*u - Sqrt(6.0d0)*v - Sqrt(3.0d0)*w)/12.0d0
               CurlBasis(1,1) = 0.0d0
               CurlBasis(1,2) = 1.0d0/Sqrt(6.0d0)
               CurlBasis(1,3) = -1.0d0/Sqrt(3.0d0)

               EdgeBasis(2,1) = (-6.0d0 - 2.0d0*Sqrt(3.0d0) + 6.0d0*(1.0d0 + Sqrt(3.0d0))*u + &
                   2.0d0*(1.0d0 + Sqrt(3.0d0))*v + Sqrt(2.0d0)*w + Sqrt(6.0d0)*w)/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               EdgeBasis(2,2) = -(-6.0d0 - 2.0d0*Sqrt(3.0d0) + 2.0d0*(1.0d0 + Sqrt(3.0d0))*u + &
                   2.0d0*(1.0d0 + Sqrt(3.0d0))*v + Sqrt(2.0d0)*w + Sqrt(6.0d0)*w)/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               EdgeBasis(2,3) = -(-3.0d0*Sqrt(2.0d0) - Sqrt(6.0d0) + (Sqrt(2.0d0) + Sqrt(6.0d0))*u + &
                   (Sqrt(2.0d0) + Sqrt(6.0d0))*v + w + Sqrt(3.0d0)*w)/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(2,1) = 0.0d0 
               CurlBasis(2,2) = (Sqrt(2.0d0) + Sqrt(6.0d0))/(6.0d0 + 2.0d0*Sqrt(3.0d0))
               CurlBasis(2,3) = -1.0d0/Sqrt(3.0d0)
             ELSE
               EdgeBasis(1,1) = -(-6.0d0 - 2.0d0*Sqrt(3.0d0) + 6.0d0*(1.0d0 + Sqrt(3.0d0))*u + &
                   2.0d0*(1.0d0 + Sqrt(3.0d0))*v + Sqrt(2.0d0)*w + Sqrt(6.0d0)*w)/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               EdgeBasis(1,2) = (-6.0d0 - 2.0d0*Sqrt(3.0d0) + 2.0d0*(1.0d0 + Sqrt(3.0d0))*u + &
                   2.0d0*(1.0d0 + Sqrt(3.0d0))*v + Sqrt(2.0d0)*w + Sqrt(6.0d0)*w)/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               EdgeBasis(1,3) = (-3.0d0*Sqrt(2.0d0) - Sqrt(6.0d0) + (Sqrt(2.0d0) + Sqrt(6.0d0))*u + &
                   (Sqrt(2.0d0) + Sqrt(6.0d0))*v + w + Sqrt(3.0d0)*w)/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(1,1) = 0.0d0
               CurlBasis(1,2) = -((Sqrt(2.0d0) + Sqrt(6.0d0))/(6.0d0 + 2.0d0*Sqrt(3.0d0)))
               CurlBasis(1,3) = 1.0d0/Sqrt(3.0d0)

               EdgeBasis(2,1) = (6.0d0 + 6.0d0*Sqrt(3.0d0)*u - 2.0d0*Sqrt(3.0d0)*v - Sqrt(6.0d0)*w)/12.0d0
               EdgeBasis(2,2) = (6.0d0 + 2.0d0*Sqrt(3.0d0)*u - 2.0d0*Sqrt(3.0d0)*v - Sqrt(6.0d0)*w)/12.0d0
               EdgeBasis(2,3) = (3.0d0*Sqrt(2.0d0) + Sqrt(6.0d0)*u - Sqrt(6.0d0)*v - Sqrt(3.0d0)*w)/12.0d0 
               CurlBasis(2,1) = 0.0d0
               CurlBasis(2,2) = -1.0d0/Sqrt(6.0d0)
               CurlBasis(2,3) = 1.0d0/Sqrt(3.0d0)
             END IF

             !-------------------------------------------------
             ! Two basis functions defined on the edge 23.
             !-------------------------------------------------
             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
               ! The sign and order of basis functions are reverted as
               ! compared with the other possibility
               EdgeBasis(3,1) = (3.0d0 + Sqrt(3.0d0))*(4.0d0*v - Sqrt(2.0d0)*w)/24.0d0
               EdgeBasis(3,2) = -(4.0d0*(-3.0d0 + Sqrt(3.0d0))*u + 8.0d0*Sqrt(3.0d0)*v + &
                   (-3.0d0 + Sqrt(3.0d0))*(4.0d0 + Sqrt(2.0d0)*w))/24.0d0
               EdgeBasis(3,3) = -(3.0d0*Sqrt(2.0d0) - Sqrt(6.0d0) - Sqrt(2.0d0)*(-3.0d0 + Sqrt(3.0d0))*u + &
                   Sqrt(2.0d0)*(3.0d0 + Sqrt(3.0d0))*v - 2.0d0*Sqrt(3.0d0)*w)/24.0d0
               CurlBasis(3,1) = -1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(3,2) = -1.0d0/(2.0d0*Sqrt(6.0d0))
               CurlBasis(3,3) = -1.0d0/Sqrt(3.0d0)

               EdgeBasis(4,1) = (-3.0d0 + Sqrt(3.0d0))*(4.0d0*v - Sqrt(2.0d0)*w)/24.0d0
               EdgeBasis(4,2) = (-4.0d0*(2.0d0 + Sqrt(3.0d0))*u + 4.0d0*(1.0d0 + Sqrt(3.0d0))*v + &
                   (2.0d0 + Sqrt(3.0d0))*(-4.0d0 + Sqrt(2.0d0)*w))/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               EdgeBasis(4,3) = -(-2.0d0*Sqrt(2.0d0) - Sqrt(6.0d0) - Sqrt(2.0d0)*(2.0d0 + Sqrt(3.0d0))*u + &
                   Sqrt(2.0d0)*v + w + Sqrt(3.0d0)*w)/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(4,1) = -1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(4,2) = -(Sqrt(2.0d0) + Sqrt(6.0d0))/(12.0d0 + 4.0d0*Sqrt(3.0d0))
               CurlBasis(4,3) = -1.0d0/Sqrt(3.0d0)
             ELSE
               EdgeBasis(3,1) = -(-3.0d0 + Sqrt(3.0d0))*(4.0d0*v - Sqrt(2.0d0)*w)/24.0d0
               EdgeBasis(3,2) = -(-4.0d0*(2.0d0 + Sqrt(3.0d0))*u + 4.0d0*(1.0d0 + Sqrt(3.0d0))*v + &
                   (2.0d0 + Sqrt(3.0d0))*(-4.0d0 + Sqrt(2.0d0)*w))/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               EdgeBasis(3,3) = (-2.0d0*Sqrt(2.0d0) - Sqrt(6.0d0) - Sqrt(2.0d0)*(2.0d0 + Sqrt(3.0d0))*u + &
                   Sqrt(2.0d0)*v + w + Sqrt(3.0d0)*w)/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(3,1) = 1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(3,2) = (Sqrt(2.0d0) + Sqrt(6.0d0))/(12.0d0 + 4.0d0*Sqrt(3.0d0))
               CurlBasis(3,3) = 1.0d0/Sqrt(3.0d0)

               EdgeBasis(4,1) = -((3.0d0 + Sqrt(3.0d0))*(4.0d0*v - Sqrt(2.0d0)*w))/24.0d0
               EdgeBasis(4,2) = (4.0d0*(-3.0d0 + Sqrt(3.0d0))*u + 8.0d0*Sqrt(3.0d0)*v + &
                   (-3.0d0 + Sqrt(3.0d0))*(4.0d0 + Sqrt(2.0d0)*w))/24.0d0
               EdgeBasis(4,3) = (3.0d0*Sqrt(2.0d0) - Sqrt(6.0d0) - Sqrt(2.0d0)*(-3.0d0 + Sqrt(3.0d0))*u + &
                   Sqrt(2.0d0)*(3.0d0 + Sqrt(3.0d0))*v - 2.0d0*Sqrt(3.0d0)*w)/24.0d0
               CurlBasis(4,1) = 1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(4,2) = 1.0d0/(2.0d0*Sqrt(6.0d0))
               CurlBasis(4,3) = 1.0d0/Sqrt(3.0d0)
             END IF

             !-------------------------------------------------
             ! Two basis functions defined on the edge 31.
             !-------------------------------------------------
             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
               ! The sign and order of basis functions are reverted as
               ! compared with the other possibility
               EdgeBasis(5,1) = ((-3.0d0 + Sqrt(3.0d0))*(4.0d0*v - Sqrt(2.0d0)*w))/24.0d0
               EdgeBasis(5,2) = -(4.0d0*(3.0d0 + Sqrt(3.0d0))*u + 8.0d0*Sqrt(3.0d0)*v + &
                   (3.0d0 + Sqrt(3.0d0))*(-4.0d0 + Sqrt(2.0d0)*w))/24.0d0
               EdgeBasis(5,3) = -(3.0d0*Sqrt(2.0d0) + Sqrt(6.0d0) - Sqrt(2.0d0)*(3.0d0 + Sqrt(3.0d0))*u + &
                   Sqrt(2.0d0)*(-3.0d0 + Sqrt(3.0d0))*v - 2.0d0*Sqrt(3.0d0)*w)/24.0d0
               CurlBasis(5,1) = 1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(5,2) = -1.0d0/(2.0d0*Sqrt(6.0d0))
               CurlBasis(5,3) = -1.0d0/Sqrt(3.0d0)

               EdgeBasis(6,1) = ((3.0d0 + 2.0d0*Sqrt(3.0d0))*(4.0d0*v - Sqrt(2.0d0)*w))/ &
                   (12.0d0*(1.0d0 + Sqrt(3.0d0)))
               EdgeBasis(6,2) = -(4.0d0 - 4.0d0*u - 4.0d0*(1.0d0 + Sqrt(3.0d0))*v + Sqrt(2.0d0)*w)/ &
                   (4.0d0*(3.0d0 + Sqrt(3.0d0)))
               EdgeBasis(6,3) = -(-Sqrt(2.0d0) + Sqrt(2.0d0)*u - Sqrt(2.0d0)*(2.0d0 + Sqrt(3.0d0))*v + &
                   w + Sqrt(3.0d0)*w)/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(6,1) = 1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(6,2) = -1.0d0/(2.0d0*Sqrt(6.0d0))
               CurlBasis(6,3) = -1.0d0/Sqrt(3.0d0)
             ELSE
               EdgeBasis(5,1) = -((3.0d0 + 2.0d0*Sqrt(3.0d0))*(4.0d0*v - Sqrt(2.0d0)*w))/ &
                   (12.0d0*(1.0d0 + Sqrt(3.0d0)))
               EdgeBasis(5,2) = (4.0d0 - 4.0d0*u - 4.0d0*(1.0d0 + Sqrt(3.0d0))*v + Sqrt(2.0d0)*w)/ &
                   (4.0d0*(3.0d0 + Sqrt(3.0d0)))
               EdgeBasis(5,3) = (-Sqrt(2.0d0) + Sqrt(2.0d0)*u - Sqrt(2.0d0)*(2.0d0 + Sqrt(3.0d0))*v + &
                   w + Sqrt(3.0d0)*w)/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(5,1) = -1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(5,2) = 1.0d0/(2.0d0*Sqrt(6.0d0))
               CurlBasis(5,3) = 1.0d0/Sqrt(3.0d0)

               EdgeBasis(6,1) = -((-3.0d0 + Sqrt(3.0d0))*(4.0d0*v - Sqrt(2.0d0)*w))/24.0d0
               EdgeBasis(6,2) = (4.0d0*(3.0d0 + Sqrt(3.0d0))*u + 8.0d0*Sqrt(3.0d0)*v + &
                   (3.0d0 + Sqrt(3.0d0))*(-4.0d0 + Sqrt(2.0d0)*w))/24.0d0
               EdgeBasis(6,3) = (3.0d0*Sqrt(2.0d0) + Sqrt(6.0d0) - Sqrt(2.0d0)*(3.0d0 + Sqrt(3.0d0))*u + &
                   Sqrt(2.0d0)*(-3.0d0 + Sqrt(3.0d0))*v - 2.0d0*Sqrt(3.0d0)*w)/24.0d0
               CurlBasis(6,1) = -1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(6,2) = 1.0d0/(2.0d0*Sqrt(6.0d0))
               CurlBasis(6,3) = 1.0d0/Sqrt(3.0d0)
             END IF

             !-------------------------------------------------
             ! Two basis functions defined on the edge 14.
             !-------------------------------------------------
             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
               ! The sign and order of basis functions are reverted as
               ! compared with the other possibility
               EdgeBasis(7,1) = -((3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(2.0d0))
               EdgeBasis(7,2) = -((3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(6.0d0))
               EdgeBasis(7,3) = -(-3.0d0*Sqrt(2.0d0) + Sqrt(6.0d0) - Sqrt(2.0d0)*(-3.0d0 + Sqrt(3.0d0))*u + &
                   Sqrt(2.0d0)*(-1.0d0 + Sqrt(3.0d0))*v + 2.0d0*Sqrt(3.0d0)*w)/8.0d0
               CurlBasis(7,1) = 1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(7,2) = -Sqrt(1.5d0)/2.0d0
               CurlBasis(7,3) = 0.0d0

               EdgeBasis(8,1) = -((-3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(2.0d0))
               EdgeBasis(8,2) = -((-3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(6.0d0))
               EdgeBasis(8,3) = -((-3.0d0 + Sqrt(3.0d0))*w - (Sqrt(2.0d0)*(3.0d0 + 2.0d0*Sqrt(3.0d0))* &
                   (-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(3.0d0 + Sqrt(3.0d0)))/ &
                   (8.0d0*Sqrt(3.0d0))
               CurlBasis(8,1) = 1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(8,2) = -(3.0d0*(Sqrt(2.0d0) + Sqrt(6.0d0)))/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(8,3) = 0.0d0 
             ELSE
               EdgeBasis(7,1) = ((-3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(2.0d0))
               EdgeBasis(7,2) = ((-3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(6.0d0))
               EdgeBasis(7,3) = ((-3.0d0 + Sqrt(3.0d0))*w - (Sqrt(2.0d0)*(3.0d0 + 2.0d0*Sqrt(3.0d0))* &
                   (-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(3.0d0 + Sqrt(3.0d0)))/ &
                   (8.0d0*Sqrt(3.0d0))
               CurlBasis(7,1) = -1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(7,2) = (3.0d0*(Sqrt(2.0d0) + Sqrt(6.0d0)))/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(7,3) = 0.0d0

               EdgeBasis(8,1) = ((3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(2.0d0))
               EdgeBasis(8,2) = ((3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(6.0d0))
               EdgeBasis(8,3) = (-3.0d0*Sqrt(2.0d0) + Sqrt(6.0d0) - Sqrt(2.0d0)*(-3.0d0 + Sqrt(3.0d0))*u + &
                   Sqrt(2.0d0)*(-1.0d0 + Sqrt(3.0d0))*v + 2.0d0*Sqrt(3.0d0)*w)/8.0d0
               CurlBasis(8,1) = -1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(8,2) = Sqrt(1.5d0)/2.0d0
               CurlBasis(8,3) = 0.0d0
             END IF

             !-------------------------------------------------
             ! Two basis functions defined on the edge 24.
             !-------------------------------------------------
             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
               ! The sign and order of basis functions are reverted as
               ! compared with the other possibility
               EdgeBasis(9,1) = ((3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(2.0d0))
               EdgeBasis(9,2) = -((3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(6.0d0))
               EdgeBasis(9,3) = -(-3.0d0*Sqrt(2.0d0) + Sqrt(6.0d0) + Sqrt(2.0d0)*(-3.0d0 + Sqrt(3.0d0))*u + &
                   Sqrt(2.0d0)*(-1.0d0 + Sqrt(3.0d0))*v + 2.0d0*Sqrt(3.0d0)*w)/8.0d0
               CurlBasis(9,1) = 1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(9,2) = Sqrt(1.5d0)/2.0d0
               CurlBasis(9,3) = 0.0d0

               EdgeBasis(10,1) = ((-3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(2.0d0))
               EdgeBasis(10,2) = -((-3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(6.0d0))
               EdgeBasis(10,3) = -((-3.0d0 + Sqrt(3.0d0))*w - (Sqrt(2.0d0)*(3.0d0 + 2.0d0*Sqrt(3.0d0))*&
                   (-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/&
                   (3.0d0 + Sqrt(3.0d0)))/(8.0d0*Sqrt(3.0d0))
               CurlBasis(10,1) = 1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(10,2) = -(-3.0d0*(Sqrt(2.0d0) + Sqrt(6.0d0)))/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(10,3) = 0.0d0
             ELSE
               EdgeBasis(9,1) = -((-3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(2.0d0))
               EdgeBasis(9,2) = ((-3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(6.0d0))
               EdgeBasis(9,3) = ((-3.0d0 + Sqrt(3.0d0))*w - (Sqrt(2.0d0)*(3.0d0 + 2.0d0*Sqrt(3.0d0))*&
                   (-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/&
                   (3.0d0 + Sqrt(3.0d0)))/(8.0d0*Sqrt(3.0d0))
               CurlBasis(9,1) = -1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(9,2) = (-3.0d0*(Sqrt(2.0d0) + Sqrt(6.0d0)))/(4.0d0*(3.0d0 + Sqrt(3.0d0)))
               CurlBasis(9,3) = 0.0d0

               EdgeBasis(10,1) = -((3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(2.0d0))
               EdgeBasis(10,2) = ((3.0d0 + Sqrt(3.0d0))*w)/(4.0d0*Sqrt(6.0d0))
               EdgeBasis(10,3) = (-3.0d0*Sqrt(2.0d0) + Sqrt(6.0d0) + Sqrt(2.0d0)*(-3.0d0 + Sqrt(3.0d0))*u + &
                   Sqrt(2.0d0)*(-1.0d0 + Sqrt(3.0d0))*v + 2.0d0*Sqrt(3.0d0)*w)/8.0d0
               CurlBasis(10,1) = -1.0d0/(2.0d0*Sqrt(2.0d0))
               CurlBasis(10,2) = -Sqrt(1.5d0)/2.0d0
               CurlBasis(10,3) = 0.0d0
             END IF

             !-------------------------------------------------
             ! Two basis functions defined on the edge 34.
             !-------------------------------------------------
             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) THEN
               ! The sign and order of basis functions are reverted as
               ! compared with the other possibility
               EdgeBasis(11,1) = 0.0d0
               EdgeBasis(11,2) = ((1.0d0 + Sqrt(3.0d0))*w)/(2.0d0*Sqrt(2.0d0))
               EdgeBasis(11,3) = -(6.0d0*Sqrt(2.0d0)*v - 4.0d0*Sqrt(6.0d0)*v - 3.0d0*w + &
                   3.0d0*Sqrt(3.0d0)*w)/(12.0d0 - 4.0d0*Sqrt(3.0d0))
               CurlBasis(11,1) = -1.0d0/Sqrt(2.0d0)
               CurlBasis(11,2) = 0.0d0
               CurlBasis(11,3) = 0.0d0

               EdgeBasis(12,1) = 0.0d0
               EdgeBasis(12,2) = ((-3.0d0 + Sqrt(3.0d0))*w)/(2.0d0*Sqrt(6.0d0))
               EdgeBasis(12,3) = -((Sqrt(2.0d0) + Sqrt(6.0d0))*v - Sqrt(3.0d0)*w)/4.0d0
               CurlBasis(12,1) = -1.0d0/Sqrt(2.0d0)
               CurlBasis(12,2) = 0.0d0
               CurlBasis(12,3) = 0.0d0
             ELSE
               EdgeBasis(11,1) = 0.0d0
               EdgeBasis(11,2) = -((-3.0d0 + Sqrt(3.0d0))*w)/(2.0d0*Sqrt(6.0d0))
               EdgeBasis(11,3) = ((Sqrt(2.0d0) + Sqrt(6.0d0))*v - Sqrt(3.0d0)*w)/4.0d0
               CurlBasis(11,1) = 1.0d0/Sqrt(2.0d0)
               CurlBasis(11,2) = 0.0d0
               CurlBasis(11,3) = 0.0d0

               EdgeBasis(12,1) = 0.0d0
               EdgeBasis(12,2) = -((1.0d0 + Sqrt(3.0d0))*w)/(2.0d0*Sqrt(2.0d0))
               EdgeBasis(12,3) = (6.0d0*Sqrt(2.0d0)*v - 4.0d0*Sqrt(6.0d0)*v - 3.0d0*w + &
                   3.0d0*Sqrt(3.0d0)*w)/(12.0d0 - 4.0d0*Sqrt(3.0d0))
               CurlBasis(12,1) = 1.0d0/Sqrt(2.0d0)
               CurlBasis(12,2) = 0.0d0
               CurlBasis(12,3) = 0.0d0
             END IF

           ELSE
             
             !-------------------------------------------------------------
             ! The optimal/Nedelec basis functions of the first kind. We employ
             ! a hierarchic basis, so the lowest-order basis functions are
             ! also utilized in the construction of the second-order basis. 
             ! The first the edge ...
             !-------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(1,1) = (6.0d0 - 2.0d0*Sqrt(3.0d0)*v - Sqrt(6.0d0)*w)/24.0d0
             EdgeBasis(1,2) = u/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(1,3) = u/(4.0d0*Sqrt(6.0d0))            
             CurlBasis(1,1) = 0.0d0
             CurlBasis(1,2) = -1.0d0/(2.0d0*Sqrt(6.0d0))
             CurlBasis(1,3) = 1.0d0/(2.0d0*Sqrt(3.0d0))
             IF (nj<ni) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,:) = -CurlBasis(1,:)
             END IF
             IF (SecondOrder) THEN
               EdgeBasis(2,1) = -(u*(-6.0d0 + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/4.0d0
               EdgeBasis(2,2) = (Sqrt(3.0d0)*u**2)/2.0d0
               EdgeBasis(2,3) = (Sqrt(1.5d0)*u**2)/2.0d0
               CurlBasis(2,1) = 0.0d0
               CurlBasis(2,2) = (-3.0d0*Sqrt(1.5d0)*u)/2.0d0
               CurlBasis(2,3) = (3.0d0*Sqrt(3.0d0)*u)/2.0d0                   
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the second edge:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 3
               EdgeBasis(4,1) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*(4.0d0*v - Sqrt(2.0d0)*w))/16.0d0
               EdgeBasis(4,2) = -((1.0d0 + u - Sqrt(3.0d0)*v)*&
                   (4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u - 3.0d0*Sqrt(2.0d0)*w))/16.0d0
               EdgeBasis(4,3) = -((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*&
                   (-1.0d0 - u + Sqrt(3.0d0)*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(4,1) = (-9.0d0*(1.0d0 + u - Sqrt(3.0d0)*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(4,2) = (-3.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(4,3) = (-3.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v))/4.0d0
             ELSE
               k = 2
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(k,1) = (-4.0d0*v + Sqrt(2.0d0)*w)/(16.0d0*Sqrt(3.0d0))
             EdgeBasis(k,2) = (4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u - 3.0d0*Sqrt(2.0d0)*w)/48.0d0
             EdgeBasis(k,3) = -(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)/(24.0d0*Sqrt(2.0d0))
             CurlBasis(k,1) = 1.0d0/(4.0d0*Sqrt(2.0d0))
             CurlBasis(k,2) = 1.0d0/(4.0d0*Sqrt(6.0d0))
             CurlBasis(k,3) = 1.0d0/(2.0d0*Sqrt(3.0d0))
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the third edge:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 5
               EdgeBasis(6,1) = ((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v)*&
                   (4.0d0*v - Sqrt(2.0d0)*w))/16.0d0
               EdgeBasis(6,2) = -((-1.0d0 + u + Sqrt(3.0d0)*v)*&
                   (-4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u + 3.0d0*Sqrt(2.0d0)*w))/16.0d0
               EdgeBasis(6,3) = ((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v)*&
                   (-1.0d0 + u + Sqrt(3.0d0)*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(6,1) = (9.0d0*(-1.0d0 + u + Sqrt(3.0d0)*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(6,2) = (-3.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(6,3) = (-3.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0
             ELSE
               k = 3
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(k,1) = (-4.0d0*v + Sqrt(2.0d0)*w)/(16.0d0*Sqrt(3.0d0))
             EdgeBasis(k,2) = (-4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u + 3.0d0*Sqrt(2.0d0)*w)/48.0d0
             EdgeBasis(k,3) = (Sqrt(6.0d0) - Sqrt(6.0d0)*u - 3.0d0*Sqrt(2.0d0)*v)/48.0d0
             CurlBasis(k,1) = -1.0d0/(4.0d0*Sqrt(2.0d0))
             CurlBasis(k,2) = 1.0d0/(4.0d0*Sqrt(6.0d0))
             CurlBasis(k,3) = 1.0d0/(2.0d0*Sqrt(3.0d0))
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the fourth edge:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 7
               EdgeBasis(8,1) = (3.0d0*w*(-Sqrt(6.0d0) + Sqrt(6.0d0)*u + Sqrt(2.0d0)*v + 4.0d0*w))/16.0d0
               EdgeBasis(8,2) = (w*(-3.0d0*Sqrt(2.0d0) + 3.0d0*Sqrt(2.0d0)*u + Sqrt(6.0d0)*v + &
                   4.0d0*Sqrt(3.0d0)*w))/16.0d0
               EdgeBasis(8,3) = -((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v)*&
                   (-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v + 2.0d0*Sqrt(6.0d0)*w))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(8,1) = (-3.0d0*(-3.0d0*Sqrt(2.0d0) + 3.0d0*Sqrt(2.0d0)*u + &
                   Sqrt(6.0d0)*v + 4.0d0*Sqrt(3.0d0)*w))/16.0d0
               CurlBasis(8,2) = (9.0d0*(-Sqrt(6.0d0) + Sqrt(6.0d0)*u + Sqrt(2.0d0)*v + 4.0d0*w))/16.0d0
               CurlBasis(8,3) = 0.0d0
             ELSE
               k = 4
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(k,1) = (Sqrt(1.5d0)*w)/8.0d0
             EdgeBasis(k,2) = w/(8.0d0*Sqrt(2.0d0))
             EdgeBasis(k,3) = (Sqrt(6.0d0) - Sqrt(6.0d0)*u - Sqrt(2.0d0)*v)/16.0d0
             CurlBasis(k,1) = -1.0d0/(4.0d0*Sqrt(2.0d0))
             CurlBasis(k,2) = Sqrt(1.5d0)/4.0d0
             CurlBasis(k,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the fifth edge:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 9
               EdgeBasis(10,1) = (3.0d0*(Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - 4.0d0*w)*w)/16.0d0
               EdgeBasis(10,2) = (w*(-3.0d0*Sqrt(2.0d0) - 3.0d0*Sqrt(2.0d0)*u + &
                   Sqrt(6.0d0)*v + 4.0d0*Sqrt(3.0d0)*w))/16.0d0
               EdgeBasis(10,3) = ((Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v)*&
                   (-3.0d0 - 3.0d0*u + Sqrt(3.0d0)*v + 2.0d0*Sqrt(6.0d0)*w))/16.0d0
               CurlBasis(10,1) = (3.0d0*(3.0d0*Sqrt(2.0d0) + 3.0d0*Sqrt(2.0d0)*u - &
                   Sqrt(6.0d0)*v - 4.0d0*Sqrt(3.0d0)*w))/16.0d0
               CurlBasis(10,2) = (9.0d0*(Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - 4.0d0*w))/16.0d0
               CurlBasis(10,3) = 0.0d0
             ELSE
               k = 5
             END IF

             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(k,1) = -(Sqrt(1.5d0)*w)/8.0d0
             EdgeBasis(k,2) = w/(8.0d0*Sqrt(2.0d0))
             EdgeBasis(k,3) = (Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v)/16.0d0
             CurlBasis(k,1) = -1.0d0/(4.0d0*Sqrt(2.0d0))
             CurlBasis(k,2) = -Sqrt(1.5d0)/4.0d0
             CurlBasis(k,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the sixth edge:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 11
               EdgeBasis(12,1) = 0.0d0
               EdgeBasis(12,2) = (Sqrt(3.0d0)*(Sqrt(2.0d0)*v - 2.0d0*w)*w)/4.0d0
               EdgeBasis(12,3) = (Sqrt(1.5d0)*v*(-v + Sqrt(2.0d0)*w))/2.0d0
               CurlBasis(12,1) = (-3.0d0*(Sqrt(6.0d0)*v - 2.0d0*Sqrt(3.0d0)*w))/4.0d0
               CurlBasis(12,2) = 0.0d0
               CurlBasis(12,3) = 0.0d0
             ELSE
               k = 6
             END IF

             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(k,1) = 0.0d0
             EdgeBasis(k,2) = -w/(4.0d0*Sqrt(2.0d0))
             EdgeBasis(k,3) = v/(4.0d0*Sqrt(2.0d0))
             CurlBasis(k,1) = 1.0d0/(2.0d0*Sqrt(2.0d0))
             CurlBasis(k,2) = 0.0d0
             CurlBasis(k,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
             END IF

             ! -------------------------------------------------------------
             ! Finally scale the lowest-order basis functions so that 
             ! (b,t) = 1 when the integration is done over the element edge.
             ! -------------------------------------------------------------
             IF (SecondOrder) THEN
               DO k=1,6
                 EdgeBasis(2*(k-1)+1,:) = 2.0d0 * EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = 2.0d0 * CurlBasis(2*(k-1)+1,:)
               END DO
             ELSE
               DO k=1,6
                 EdgeBasis(k,:) = 2.0d0 * EdgeBasis(k,:)
                 CurlBasis(k,:) = 2.0d0 * CurlBasis(k,:)
               END DO
             END IF

             IF (SecondOrder) THEN
               !-------------------------------------------------
               ! Two basis functions defined on the face 213:
               !-------------------------------------------------
               TriangleFaceMap(:) = (/ 2,1,3 /)          
               Ind => Element % Nodeindexes

               DO j=1,3
                 FaceIndeces(j) = Ind(TriangleFaceMap(j))
               END DO
               IF (Parallel) THEN
                 DO j=1,3
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                 END DO
               END IF
               CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

               WorkBasis(1,1) = ((4.0d0*v - Sqrt(2.0d0)*w)*&
                   (-6.0d0 + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(48.0d0*Sqrt(3.0d0))
               WorkBasis(1,2) = -(u*(4.0d0*v - Sqrt(2.0d0)*w))/24.0d0
               WorkBasis(1,3) = (u*(-2.0d0*Sqrt(2.0d0)*v + w))/24.0d0
               WorkCurlBasis(1,1) = -u/(4.0d0*Sqrt(2.0d0))
               WorkCurlBasis(1,2) = (Sqrt(6.0d0) + 3.0d0*Sqrt(2.0d0)*v - 3.0d0*w)/24.0d0
               WorkCurlBasis(1,3) = (Sqrt(3.0d0) - 3.0d0*v)/6.0d0

               WorkBasis(2,1) = ((4.0d0*v - Sqrt(2.0d0)*w)*(-6.0d0 + 6.0d0*u + &
                   2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(96.0d0*Sqrt(3.0d0))
               WorkBasis(2,2) = -((4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u - 3.0d0*Sqrt(2.0d0)*w)*&
                   (-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/288.0d0
               WorkBasis(2,3) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*&
                   (-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(144.0d0*Sqrt(2.0d0))
               WorkCurlBasis(2,1) = -(-6.0d0 + 2.0d0*u + 2.0d0*Sqrt(3.0d0)*v + &
                   Sqrt(6.0d0)*w)/(16.0d0*Sqrt(2.0d0))
               WorkCurlBasis(2,2) = (2.0d0*Sqrt(3.0d0) - 6.0d0*Sqrt(3.0d0)*u + &
                   6.0d0*v - 3.0d0*Sqrt(2.0d0)*w)/(48.0d0*Sqrt(2.0d0))
               WorkCurlBasis(2,3) = (Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u - 3.0d0*v)/12.0d0

               WorkBasis(3,1) = -((4.0d0*v - Sqrt(2.0d0)*w)*(-6.0d0 - 6.0d0*u + &
                   2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(96.0d0*Sqrt(3.0d0))
               WorkBasis(3,2) = ((-4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u + 3.0d0*Sqrt(2.0d0)*w)* &
                   (-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/288.0d0
               WorkBasis(3,3) = -((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v)* &
                   (-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(144.0d0*Sqrt(2.0d0))
               WorkCurlBasis(3,1) = -(-6.0d0 - 2.0d0*u + 2.0d0*Sqrt(3.0d0)*v + &
                   Sqrt(6.0d0)*w)/(16.0d0*Sqrt(2.0d0))
               WorkCurlBasis(3,2) = (-2.0d0*Sqrt(3.0d0) - 6.0d0*Sqrt(3.0d0)*u - 6.0d0*v + &
                   3.0d0*Sqrt(2.0d0)*w)/(48.0d0*Sqrt(2.0d0))
               WorkCurlBasis(3,3) = (-Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u + 3.0d0*v)/12.0d0

               EdgeBasis(13,:) = D1 * WorkBasis(I1,:)
               CurlBasis(13,:) = D1 * WorkCurlBasis(I1,:)
               EdgeBasis(14,:) = D2 * WorkBasis(I2,:)
               CurlBasis(14,:) = D2 * WorkCurlBasis(I2,:)  

               !-------------------------------------------------
               ! Two basis functions defined on the face 124:
               !-------------------------------------------------
               TriangleFaceMap(:) = (/ 1,2,4 /)          
               Ind => Element % Nodeindexes

               DO j=1,3
                 FaceIndeces(j) = Ind(TriangleFaceMap(j))
               END DO
               IF (Parallel) THEN
                 DO j=1,3
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                 END DO
               END IF
               CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

               WorkBasis(1,1) = -(w*(-6.0d0 + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(8.0d0*Sqrt(6.0d0))
               WorkBasis(1,2) = (u*w)/(4.0d0*Sqrt(2.0d0))
               WorkBasis(1,3) = (u*w)/8.0d0
               WorkCurlBasis(1,1) = -u/(4.0d0*Sqrt(2.0d0))
               WorkCurlBasis(1,2) = (Sqrt(6.0d0) - Sqrt(2.0d0)*v - 3.0d0*w)/8.0d0
               WorkCurlBasis(1,3) = w/(2.0d0*Sqrt(2.0d0))

               WorkBasis(2,1) = -(w*(-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + &
                   Sqrt(6.0d0)*w))/(16.0d0*Sqrt(6.0d0))
               WorkBasis(2,2) = (w*(1.0d0 + u - v/Sqrt(3.0d0) - w/Sqrt(6.0d0)))/(8.0d0*Sqrt(2.0d0))
               WorkBasis(2,3) = ((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v)* &
                   (-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(48.0d0*Sqrt(2.0d0))
               WorkCurlBasis(2,1) = (-3.0d0*Sqrt(2.0d0) - Sqrt(2.0d0)*u + Sqrt(6.0d0)*v + Sqrt(3.0d0)*w)/16.0d0
               WorkCurlBasis(2,2) = (Sqrt(6.0d0) + 3.0d0*Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - 3.0d0*w)/16.0d0
               WorkCurlBasis(2,3) =  w/(4.0d0*Sqrt(2.0d0))

               WorkBasis(3,1) = (w*(-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(16.0d0*Sqrt(6.0d0))
               WorkBasis(3,2) = -(w*(-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(48.0d0*Sqrt(2.0d0))
               WorkBasis(3,3) = -((Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v)*&
                   (-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/96.0d0
               WorkCurlBasis(3,1) = (-3.0d0*Sqrt(2.0d0) + Sqrt(2.0d0)*u + Sqrt(6.0d0)*v + Sqrt(3.0d0)*w)/16.0d0
               WorkCurlBasis(3,2) = (-Sqrt(6.0d0) + 3.0d0*Sqrt(6.0d0)*u + Sqrt(2.0d0)*v + 3.0d0*w)/16.0d0
               WorkCurlBasis(3,3) = -w/(4.0d0*Sqrt(2.0d0))

               EdgeBasis(15,:) = D1 * WorkBasis(I1,:)
               CurlBasis(15,:) = D1 * WorkCurlBasis(I1,:)
               EdgeBasis(16,:) = D2 * WorkBasis(I2,:)
               CurlBasis(16,:) = D2 * WorkCurlBasis(I2,:)  

               !-------------------------------------------------
               ! Two basis functions defined on the face 234:
               !-------------------------------------------------
               TriangleFaceMap(:) = (/ 2,3,4 /)          
               Ind => Element % Nodeindexes

               DO j=1,3
                 FaceIndeces(j) = Ind(TriangleFaceMap(j))
               END DO
               IF (Parallel) THEN
                 DO j=1,3
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                 END DO
               END IF
               CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

               WorkBasis(1,1) = (w*(-2.0d0*Sqrt(2.0d0)*v + w))/16.0d0
               WorkBasis(1,2) = (w*(4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u - &
                   3.0d0*Sqrt(2.0d0)*w))/(16.0d0*Sqrt(6.0d0))
               WorkBasis(1,3) = -((1.0d0 + u - Sqrt(3.0d0)*v)*w)/16.0d0
               WorkCurlBasis(1,1) = (-2.0d0*Sqrt(2.0d0) - 2.0d0*Sqrt(2.0d0)*u + 3.0d0*Sqrt(3.0d0)*w)/16.0d0
               WorkCurlBasis(1,2) = (-2.0d0*Sqrt(2.0d0)*v + 3.0d0*w)/16.0d0
               WorkCurlBasis(1,3) = w/(2.0d0*Sqrt(2.0d0))

               WorkBasis(2,1) = (w*(-2.0d0*Sqrt(2.0d0)*v + w))/16.0d0
               WorkBasis(2,2) = -(w*(-4.0d0*v + Sqrt(2.0d0)*w))/(16.0d0*Sqrt(6.0d0))
               WorkBasis(2,3) = -((Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v)*&
                   (-4.0d0*v + Sqrt(2.0d0)*w))/(32.0d0*Sqrt(3.0d0))
               WorkCurlBasis(2,1) = (2.0d0*Sqrt(2.0d0) + 2.0d0*Sqrt(2.0d0)*u - &
                   2.0d0*Sqrt(6.0d0)*v + Sqrt(3.0d0)*w)/16.0d0
               WorkCurlBasis(2,2) = (-4.0d0*Sqrt(2.0d0)*v + 3.0d0*w)/16.0d0
               WorkCurlBasis(2,3) = w/(4.0d0*Sqrt(2.0d0))

               WorkBasis(3,1) = 0.0d0
               WorkBasis(3,2) = (w*(-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(24.0d0*Sqrt(2.0d0))
               WorkBasis(3,3) = -(v*(-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(24.0d0*Sqrt(2.0d0))
               WorkCurlBasis(3,1) = (2.0d0*Sqrt(2.0d0) + 2.0d0*Sqrt(2.0d0)*u - Sqrt(6.0d0)*v - Sqrt(3.0d0)*w)/8.0d0
               WorkCurlBasis(3,2) = -v/(4.0d0*Sqrt(2.0d0))
               WorkCurlBasis(3,3) = -w/(4.0d0*Sqrt(2.0d0))

               EdgeBasis(17,:) = D1 * WorkBasis(I1,:)
               CurlBasis(17,:) = D1 * WorkCurlBasis(I1,:)
               EdgeBasis(18,:) = D2 * WorkBasis(I2,:)
               CurlBasis(18,:) = D2 * WorkCurlBasis(I2,:)  

               !-------------------------------------------------
               ! Two basis functions defined on the face 314:
               !-------------------------------------------------
               TriangleFaceMap(:) = (/ 3,1,4 /)          
               Ind => Element % Nodeindexes

               DO j=1,3
                 FaceIndeces(j) = Ind(TriangleFaceMap(j))
               END DO
               IF (Parallel) THEN
                 DO j=1,3
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                 END DO
               END IF
               CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

               WorkBasis(1,1) = (w*(-2.0d0*Sqrt(2.0d0)*v + w))/16.0d0
               WorkBasis(1,2) = (w*(-4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u + &
                   3.0d0*Sqrt(2.0d0)*w))/(16.0d0*Sqrt(6.0d0))
               WorkBasis(1,3) = -((-1.0d0 + u + Sqrt(3.0d0)*v)*w)/16.0d0
               WorkCurlBasis(1,1) = (2.0d0*Sqrt(2.0d0) - 2.0d0*Sqrt(2.0d0)*u - 3.0d0*Sqrt(3.0d0)*w)/16.0d0
               WorkCurlBasis(1,2) = (-2.0d0*Sqrt(2.0d0)*v + 3.0d0*w)/16.0d0
               WorkCurlBasis(1,3) = w/(2.0d0*Sqrt(2.0d0))

               WorkBasis(2,1) = 0.0d0
               WorkBasis(2,2) = (w*(-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(24.0d0*Sqrt(2.0d0))
               WorkBasis(2,3) = -(v*(-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(24.0d0*Sqrt(2.0d0))
               WorkCurlBasis(2,1) = (2.0d0*Sqrt(2.0d0) - 2.0d0*Sqrt(2.0d0)*u - Sqrt(6.0d0)*v - Sqrt(3.0d0)*w)/8.0d0
               WorkCurlBasis(2,2) = v/(4.0d0*Sqrt(2.0d0))
               WorkCurlBasis(2,3) =  w/(4.0d0*Sqrt(2.0d0))

               WorkBasis(3,1) = ((2.0d0*Sqrt(2.0d0)*v - w)*w)/16.0d0
               WorkBasis(3,2) = -(w*(-4.0d0*v + Sqrt(2.0d0)*w))/(16.0d0*Sqrt(6.0d0))
               WorkBasis(3,3) = ((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v)*&
                   (-4.0d0*v + Sqrt(2.0d0)*w))/(16.0d0*Sqrt(6.0d0))
               WorkCurlBasis(3,1) = (2.0d0*Sqrt(2.0d0) - 2.0d0*Sqrt(2.0d0)*u - &
                   2.0d0*Sqrt(6.0d0)*v + Sqrt(3.0d0)*w)/16.0d0
               WorkCurlBasis(3,2) = (4.0d0*Sqrt(2.0d0)*v - 3.0d0*w)/16.0d0
               WorkCurlBasis(3,3) =  -w/(4.0d0*Sqrt(2.0d0))

               EdgeBasis(19,:) = D1 * WorkBasis(I1,:)
               CurlBasis(19,:) = D1 * WorkCurlBasis(I1,:)
               EdgeBasis(20,:) = D2 * WorkBasis(I2,:)
               CurlBasis(20,:) = D2 * WorkCurlBasis(I2,:)                  
             END IF
           END IF
           
         CASE(6)
           !--------------------------------------------------------------
           ! This branch is for handling pyramidic elements
           !--------------------------------------------------------------         
           EdgeMap => LGetEdgeMap(6)
           Ind => Element % Nodeindexes

           IF (SecondOrder) THEN
             EdgeSign = 1.0d0

             LBasis(1) = 0.1D1 / 0.4D1 - u / 0.4D1 - v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1 + &
                 u * v / ( (0.1D1 - w * sqrt(0.2D1) / 0.2D1) * 0.4D1 )
             LBasis(2) = 0.1D1 / 0.4D1 + u / 0.4D1 - v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1 - &
                 u * v / ( (0.1D1 - w * sqrt(0.2D1) / 0.2D1) * 0.4D1 )
             LBasis(3) = 0.1D1 / 0.4D1 + u / 0.4D1 + v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1 + &
                 u * v / ( (0.1D1 - w * sqrt(0.2D1) / 0.2D1) * 0.4D1 )
             LBasis(4) = 0.1D1 / 0.4D1 - u / 0.4D1 + v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1 - &
                 u * v / ( (0.1D1 - w * sqrt(0.2D1) / 0.2D1) * 0.4D1 )
             LBasis(5) = w * sqrt(0.2D1) / 0.2D1

             Beta(1) = 0.1D1 / 0.2D1 - u / 0.2D1 - w * sqrt(0.2D1) / 0.4D1
             Beta(2) = 0.1D1 / 0.2D1 - v / 0.2D1 - w * sqrt(0.2D1) / 0.4D1
             Beta(3) = 0.1D1 / 0.2D1 + u / 0.2D1 - w * sqrt(0.2D1) / 0.4D1
             Beta(4) = 0.1D1 / 0.2D1 + v / 0.2D1 - w * sqrt(0.2D1) / 0.4D1

             ! Edge 12:
             !--------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(1,1) = 0.1D1 / 0.4D1 - v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1
             EdgeBasis(1,2) = 0.0d0
             EdgeBasis(1,3) = sqrt(0.2D1) * u * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 ((w * sqrt(0.2D1) - 0.2D1) * 0.8D1)
             CurlBasis(1,1) = sqrt(0.2D1) * u / ((w * sqrt(0.2D1) - 0.2D1) * 0.4D1)
             CurlBasis(1,2) = -sqrt(0.2D1) / 0.8D1 - sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             CurlBasis(1,3) = 0.1D1 / 0.4D1
             IF (nj<ni) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,:) = -CurlBasis(1,:)
               EdgeSign(1) = -1.0d0
             END IF

             EdgeBasis(2,1:3) = 3.0d0 * u * EdgeBasis(1,1:3)
             CurlBasis(2,1) = 0.3D1 / 0.4D1 * u ** 2 * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(2,2) = -0.3D1 / 0.8D1 * u * sqrt(0.2D1) * (0.3D1 * w * sqrt(0.2D1) + &
                 4.0D0 * v - 0.6D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(2,3) = 0.3D1 / 0.4D1 * u

             ! Edge 23:
             !--------------------------------------------------------------
             k = 3 ! k=2 for first-order
             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)

             EdgeBasis(k,1) = 0.0d0
             EdgeBasis(k,2) = 0.1D1 / 0.4D1 + u / 0.4D1 - w * sqrt(0.2D1) / 0.8D1
             EdgeBasis(k,3) = sqrt(0.2D1) * v * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             CurlBasis(k,1) = sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 ) + sqrt(0.2D1) /  0.8D1
             CurlBasis(k,2) = sqrt(0.2D1) * v / ( (w * sqrt(0.2D1) - 0.2D1) * 0.4D1 )
             CurlBasis(k,3) = 0.1D1 / 0.4D1
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * v * EdgeBasis(k,1:3)
             CurlBasis(k+1,1) = 0.3D1 / 0.8D1 * v * sqrt(0.2D1) * (0.3D1 * w * sqrt(0.2D1) - & 
                 4.0D0 * u - 0.6D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,2) = 0.3D1 / 0.4D1 * v ** 2 * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,3) = 0.3D1 / 0.4D1 * v

             ! Edge 43:
             !--------------------------------------------------------------
             k = 5 ! k=3 for first-order
             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)

             EdgeBasis(k,1) = 0.1D1 / 0.4D1 + v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1
             EdgeBasis(k,2) = 0.0d0
             EdgeBasis(k,3) = sqrt(0.2D1) * u * (w * sqrt(0.2D1) - 2.0D0 * v - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )

             CurlBasis(k,1) = -sqrt(0.2D1) * u / ( (w * sqrt(0.2D1) - 0.2D1) * 0.4D1 )
             CurlBasis(k,2) = -sqrt(0.2D1) / 0.8D1 - sqrt(0.2D1) * (w * sqrt(0.2D1) - &
                 2.0D0 * v - 0.2D1) / ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             CurlBasis(k,3) = -0.1D1 / 0.4D1
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * u * EdgeBasis(k,1:3)
             CurlBasis(k+1,1) = -0.3D1 / 0.4D1 * u ** 2 * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,2) = -0.3D1 / 0.8D1 * u * sqrt(0.2D1) * (0.3D1 * w * sqrt(0.2D1) - &
                 4.0D0 * v - 0.6D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,3) = -0.3D1 / 0.4D1 * u


             ! Edge 14:
             !--------------------------------------------------------------
             k = 7 ! k=4 for first-order
             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)

             EdgeBasis(k,1) = 0.0d0 
             EdgeBasis(k,2) = 0.1D1 / 0.4D1 - u / 0.4D1 - w * sqrt(0.2D1) / 0.8D1
             EdgeBasis(k,3) = sqrt(0.2D1) * v * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / & 
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )

             CurlBasis(k,1) = sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / ( (w * &
                 sqrt(0.2D1) - 0.2D1) * 0.8D1 ) + sqrt(0.2D1) / 0.8D1
             CurlBasis(k,2) = -sqrt(0.2D1) * v / ( (w * sqrt(0.2D1) - 0.2D1) * 0.4D1 )
             CurlBasis(k,3) = -0.1D1 / 0.4D1
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * v * EdgeBasis(k,1:3)
             CurlBasis(k+1,1) = 0.3D1 / 0.8D1 * v * sqrt(0.2D1) * (0.3D1 * w * sqrt(0.2D1) + &
                 4.0D0 * u - 0.6D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,2) = -0.3D1 / 0.4D1 * v ** 2 * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,3) = -0.3D1 / 0.4D1 * v


             ! Edge 15:
             !--------------------------------------------------------------
             k = 9 ! k=5 for first-order             
             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)

             EdgeBasis(k,1) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             EdgeBasis(k,2) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             EdgeBasis(k,3) = -sqrt(0.2D1)/ 0.4D1 * (0.2D1 * sqrt(0.2D1) * u * v * w - &
                 0.2D1 * sqrt(0.2D1) * u * w - &
                 0.2D1 * sqrt(0.2D1) * v * w + u * w ** 2 + v * w ** 2 + 0.2D1 * w * sqrt(0.2D1) - &
                 0.2D1 * u * v - w ** 2 + 0.2D1 * u + 0.2D1 * v - 0.2D1) / (w * sqrt(0.2D1) - 0.2D1) ** 2 

             CurlBasis(k,1) = (-sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - 0.2D1 * &
                 u * w - 0.2D1 * sqrt(0.2D1) + 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,2) = -(-sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - 0.2D1 * &
                 v * w - 0.2D1 * sqrt(0.2D1) + 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,3) = 0.0d0 
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * EdgeSign(k) * EdgeBasis(k,1:3) * ( LBasis(5)-LBasis(1)+LBasis(3) )

             CurlBasis(k+1,1) = 0.3D1 / 0.8D1 * (-0.9D1 * sqrt(0.2D1) * u * w ** 2 - &
                 0.3D1 * sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 + &
                 0.6D1 * u * v * sqrt(0.2D1) + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w - &
                 0.8D1 * u * v * w - 0.6D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) - 0.6D1 * v * sqrt(0.2D1) + &
                 0.24D2 * u * w + 0.12D2 * v * w + 0.2D1 * sqrt(0.2D1) - 0.16D2 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1)**2
             CurlBasis(k+1,2) = -0.3D1 / 0.8D1 * (-0.3D1 * sqrt(0.2D1) * u * w ** 2 - &
                 0.9D1 * sqrt(0.2D1) * v * w ** 2 + 0.6D1 * u * v * sqrt(0.2D1) + &
                 0.4D1 * sqrt(0.2D1) * v ** 2 + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u* v * w - &
                 0.4D1 * v ** 2 * w - 0.6D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) - 0.6D1 * v * sqrt(0.2D1) + &
                 0.12D2 * u * w + 0.24D2 * v * w + 0.2D1 * sqrt(0.2D1) - 0.16D2 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1)**2
             CurlBasis(k+1,3) = 0.3D1 / 0.8D1 * w * sqrt(0.2D1) * (u - v) / (w * sqrt(0.2D1) - 0.2D1)


             ! Edge 25:
             !--------------------------------------------------------------
             k = 11 ! k=6 for first-order  
             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)

             EdgeBasis(k,1) = -w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             EdgeBasis(k,2) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             EdgeBasis(k,3) = sqrt(0.2D1)/ 0.4D1 * (0.2D1 * sqrt(0.2D1) * u * v * w - 0.2D1 * &
                 sqrt(0.2D1) * u * w + 0.2D1 * sqrt(0.2D1) * v * w + u * w ** 2 - v * w ** 2 - &
                 0.2D1 * w * sqrt(0.2D1) - 0.2D1 * u * v + w ** 2 + 0.2D1 * u - 0.2D1 * v + 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2 
             CurlBasis(k,1) = -(sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - 0.2D1 * u * w + &
                 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,2) = (-sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - 0.2D1 * & 
                 v * w - 0.2D1 * sqrt(0.2D1) + 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,3) = 0.0d0 
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * EdgeSign(k) * EdgeBasis(k,1:3) * ( LBasis(5)-LBasis(2)+LBasis(4) )

             CurlBasis(k+1,1) = 0.3D1 / 0.8D1 * (0.9D1 * sqrt(0.2D1) * u * w ** 2 - &
                 0.3D1 * sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.6D1 * u * v * sqrt(0.2D1) + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u** 2 * w + &
                 0.8D1 * u * v * w - 0.6D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) - &
                 0.6D1 * v * sqrt(0.2D1) - 0.24D2 * u * w + 0.12D2 * v * w + 0.2D1 * sqrt(0.2D1) - &
                 0.16D2 * w) / (w * sqrt(0.2D1) - 0.2D1)**2
             CurlBasis(k+1,2) = -0.3D1 / 0.8D1 * (-0.3D1 * sqrt(0.2D1) * u * w ** 2 + &
                 0.9D1 * sqrt(0.2D1) * v * w ** 2 + 0.6D1 * u * v * sqrt(0.2D1) - &
                 0.4D1 * sqrt(0.2D1) * v ** 2 - 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w + &
                 0.4D1 * v ** 2 * w + 0.6D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) + &
                 0.6D1 * v * sqrt(0.2D1) + 0.12D2 * u * w - 0.24D2 * v * w - 0.2D1 * sqrt(0.2D1) + &
                 0.16D2 * w) / (w * sqrt(0.2D1) - 0.2D1)** 2
             CurlBasis(k+1,3) = 0.3D1 / 0.8D1 * w * sqrt(0.2D1) * (u + v) / (w * sqrt(0.2D1) - 0.2D1)


             ! Edge 35:
             !--------------------------------------------------------------
             k = 13 ! k=7 for first-order  
             i = EdgeMap(7,1)
             j = EdgeMap(7,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)

             EdgeBasis(k,1) = -w * sqrt(0.2D1)/ 0.8D1 * (w * sqrt(0.2D1) - 2.0D0 * v - 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1) 
             EdgeBasis(k,2) = -w * sqrt(0.2D1) / 0.8D1 * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / & 
                 (w * sqrt(0.2D1) - 0.2D1)
             EdgeBasis(k,3) = -sqrt(0.2D1)/ 0.4D1 * (0.2D1 * sqrt(0.2D1) * u * v * w + 0.2D1 * &
                 sqrt(0.2D1) * u * w + 0.2D1 * sqrt(0.2D1) * v * w - u * w ** 2 - v * w ** 2 + &
                 0.2D1 * w * sqrt(0.2D1) - 0.2D1 * u * v - w ** 2 - 0.2D1 * u - 0.2D1 * v - 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2 
             CurlBasis(k,1) = (sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - 0.2D1 * u * w + &
                 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,2) = -(sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - 0.2D1 * &
                 v * w + 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,3) = 0.0d0 
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * EdgeSign(k) * EdgeBasis(k,1:3) * ( LBasis(5)-LBasis(3)+LBasis(1) )

             CurlBasis(k+1,1) = -0.3D1 / 0.8D1 * (0.9D1 * sqrt(0.2D1) * u * w ** 2 + &
                 0.3D1 * sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 + &
                 0.6D1 * u * v * sqrt(0.2D1) + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w - &
                 0.8D1 * u * v * w - 0.6D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) + &
                 0.6D1 * v * sqrt(0.2D1) - 0.24D2 * u * w - 0.12D2 * v * w + 0.2D1 * sqrt(0.2D1) - &
                 0.16D2 * w) / (w * sqrt(0.2D1) - 0.2D1)**2
             CurlBasis(k+1,2) = 0.3D1 / 0.8D1 * (0.3D1 * sqrt(0.2D1) * u * w ** 2 + &
                 0.9D1 * sqrt(0.2D1) * v * w ** 2 + 0.6D1 * u * v * sqrt(0.2D1) + &
                 0.4D1 * sqrt(0.2D1) * v ** 2 + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u *v * w - &
                 0.4D1 * v ** 2 * w - 0.6D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) + 0.6D1 * v * sqrt(0.2D1) - &
                 0.12D2 * u * w - 0.24D2 * v * w + 0.2D1 * sqrt(0.2D1) - 0.16D2 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2
             CurlBasis(k+1,3) = -0.3D1 / 0.8D1 * w * sqrt(0.2D1) * (u - v) / (w * sqrt(0.2D1) - 0.2D1)


             ! Edge 45:
             !--------------------------------------------------------------
             k = 15 ! k=8 for first-order  
             i = EdgeMap(8,1)
             j = EdgeMap(8,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)

             EdgeBasis(k,1) = w * sqrt(0.2D1) / 0.8D1 * (w * sqrt(0.2D1) - 2.0D0 * v - 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1) 
             EdgeBasis(k,2) = -w * sqrt(0.2D1) / 0.8D1 * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1)
             EdgeBasis(k,3) = sqrt(0.2D1) / 0.4D1 * (0.2D1 * sqrt(0.2D1) * u * v * w + &
                 0.2D1 * sqrt(0.2D1) * u * w - 0.2D1 * sqrt(0.2D1) * v * w - u * w ** 2 + v * w ** 2 - &
                 0.2D1 * w * sqrt(0.2D1) - 0.2D1 * u * v + w ** 2 - 0.2D1 * u + 0.2D1 * v + 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2 
             CurlBasis(k,1) = -(-sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - 0.2D1 * u * w - &
                 0.2D1 * sqrt(0.2D1) + 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1)** 2 * 0.2D1 )
             CurlBasis(k,2) = (sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - 0.2D1 * v * w + &
                 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1)** 2 * 0.2D1 )
             CurlBasis(k,3) = 0.0d0 
             IF (nj<ni) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * EdgeSign(k) * EdgeBasis(k,1:3) * ( LBasis(5)-LBasis(4)+LBasis(2) )

             CurlBasis(k+1,1) = -0.3D1 / 0.8D1 * (-0.9D1 * sqrt(0.2D1) * u * w ** 2 + &
                 0.3D1 * sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.6D1 * u * v * sqrt(0.2D1) + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u** 2 * w + &
                 0.8D1 * u * v * w - 0.6D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) + &
                 0.6D1 * v * sqrt(0.2D1) + 0.24D2 * u * w - 0.12D2 * v * w + 0.2D1 * sqrt(0.2D1) - &
                 0.16D2 * w) / (w * sqrt(0.2D1) - 0.2D1) ** 2
             CurlBasis(k+1,2) = 0.3D1 / 0.8D1 * (0.3D1 * sqrt(0.2D1) * u * w ** 2 - &
                 0.9D1 * sqrt(0.2D1) * v * w ** 2 + 0.6D1 * u * v * sqrt(0.2D1) - &
                 0.4D1 * sqrt(0.2D1) * v ** 2 - 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u *v * w + &
                 0.4D1 * v ** 2 * w + 0.6D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) - &
                 0.6D1 * v * sqrt(0.2D1) - 0.12D2 * u * w + 0.24D2 * v * w - 0.2D1 * sqrt(0.2D1) + &
                 0.16D2 * w) / (w * sqrt(0.2D1) - 0.2D1)**2
             CurlBasis(k+1,3) = -0.3D1 / 0.8D1 * w * sqrt(0.2D1) * (u + v) / (w * sqrt(0.2D1) - 0.2D1)


             ! Square face:
             ! ------------------------------------------------------------------
             SquareFaceMap(:) = (/ 1,2,3,4 /)

             WorkBasis(1,1:3) = 2.0d0 * ( EdgeSign(1) * EdgeBasis(1,1:3) * Beta(4) + &
                 EdgeSign(5) * EdgeBasis(5,1:3) * Beta(2) ) / (1.0d0 - LBasis(5))
             WorkCurlBasis(1,1) = -0.2D1 * u * v * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(1,2) = -(sqrt(0.2D1) * w ** 2 + 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / & 
                 (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(1,3) = -0.2D1 * v / (w * sqrt(0.2D1) - 0.2D1)

             WorkBasis(2,1:3) = 3.0d0 * WorkBasis(1,1:3) * u
             WorkCurlBasis(2,1) = -0.6D1 * u ** 2 * sqrt(0.2D1) * v / (w * sqrt(0.2D1) - 0.2D1)** 2
             WorkCurlBasis(2,2) = 0.3D1 / 0.2D1 * u * (0.2D1 * sqrt(0.2D1) * v ** 2 - &
                 0.3D1 * sqrt(0.2D1) * w ** 2 - 0.6D1 * sqrt(0.2D1) + 0.12D2 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(2,3) = -0.6D1 * u * v / (w * sqrt(0.2D1) - 0.2D1)

             WorkBasis(3,1:3) = 2.0d0 * ( EdgeSign(3) * EdgeBasis(3,1:3) * Beta(1) + &
                 EdgeSign(7) * EdgeBasis(7,1:3) * Beta(3) ) / (1.0d0 - LBasis(5))
             WorkCurlBasis(3,1) = (sqrt(0.2D1) * w ** 2 + 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(3,2) = 0.2D1 * u * v * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(3,3) = 0.2D1 * u / (w * sqrt(0.2D1) - 0.2D1)

             WorkBasis(4,1:3) = 3.0d0 * WorkBasis(3,1:3) * v
             WorkCurlBasis(4,1) = -0.3D1 / 0.2D1 * v * (0.2D1 * sqrt(0.2D1) * u ** 2 - &
                 0.3D1 * sqrt(0.2D1) * w ** 2 - 0.6D1 * sqrt(0.2D1) + 0.12D2 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(4,2) = 0.6D1 * sqrt(0.2D1) * v ** 2 * u / (w * sqrt(0.2D1) - 0.2D1)**2
             WorkCurlBasis(4,3) = 0.6D1 * u * v / (w * sqrt(0.2D1) - 0.2D1)

             ! -------------------------------------------------------------------
             ! Finally apply an order change and sign reversions if needed. 
             ! -------------------------------------------------------------------
             DO j=1,4
               FaceIndeces(j) = Ind(SquareFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(17,:) = D1 * WorkBasis(2*(I1-1)+1,:)
             CurlBasis(17,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
             EdgeBasis(18,:) = WorkBasis(2*(I1-1)+2,:)
             CurlBasis(18,:) = WorkCurlBasis(2*(I1-1)+2,:)
             EdgeBasis(19,:) = D2 * WorkBasis(2*(I2-1)+1,:)
             CurlBasis(19,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
             EdgeBasis(20,:) = WorkBasis(2*(I2-1)+2,:)
             CurlBasis(20,:) = WorkCurlBasis(2*(I2-1)+2,:) 

             
             !-------------------------------------------------
             ! Two basis functions defined on the face 125:
             !-------------------------------------------------
             TriangleFaceMap(:) = (/ 1,2,5 /)          
 
             DO j=1,3
               FaceIndeces(j) = Ind(TriangleFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,3
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             WorkBasis(1,1:3) = LBasis(5) * EdgeSign(1) * EdgeBasis(1,1:3)
             WorkCurlBasis(1,1) = w * u / (w * sqrt(0.2D1) - 0.2D1) / 0.4D1
             WorkCurlBasis(1,2) = (-0.3D1 * sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - & 
                 0.4D1 * v * w - 0.2D1 * sqrt(0.2D1) + 0.8D1 * w) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             WorkCurlBasis(1,3) = w * sqrt(0.2D1) / 0.8D1

             WorkBasis(2,1:3) = Beta(3) * EdgeSign(9) * EdgeBasis(9,1:3)
             WorkCurlBasis(2,1) = (sqrt(0.2D1) * u * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w + 0.3D1 * w ** 3 - &
                 0.2D1 * u * w - 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(2,2) = -(-0.3D1 * sqrt(0.2D1) * u * w ** 2 + 0.2D1 * sqrt(0.2D1) * &
                 v * w ** 2 + 0.6D1 * u * v * sqrt(0.2D1) - 0.7D1 * sqrt(0.2D1) * w ** 2 - &
                 0.8D1 * u * v * w + 0.3D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) + 0.2D1 * v * sqrt(0.2D1) + &
                 0.12D2 * u * w - 0.6D1 * v * w - 0.2D1 * sqrt(0.2D1) + 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1)**2 )
             WorkCurlBasis(2,3) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.16D2 )

             WorkBasis(3,1:3) = Beta(1) * EdgeSign(11) * EdgeBasis(11,1:3)
             WorkCurlBasis(3,1) = (-sqrt(0.2D1) * u * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w + 0.3D1 * w ** 3 + &
                 0.2D1 * u * w - 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1)** 2 ) 
             WorkCurlBasis(3,2) = -(-0.3D1 * sqrt(0.2D1) * u * w ** 2 - 0.2D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * u * v * sqrt(0.2D1) + 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w - &
                 0.3D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) - 0.2D1 * v * sqrt(0.2D1) + 0.12D2 * u * w + &
                 0.6D1 * v * w + 0.2D1 * sqrt(0.2D1) - 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1)**2 ) 
             WorkCurlBasis(3,3) = -w * sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             EdgeBasis(21,:) = D1 * WorkBasis(I1,:)
             CurlBasis(21,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(22,:) = D2 * WorkBasis(I2,:)
             CurlBasis(22,:) = D2 * WorkCurlBasis(I2,:)              

             !-------------------------------------------------
             ! Two basis functions defined on the face 235:
             !-------------------------------------------------
             TriangleFaceMap(:) = (/ 2,3,5 /)          
 
             DO j=1,3
               FaceIndeces(j) = Ind(TriangleFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,3
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             WorkBasis(1,1:3) = LBasis(5) * EdgeSign(3) * EdgeBasis(3,1:3)
             WorkCurlBasis(1,1) = (0.3D1 * sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - 0.4D1 * u * w + &
                 0.2D1 * sqrt(0.2D1) - 0.8D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             WorkCurlBasis(1,2) = w * v / (w * sqrt(0.2D1) - 0.2D1) / 0.4D1
             WorkCurlBasis(1,3) = w * sqrt(0.2D1) / 0.8D1

             WorkBasis(2,1:3) = Beta(4) * EdgeSign(11) * EdgeBasis(11,1:3)
             WorkCurlBasis(2,1) = -(0.2D1 * sqrt(0.2D1) * u * w ** 2 + 0.3D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v + 0.7D1 * sqrt(0.2D1) * w** 2 - 0.8D1 * u * v * w - &
                 0.3D1 * w ** 3 + 0.2D1 * u * sqrt(0.2D1) + 0.6D1 * v * sqrt(0.2D1) - 0.6D1 * u * w - &
                 0.12D2 * w * v + 0.2D1 * sqrt(0.2D1) - 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2) 
             WorkCurlBasis(2,2) = (sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * v ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * v ** 2 * w + 0.3D1 * w ** 3 - 0.2D1 * w * v - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 )
             WorkCurlBasis(2,3) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             WorkBasis(3,1:3) = Beta(2) * EdgeSign(13) * EdgeBasis(13,1:3)
             WorkCurlBasis(3,1) = -(-0.2D1 * sqrt(0.2D1) * u * w ** 2 + 0.3D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v - 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w + &
                 0.3D1 * w ** 3 - 0.2D1 * u * sqrt(0.2D1) + 0.6D1 * v * sqrt(0.2D1) + 0.6D1 * u * w - &
                 0.12D2 * w * v - 0.2D1 * sqrt(0.2D1) + 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,2) = (-sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * v ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * v ** 2 * w + 0.3D1 * w ** 3 + 0.2D1 * w * v - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,3) = -w * sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * v - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.16D2 )

             EdgeBasis(23,:) = D1 * WorkBasis(I1,:)
             CurlBasis(23,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(24,:) = D2 * WorkBasis(I2,:)
             CurlBasis(24,:) = D2 * WorkCurlBasis(I2,:)              

             !-------------------------------------------------
             ! Two basis functions defined on the face 345:
             !-------------------------------------------------
             TriangleFaceMap(:) = (/ 3,4,5 /)          
 
             DO j=1,3
               FaceIndeces(j) = Ind(TriangleFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,3
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             WorkBasis(1,1:3) = -LBasis(5) * EdgeSign(5) * EdgeBasis(5,1:3)
             WorkCurlBasis(1,1) = w * u / (w * sqrt(0.2D1) - 0.2D1) / 0.4D1
             WorkCurlBasis(1,2) = (0.3D1 * sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - 0.4D1 * w * v + &
                 0.2D1 * sqrt(0.2D1) - 0.8D1 * w) / (0.8D1 * (w * sqrt(0.2D1)- 0.2D1) )
             WorkCurlBasis(1,3) = w * sqrt(0.2D1) / 0.8D1

             WorkBasis(2,1:3) = Beta(1) * EdgeSign(13) * EdgeBasis(13,1:3)
             WorkCurlBasis(2,1) = -(-sqrt(0.2D1) * u * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w + 0.3D1 * w ** 3 + 0.2D1 * u * w - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(2,2) = (0.3D1 * sqrt(0.2D1) * u * w ** 2 - 0.2D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v - 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w + &
                 0.3D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) - 0.2D1 * v * sqrt(0.2D1) - 0.12D2 * u * w + &
                 0.6D1 * w * v - 0.2D1 * sqrt(0.2D1) + 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(2,3) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             WorkBasis(3,1:3) = Beta(3) * EdgeSign(15) * EdgeBasis(15,1:3)
             WorkCurlBasis(3,1) = -(sqrt(0.2D1) * u * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w + 0.3D1 * w ** 3 - 0.2D1 * u * w - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,2) = (0.3D1 * sqrt(0.2D1) * u * w ** 2 + 0.2D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v + 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w - &
                 0.3D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) + 0.2D1 * v * sqrt(0.2D1) - 0.12D2 * u * w - &
                 0.6D1 * w * v + 0.2D1 * sqrt(0.2D1) - 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,3) = -w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             EdgeBasis(25,:) = D1 * WorkBasis(I1,:)
             CurlBasis(25,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(26,:) = D2 * WorkBasis(I2,:)
             CurlBasis(26,:) = D2 * WorkCurlBasis(I2,:)              

             !-------------------------------------------------
             ! Two basis functions defined on the face 415:
             !-------------------------------------------------
             TriangleFaceMap(:) = (/ 4,1,5 /)          
 
             DO j=1,3
               FaceIndeces(j) = Ind(TriangleFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,3
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             WorkBasis(1,1:3) = -LBasis(5) * EdgeSign(7) * EdgeBasis(7,1:3)
             WorkCurlBasis(1,1) = (-0.3D1 * sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - &
                 0.4D1 * u * w - 0.2D1 * sqrt(0.2D1) + 0.8D1 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) )
             WorkCurlBasis(1,2) = w * v / (w * sqrt(0.2D1) - 0.2D1) / 0.4D1
             WorkCurlBasis(1,3) = w * sqrt(0.2D1) / 0.8D1

             WorkBasis(2,1:3) = Beta(2) * EdgeSign(15) * EdgeBasis(15,1:3)
             WorkCurlBasis(2,1) = (-0.2D1 * sqrt(0.2D1) * u * w ** 2 - 0.3D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v + 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w - &
                 0.3D1 * w ** 3 - 0.2D1 * u * sqrt(0.2D1) - 0.6D1 * v * sqrt(0.2D1) + 0.6D1 * u * w + &
                 0.12D2 * w * v + 0.2D1 * sqrt(0.2D1) - 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 )
             WorkCurlBasis(2,2) = -(-sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * v ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * v ** 2 * w + 0.3D1 * w ** 3 + 0.2D1 * w * v - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(2,3) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * v - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             WorkBasis(3,1:3) = Beta(4) * EdgeSign(9) * EdgeBasis(9,1:3)
             WorkCurlBasis(3,1) = (0.2D1 * sqrt(0.2D1) * u * w ** 2 - 0.3D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v - 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w + &
                 0.3D1 * w ** 3 + 0.2D1 * u * sqrt(0.2D1) - 0.6D1 * v * sqrt(0.2D1) - 0.6D1 * u * w + &
                 0.12D2 * w * v - 0.2D1 * sqrt(0.2D1) + 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,2) = -(sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * v ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * v ** 2 * w + 0.3D1 * w ** 3 - 0.2D1 * w * v - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,3) = -w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             EdgeBasis(27,:) = D1 * WorkBasis(I1,:)
             CurlBasis(27,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(28,:) = D2 * WorkBasis(I2,:)
             CurlBasis(28,:) = D2 * WorkCurlBasis(I2,:)              


             ! Finally three interior basis functions:
             ! -----------------------------------------------------------------------------------
             EdgeBasis(29,1:3) = LBasis(5) * Beta(4) * EdgeSign(1) * EdgeBasis(1,1:3)
             CurlBasis(29,1) = u * v * w / (0.4D1 * (w * sqrt(0.2D1) - 0.2D1) ) 
             CurlBasis(29,2) = (0.2D1 * sqrt(0.2D1) * v ** 2 - 0.9D1 * sqrt(0.2D1) * w ** 2 - &
                 0.4D1 * v ** 2 * w + 0.4D1 * w ** 3 - 0.2D1 * sqrt(0.2D1) + 0.12D2 * w) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 
             CurlBasis(29,3) = sqrt(0.2D1) * v * w / 0.8D1

             EdgeBasis(30,1:3) = LBasis(5) * Beta(3) * EdgeSign(7) * EdgeBasis(7,1:3)
             CurlBasis(30,1) = -(0.2D1 * sqrt(0.2D1) * u ** 2 - 0.9D1 * sqrt(0.2D1) * w **2 - &
                 0.4D1 * u ** 2 * w + 0.4D1 * w ** 3 - 0.2D1 * sqrt(0.2D1) + 0.12D2 * w) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 
             CurlBasis(30,2) = -u * v * w / (0.4D1* (w * sqrt(0.2D1) - 0.2D1) ) 
             CurlBasis(30,3) = -sqrt(0.2D1) * u * w / 0.8D1

             EdgeBasis(31,1:3) = Beta(3) * Beta(4) * EdgeSign(9) * EdgeBasis(9,1:3)
             CurlBasis(31,1) = (0.2D1 * sqrt(0.2D1) * u ** 2 * w ** 2 + 0.2D1 * sqrt(0.2D1) * u * v * w ** 2 -&
                 0.2D1 * sqrt(0.2D1) * w ** 4 + 0.6D1 * sqrt(0.2D1) * u ** 2 * v - &
                 0.11D2 * sqrt(0.2D1) * v * w ** 2 - 0.8D1 * u ** 2 * v * w + 0.4D1 * v * w ** 3 + &
                 0.2D1 * sqrt(0.2D1) * u ** 2 - 0.15D2 * sqrt(0.2D1) * w ** 2 - 0.6D1 * u ** 2 * w - &
                 0.4D1 * u * v * w + 0.13D2 * w ** 3 - 0.6D1 * v * sqrt(0.2D1) + 0.20D2 * w * v - &
                 0.2D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             CurlBasis(31,2) = -(0.2D1 * sqrt(0.2D1) * u * v * w ** 2 + 0.2D1 * sqrt(0.2D1) * v ** 2 * w**2 - &
                 0.2D1 * sqrt(0.2D1) * w ** 4 + 0.6D1 * sqrt(0.2D1) * u * v ** 2 - &
                 0.11D2 * sqrt(0.2D1) * u * w ** 2 - 0.8D1 * u * v ** 2 * w + 0.4D1 * u * w ** 3 + &
                 0.2D1 * sqrt(0.2D1) * v ** 2 - 0.15D2 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u * v * w - &
                 0.6D1 * v ** 2 * w + 0.13D2 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) + 0.20D2 * u *w - &
                 0.2D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             CurlBasis(31,3) = -(u - v) * w * sqrt(0.2D1) / 0.16D2

           ELSE
             !-----------------------------------------------------------------------------------------
             ! The lowest-order pyramid from the optimal family. Now these basis functions are 
             ! also contained in the set of hierarchic basis functions, so this branch could be
             ! removed by making some code modifications (to do?).
             !-----------------------------------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(1,1) = (v*(-1 + (2*v)/(2 - Sqrt(2.0d0)*w)))/4.0d0
             EdgeBasis(1,2) = 0.0d0
             EdgeBasis(1,3) = (u*v*(-Sqrt(2.0d0) + Sqrt(2.0d0)*v + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(1,1) = (u*(-Sqrt(2.0d0) + 2*Sqrt(2.0d0)*v + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(1,2) = (v*(Sqrt(2.0d0) - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(1,3) = (-2 + 4*v + Sqrt(2.0d0)*w)/(-8 + 4*Sqrt(2.0d0)*w)
             IF (nj<ni) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,:) = -CurlBasis(1,:)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(2,1) = 0.0d0
             EdgeBasis(2,2) = (u*(1 + (2*u)/(2 - Sqrt(2.0d0)*w)))/4.0d0
             EdgeBasis(2,3) = (u*v*(Sqrt(2.0d0) + Sqrt(2.0d0)*u - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(2,1) = (u*(Sqrt(2.0d0) - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(2,2) = -(v*(Sqrt(2.0d0) + 2*Sqrt(2.0d0)*u - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(2,3) = (2 + 4*u - Sqrt(2.0d0)*w)/(8 - 4*Sqrt(2.0d0)*w)
             IF (nj<ni) THEN
               EdgeBasis(2,:) = -EdgeBasis(2,:)
               CurlBasis(2,:) = -CurlBasis(2,:)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(3,1) = (v*(1 + (2*v)/(2 - Sqrt(2.0d0)*w)))/4.0d0
             EdgeBasis(3,2) = 0.0d0
             EdgeBasis(3,3) = (u*v*(Sqrt(2.0d0) + Sqrt(2.0d0)*v - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(3,1) = (u*(Sqrt(2.0d0) + 2*Sqrt(2.0d0)*v - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(3,2) = (v*(-Sqrt(2.0d0) + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(3,3) = (2 + 4*v - Sqrt(2.0d0)*w)/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             IF (nj<ni) THEN
               EdgeBasis(3,:) = -EdgeBasis(3,:)
               CurlBasis(3,:) = -CurlBasis(3,:)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(4,1) = 0.0d0
             EdgeBasis(4,2) = (u*(-1 + (2*u)/(2 - Sqrt(2.0d0)*w)))/4.0d0
             EdgeBasis(4,3) = (u*v*(-Sqrt(2.0d0) + Sqrt(2.0d0)*u + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(4,1) = (u*(-Sqrt(2.0d0) + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(4,2) = -(v*(-Sqrt(2.0d0) + 2*Sqrt(2.0d0)*u + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(4,3) = (2 - 4*u - Sqrt(2.0d0)*w)/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             IF (nj<ni) THEN
               EdgeBasis(4,:) = -EdgeBasis(4,:)
               CurlBasis(4,:) = -CurlBasis(4,:)
             END IF

             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(5,1) = (w*(-Sqrt(2.0d0) + Sqrt(2.0d0)*v + w))/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(5,2) = (w*(-Sqrt(2.0d0) + Sqrt(2.0d0)*u + w))/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(5,3) = (u*(-2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - 2*w) + 4*w - Sqrt(2.0d0)*w**2) - &
                 (-1 + v)*(2*Sqrt(2.0d0) - 4*w + Sqrt(2.0d0)*w**2))/(4.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(5,1) = (-2*Sqrt(2.0d0) + 2*u*(Sqrt(2.0d0) - w) + 4*w - Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(5,2) = (2*Sqrt(2.0d0) - 2*Sqrt(2.0d0)*v - 4*w + 2*v*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(5,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(5,:) = -EdgeBasis(5,:)
               CurlBasis(5,:) = -CurlBasis(5,:)
             END IF

             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(6,1) = (w*(-Sqrt(2.0d0) + Sqrt(2.0d0)*v + w))/(8.0d0 - 4*Sqrt(2.0d0)*w)
             EdgeBasis(6,2) = (w*(-Sqrt(2.0d0) - Sqrt(2.0d0)*u + w))/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(6,3) = (-((-1 + v)*(2*Sqrt(2.0d0) - 4*w + Sqrt(2.0d0)*w**2)) + & 
                 u*(2*Sqrt(2.0d0) - 2*Sqrt(2.0d0)*v - 4*w + 4*v*w + Sqrt(2.0d0)*w**2))/ &
                 (4.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(6,1) = -(2*Sqrt(2.0d0) + 2*u*(Sqrt(2.0d0) - w) - 4*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(6,2) = (-2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - w) + 4*w - Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2) 
             CurlBasis(6,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(6,:) = -EdgeBasis(6,:)
               CurlBasis(6,:) = -CurlBasis(6,:)
             END IF

             i = EdgeMap(7,1)
             j = EdgeMap(7,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(7,1) = ((Sqrt(2.0d0) + Sqrt(2.0d0)*v - w)*w)/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(7,2) = ((Sqrt(2.0d0) + Sqrt(2.0d0)*u - w)*w)/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(7,3) = ((1 + v)*(2*Sqrt(2.0d0) - 4*w + Sqrt(2.0d0)*w**2) + &
                 u*(2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - 2*w) - 4*w + Sqrt(2.0d0)*w**2))/ &
                 (4.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(7,1) = (2*Sqrt(2.0d0) + 2*u*(Sqrt(2.0d0) - w) - 4*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(7,2) = -(2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - w) - 4*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(7,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(7,:) = -EdgeBasis(7,:)
               CurlBasis(7,:) = -CurlBasis(7,:)
             END IF

             i = EdgeMap(8,1)
             j = EdgeMap(8,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(8,1) = (w*(-Sqrt(2.0d0) - Sqrt(2.0d0)*v + w))/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(8,2) = (w*(-Sqrt(2.0d0) + Sqrt(2.0d0)*u + w))/(8.0d0 - 4*Sqrt(2.0d0)*w)
             EdgeBasis(8,3) = ((1 + v)*(2*Sqrt(2.0d0) - 4*w + Sqrt(2.0d0)*w**2) - &
                 u*(2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - 2*w) - 4*w + Sqrt(2.0d0)*w**2))/ &
                 (4.0d0*(-2.0d0 + Sqrt(2.0d0)*w)**2)
             CurlBasis(8,1) = (2*Sqrt(2.0d0) - 2*Sqrt(2.0d0)*u - 4*w + 2*u*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2.0d0 + Sqrt(2.0d0)*w)**2)
             CurlBasis(8,2) = (2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - w) - 4*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2.0d0 + Sqrt(2.0d0)*w)**2)
             CurlBasis(8,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(8,:) = -EdgeBasis(8,:)
               CurlBasis(8,:) = -CurlBasis(8,:)
             END IF

             ! ------------------------------------------------------------------
             ! The last two basis functions are associated with the square face.
             ! We first create the basis function in the default order without
             ! sign reversions.
             ! ------------------------------------------------------------------
             SquareFaceMap(:) = (/ 1,2,3,4 /)
             Ind => Element % Nodeindexes          

             WorkBasis(1,1) = (2.0d0 - 2*v**2 - 2*Sqrt(2.0d0)*w + w**2)/(4.0d0 - 2*Sqrt(2.0d0)*w)
             WorkBasis(1,2) = 0.0d0
             WorkBasis(1,3) = (u*(1.0d0 - (4*v**2)/(-2.0d0 + Sqrt(2.0d0)*w)**2))/(2.0d0*Sqrt(2.0d0))
             WorkCurlBasis(1,1) = (-2*Sqrt(2.0d0)*u*v)/(-2.0d0 + Sqrt(2.0d0)*w)**2
             WorkCurlBasis(1,2) = (-2*Sqrt(2.0d0) + 4*w - Sqrt(2.0d0)*w**2)/(-2.0d0 + Sqrt(2.0d0)*w)**2
             WorkCurlBasis(1,3) = (2.0d0*v)/(2.0d0 - Sqrt(2.0d0)*w)

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = (2.0d0 - 2*u**2 - 2*Sqrt(2.0d0)*w + w**2)/(4.0d0 - 2*Sqrt(2.0d0)*w)
             WorkBasis(2,3) = (v*(1.0d0 - (4*u**2)/(-2.0d0 + Sqrt(2.0d0)*w)**2))/(2.0d0*Sqrt(2.0d0))
             WorkCurlBasis(2,1) = (2*Sqrt(2.0d0) - 4*w + Sqrt(2.0d0)*w**2)/(-2.0d0 + Sqrt(2.0d0)*w)**2
             WorkCurlBasis(2,2) = (2*Sqrt(2.0d0)*u*v)/(-2.0d0 + Sqrt(2.0d0)*w)**2
             WorkCurlBasis(2,3) = (2*u)/(-2.0d0 + Sqrt(2.0d0)*w)

             ! -------------------------------------------------------------------
             ! Finally apply an order change and sign reversions if needed. 
             ! -------------------------------------------------------------------
             DO j=1,4
               FaceIndeces(j) = Ind(SquareFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(9,:) = D1 * WorkBasis(I1,:)
             CurlBasis(9,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(10,:) = D2 * WorkBasis(I2,:)
             CurlBasis(10,:) = D2 * WorkCurlBasis(I2,:)          
           END IF

         CASE(7)
           !--------------------------------------------------------------
           ! This branch is for handling prismatic (or wedge) elements
           !--------------------------------------------------------------           
           EdgeMap => LGetEdgeMap(7)
           Ind => Element % Nodeindexes

           IF (SecondOrder) THEN
             !---------------------------------------------------------------
             ! The second-order element from the Nedelec's first family 
             ! (note that the lowest-order prism element is from a different 
             ! family). This element may not be optimally accurate if 
             ! the physical element is not affine.
             !--------------------------------------------------------------             
             h1 = 0.5d0 * (1-w)
             dh1 = -0.5d0
             h2 = 0.5d0 * (1+w)
             dh2 = 0.5d0
             h3 = h1 * h2
             dh3 = -0.5d0 * w

             ! ---------------------------------------------------------
             ! The first and fourth edges ...
             !--------------------------------------------------------
             ! The corresponding basis functions for the triangle:
             !--------------------------------------------------------
             WorkBasis(1,1) = (3.0d0 - Sqrt(3.0d0)*v)/6.0d0
             WorkBasis(1,2) = u/(2.0d0*Sqrt(3.0d0))
             WorkCurlBasis(1,3) = 1.0d0/Sqrt(3.0d0)
             WorkBasis(2,1) = -(u*(-3.0d0 + Sqrt(3.0d0)*v))/2.0d0
             WorkBasis(2,2) = (Sqrt(3.0d0)*u**2)/2.0d0
             WorkCurlBasis(2,3) = (3.0d0*Sqrt(3.0d0)*u)/2.0d0

             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(1,1:2) = WorkBasis(1,1:2) * h1
             CurlBasis(1,1) = -WorkBasis(1,2) * dh1
             CurlBasis(1,2) = WorkBasis(1,1) * dh1
             CurlBasis(1,3) = WorkCurlBasis(1,3) * h1
             EdgeBasis(2,1:2) = WorkBasis(2,1:2) * h1
             CurlBasis(2,1) = -WorkBasis(2,2) * dh1
             CurlBasis(2,2) = WorkBasis(2,1) * dh1
             CurlBasis(2,3) = WorkCurlBasis(2,3) * h1
             IF (nj<ni) THEN
               EdgeBasis(1,1:2) = -EdgeBasis(1,1:2)
               CurlBasis(1,1:3) = -CurlBasis(1,1:3)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(7,1:2) = WorkBasis(1,1:2) * h2
             CurlBasis(7,1) = -WorkBasis(1,2) * dh2
             CurlBasis(7,2) = WorkBasis(1,1) * dh2
             CurlBasis(7,3) = WorkCurlBasis(1,3) * h2
             EdgeBasis(8,1:2) = WorkBasis(2,1:2) * h2
             CurlBasis(8,1) = -WorkBasis(2,2) * dh2
             CurlBasis(8,2) = WorkBasis(2,1) * dh2
             CurlBasis(8,3) = WorkCurlBasis(2,3) * h2
             IF (nj<ni) THEN
               EdgeBasis(7,1:2) = -EdgeBasis(7,1:2)
               CurlBasis(7,1:3) = -CurlBasis(7,1:3)
             END IF

             ! ---------------------------------------------------------
             ! The second and fifth edges ...
             !--------------------------------------------------------
             ! The corresponding basis functions for the triangle:
             !--------------------------------------------------------
             WorkBasis(1,1) = -v/(2.0d0*Sqrt(3.0d0))
             WorkBasis(1,2) = (1 + u)/(2.0d0*Sqrt(3.0d0))
             WorkCurlBasis(1,3) = 1.0d0/Sqrt(3.0d0)
             WorkBasis(2,1) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*v)/4.0d0
             WorkBasis(2,2) = (Sqrt(3.0d0)*(1.0d0 + u)*(-1.0d0 - u + Sqrt(3.0d0)*v))/4.0d0
             WorkCurlBasis(2,3) = (-3.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v))/4.0d0

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(3,1:2) = WorkBasis(1,1:2) * h1
             CurlBasis(3,1) = -WorkBasis(1,2) * dh1
             CurlBasis(3,2) = WorkBasis(1,1) * dh1
             CurlBasis(3,3) = WorkCurlBasis(1,3) * h1
             EdgeBasis(4,1:2) = WorkBasis(2,1:2) * h1
             CurlBasis(4,1) = -WorkBasis(2,2) * dh1
             CurlBasis(4,2) = WorkBasis(2,1) * dh1
             CurlBasis(4,3) = WorkCurlBasis(2,3) * h1
             IF (nj<ni) THEN
               EdgeBasis(3,1:2) = -EdgeBasis(3,1:2)
               CurlBasis(3,1:3) = -CurlBasis(3,1:3)
             END IF

             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(9,1:2) = WorkBasis(1,1:2) * h2
             CurlBasis(9,1) = -WorkBasis(1,2) * dh2
             CurlBasis(9,2) = WorkBasis(1,1) * dh2
             CurlBasis(9,3) = WorkCurlBasis(1,3) * h2
             EdgeBasis(10,1:2) = WorkBasis(2,1:2) * h2
             CurlBasis(10,1) = -WorkBasis(2,2) * dh2
             CurlBasis(10,2) = WorkBasis(2,1) * dh2
             CurlBasis(10,3) = WorkCurlBasis(2,3) * h2
             IF (nj<ni) THEN
               EdgeBasis(9,1:2) = -EdgeBasis(9,1:2)
               CurlBasis(9,1:3) = -CurlBasis(9,1:3)
             END IF

             ! ---------------------------------------------------------
             ! The third and sixth edges ...
             !--------------------------------------------------------
             ! The corresponding basis functions for the triangle:
             !--------------------------------------------------------
             WorkBasis(1,1) = -v/(2.0d0*Sqrt(3.0d0))
             WorkBasis(1,2) = (-1 + u)/(2.0d0*Sqrt(3.0d0))
             WorkCurlBasis(1,3) =  1.0d0/Sqrt(3.0d0)
             WorkBasis(2,1) = (v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0
             WorkBasis(2,2) = -(Sqrt(3.0d0)*(-1.0d0 + u)*(-1.0d0 + u + Sqrt(3.0d0)*v))/4.0d0
             WorkCurlBasis(2,3) = (-3.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(5,1:2) = WorkBasis(1,1:2) * h1
             CurlBasis(5,1) = -WorkBasis(1,2) * dh1
             CurlBasis(5,2) = WorkBasis(1,1) * dh1
             CurlBasis(5,3) = WorkCurlBasis(1,3) * h1
             EdgeBasis(6,1:2) = WorkBasis(2,1:2) * h1
             CurlBasis(6,1) = -WorkBasis(2,2) * dh1
             CurlBasis(6,2) = WorkBasis(2,1) * dh1
             CurlBasis(6,3) = WorkCurlBasis(2,3) * h1
             IF (nj<ni) THEN
               EdgeBasis(5,1:2) = -EdgeBasis(5,1:2)
               CurlBasis(5,1:3) = -CurlBasis(5,1:3)
             END IF

             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(11,1:2) = WorkBasis(1,1:2) * h2
             CurlBasis(11,1) = -WorkBasis(1,2) * dh2
             CurlBasis(11,2) = WorkBasis(1,1) * dh2
             CurlBasis(11,3) = WorkCurlBasis(1,3) * h2
             EdgeBasis(12,1:2) = WorkBasis(2,1:2) * h2
             CurlBasis(12,1) = -WorkBasis(2,2) * dh2
             CurlBasis(12,2) = WorkBasis(2,1) * dh2
             CurlBasis(12,3) = WorkCurlBasis(2,3) * h2
             IF (nj<ni) THEN
               EdgeBasis(11,1:2) = -EdgeBasis(11,1:2)
               CurlBasis(11,1:3) = -CurlBasis(11,1:3)
             END IF

             ! -------------------------------------------------------
             ! The edges 14, 25 and 36
             !--------------------------------------------------------
             DO q = 1,3
               i = EdgeMap(6+q,1)
               j = EdgeMap(6+q,2)
               ni = Element % NodeIndexes(i)
               IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
               nj = Element % NodeIndexes(j)
               IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)

               grad(1:2) = dTriangleNodalPBasis(q, u, v)
               EdgeBasis(12+(q-1)*2+1,3) = 0.5d0 * TriangleNodalPBasis(q, u, v)
               CurlBasis(12+(q-1)*2+1,1) = 0.5d0* grad(2)
               CurlBasis(12+(q-1)*2+1,2) = -0.5d0* grad(1)
               EdgeBasis(12+(q-1)*2+2,3) = 3.0d0 * EdgeBasis(12+(q-1)*2+1,3) * w
               CurlBasis(12+(q-1)*2+2,1) = 1.5d0 * grad(2) * w
               CurlBasis(12+(q-1)*2+2,2) = -1.5d0 * grad(1) * w

               IF (nj<ni) THEN
                 EdgeBasis(12+(q-1)*2+1,3) = -EdgeBasis(12+(q-1)*2+1,3)
                 CurlBasis(12+(q-1)*2+1,1:2) = -CurlBasis(12+(q-1)*2+1,1:2)
               END IF
             END DO

             !-------------------------------------------------
             ! Two basis functions defined on the face 123:
             !-------------------------------------------------
             TriangleFaceMap(:) = (/ 1,2,3 /)

             DO j=1,3
               FaceIndeces(j) = Ind(TriangleFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,3
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             WorkBasis(1,1) = ((Sqrt(3.0d0) - v)*v)/6.0d0
             WorkBasis(1,2) = (u*v)/6.0d0
             WorkCurlBasis(1,3) = (-Sqrt(3.0d0) + 3.0d0*v)/6.0d0
             WorkBasis(2,1) = (v*(1.0d0 + u - v/Sqrt(3.0d0)))/(4.0d0*Sqrt(3.0d0))
             WorkBasis(2,2) = ((-1.0d0 + u)*(-3.0d0 - 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
             WorkCurlBasis(2,3) =(-Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u + 3.0d0*v)/12.0d0
             WorkBasis(3,1) = (v*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
             WorkBasis(3,2) = -((1.0d0 + u)*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
             WorkCurlBasis(3,3) = (Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u - 3.0d0*v)/12.0d0

             EdgeBasis(19,1:2) = D1 * WorkBasis(I1,1:2) * h1
             CurlBasis(19,1) = -D1 * WorkBasis(I1,2) * dh1
             CurlBasis(19,2) = D1 * WorkBasis(I1,1) * dh1
             CurlBasis(19,3) = D1 * WorkCurlBasis(I1,3) * h1

             EdgeBasis(20,1:2) = D2 * WorkBasis(I2,1:2) * h1
             CurlBasis(20,1) = -D2 * WorkBasis(I2,2) * dh1
             CurlBasis(20,2) = D2 * WorkBasis(I2,1) * dh1
             CurlBasis(20,3) = D2 * WorkCurlBasis(I2,3) * h1

             !-------------------------------------------------
             ! Two basis functions defined on the face 456:
             !-------------------------------------------------
             TriangleFaceMap(:) = (/ 4,5,6 /)

             DO j=1,3
               FaceIndeces(j) = Ind(TriangleFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,3
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(21,1:2) = D1 * WorkBasis(I1,1:2) * h2
             CurlBasis(21,1) = -D1 * WorkBasis(I1,2) * dh2
             CurlBasis(21,2) = D1 * WorkBasis(I1,1) * dh2
             CurlBasis(21,3) = D1 * WorkCurlBasis(I1,3) * h2

             EdgeBasis(22,1:2) = D2 * WorkBasis(I2,1:2) * h2
             CurlBasis(22,1) = -D2 * WorkBasis(I2,2) * dh2
             CurlBasis(22,2) = D2 * WorkBasis(I2,1) * dh2
             CurlBasis(22,3) = D2 * WorkCurlBasis(I2,3) * h2

             !-------------------------------------------------
             ! Four basis functions defined on the face 1254:
             !-------------------------------------------------              
             SquareFaceMap(:) = (/ 1,2,5,4 /)          
             WorkBasis = 0.0d0
             WorkCurlBasis = 0.0d0

             WorkBasis(1,1) = (3.0d0 - Sqrt(3.0d0)*v)/6.0d0 * 4.0d0 * h3
             WorkBasis(1,2) = u/(2.0d0*Sqrt(3.0d0)) * 4.0d0 * h3
             WorkCurlBasis(1,1) = -WorkBasis(1,2)/h3 * dh3 
             WorkCurlBasis(1,2) = WorkBasis(1,1)/h3 * dh3 
             WorkCurlBasis(1,3) = 1.0d0/Sqrt(3.0d0) * 4.0d0 * h3
             WorkBasis(2,1) = -(u*(-3.0d0 + Sqrt(3.0d0)*v))/2.0d0 * 4.0d0 * h3
             WorkBasis(2,2) = (Sqrt(3.0d0)*u**2)/2.0d0 * 4.0d0 * h3
             WorkCurlBasis(2,1) = -WorkBasis(2,2)/h3 * dh3 
             WorkCurlBasis(2,2) = WorkBasis(2,1)/h3 * dh3
             WorkCurlBasis(2,3) = (3.0d0*Sqrt(3.0d0)*u)/2.0d0 * 4.0d0 * h3

             WorkBasis(3,3) = 2.0d0 * TriangleNodalPBasis(1, u, v) * TriangleNodalPBasis(2, u, v)
             grad(1:2) = dTriangleNodalPBasis(1, u, v) * TriangleNodalPBasis(2, u, v) + &
                 TriangleNodalPBasis(1, u, v) * dTriangleNodalPBasis(2, u, v)
             WorkCurlBasis(3,1) = 2.0d0 * grad(2)
             WorkCurlBasis(3,2) = -2.0d0 * grad(1)
             WorkBasis(4,3) = 3.0d0 * WorkBasis(3,3) * w
             WorkCurlBasis(4,1) = 6.0d0 * grad(2) * w
             WorkCurlBasis(4,2) = -6.0d0 * grad(1) * w

             DO j=1,4
               FaceIndeces(j) = Ind(SquareFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(23,:) = D1 * WorkBasis(2*(I1-1)+1,:)
             CurlBasis(23,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
             EdgeBasis(24,:) = WorkBasis(2*(I1-1)+2,:)
             CurlBasis(24,:) = WorkCurlBasis(2*(I1-1)+2,:)
             EdgeBasis(25,:) = D2 * WorkBasis(2*(I2-1)+1,:)
             CurlBasis(25,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
             EdgeBasis(26,:) = WorkBasis(2*(I2-1)+2,:)
             CurlBasis(26,:) = WorkCurlBasis(2*(I2-1)+2,:)            

             !-------------------------------------------------
             ! Four basis functions defined on the face 2365:
             !-------------------------------------------------              
             SquareFaceMap(:) = (/ 2,3,6,5 /)          
             WorkBasis = 0.0d0
             WorkCurlBasis = 0.0d0

             WorkBasis(1,1) = -v/(2.0d0*Sqrt(3.0d0)) * 4.0d0 * h3
             WorkBasis(1,2) = (1 + u)/(2.0d0*Sqrt(3.0d0)) * 4.0d0 * h3
             WorkCurlBasis(1,1) = -WorkBasis(1,2)/h3 * dh3 
             WorkCurlBasis(1,2) = WorkBasis(1,1)/h3 * dh3 
             WorkCurlBasis(1,3) = 1.0d0/Sqrt(3.0d0) * 4.0d0 * h3
             WorkBasis(2,1) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*v)/4.0d0 * 4.0d0 * h3
             WorkBasis(2,2) = (Sqrt(3.0d0)*(1.0d0 + u)*(-1.0d0 - u + Sqrt(3.0d0)*v))/4.0d0 * 4.0d0 * h3
             WorkCurlBasis(2,1) = -WorkBasis(2,2)/h3 * dh3 
             WorkCurlBasis(2,2) = WorkBasis(2,1)/h3 * dh3
             WorkCurlBasis(2,3) = (-3.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v))/4.0d0 * 4.0d0 * h3

             WorkBasis(3,3) = 2.0d0 * TriangleNodalPBasis(2, u, v) * TriangleNodalPBasis(3, u, v)
             grad(1:2) = dTriangleNodalPBasis(2, u, v) * TriangleNodalPBasis(3, u, v) + &
                 TriangleNodalPBasis(2, u, v) * dTriangleNodalPBasis(3, u, v)
             WorkCurlBasis(3,1) = 2.0d0 * grad(2)
             WorkCurlBasis(3,2) = -2.0d0 * grad(1)
             WorkBasis(4,3) = 3.0d0 * WorkBasis(3,3) * w
             WorkCurlBasis(4,1) = 6.0d0 * grad(2) * w
             WorkCurlBasis(4,2) = -6.0d0 * grad(1) * w

             DO j=1,4
               FaceIndeces(j) = Ind(SquareFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(27,:) = D1 * WorkBasis(2*(I1-1)+1,:)
             CurlBasis(27,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
             EdgeBasis(28,:) = WorkBasis(2*(I1-1)+2,:)
             CurlBasis(28,:) = WorkCurlBasis(2*(I1-1)+2,:)
             EdgeBasis(29,:) = D2 * WorkBasis(2*(I2-1)+1,:)
             CurlBasis(29,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
             EdgeBasis(30,:) = WorkBasis(2*(I2-1)+2,:)
             CurlBasis(30,:) = WorkCurlBasis(2*(I2-1)+2,:)  

             !-------------------------------------------------
             ! Four basis functions defined on the face 3146:
             !-------------------------------------------------              
             SquareFaceMap(:) = (/ 3,1,4,6 /)          
             WorkBasis = 0.0d0
             WorkCurlBasis = 0.0d0

             WorkBasis(1,1) = -v/(2.0d0*Sqrt(3.0d0)) * 4.0d0 * h3
             WorkBasis(1,2) = (-1 + u)/(2.0d0*Sqrt(3.0d0)) * 4.0d0 * h3
             WorkCurlBasis(1,1) = -WorkBasis(1,2)/h3 * dh3 
             WorkCurlBasis(1,2) = WorkBasis(1,1)/h3 * dh3 
             WorkCurlBasis(1,3) = 1.0d0/Sqrt(3.0d0) * 4.0d0 * h3
             WorkBasis(2,1) = (v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0 * 4.0d0 * h3
             WorkBasis(2,2) =  -(Sqrt(3.0d0)*(-1.0d0 + u)*(-1.0d0 + u + Sqrt(3.0d0)*v))/4.0d0 * 4.0d0 * h3
             WorkCurlBasis(2,1) = -WorkBasis(2,2)/h3 * dh3 
             WorkCurlBasis(2,2) = WorkBasis(2,1)/h3 * dh3
             WorkCurlBasis(2,3) = (-3.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0 * 4.0d0 * h3

             WorkBasis(3,3) = 2.0d0 * TriangleNodalPBasis(3, u, v) * TriangleNodalPBasis(1, u, v)
             grad(1:2) = dTriangleNodalPBasis(3, u, v) * TriangleNodalPBasis(1, u, v) + &
                 TriangleNodalPBasis(3, u, v) * dTriangleNodalPBasis(1, u, v)
             WorkCurlBasis(3,1) = 2.0d0 * grad(2)
             WorkCurlBasis(3,2) = -2.0d0 * grad(1)
             WorkBasis(4,3) = 3.0d0 * WorkBasis(3,3) * w
             WorkCurlBasis(4,1) = 6.0d0 * grad(2) * w
             WorkCurlBasis(4,2) = -6.0d0 * grad(1) * w

             DO j=1,4
               FaceIndeces(j) = Ind(SquareFaceMap(j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(31,:) = D1 * WorkBasis(2*(I1-1)+1,:)
             CurlBasis(31,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
             EdgeBasis(32,:) = WorkBasis(2*(I1-1)+2,:)
             CurlBasis(32,:) = WorkCurlBasis(2*(I1-1)+2,:)
             EdgeBasis(33,:) = D2 * WorkBasis(2*(I2-1)+1,:)
             CurlBasis(33,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
             EdgeBasis(34,:) = WorkBasis(2*(I2-1)+2,:)
             CurlBasis(34,:) = WorkCurlBasis(2*(I2-1)+2,:)  

             !-------------------------------------------------
             ! Two basis functions associated with the interior
             !-------------------------------------------------    
             EdgeBasis(35,1) = (v*(1.0d0 + u - v/Sqrt(3.0d0)))/(4.0d0*Sqrt(3.0d0)) * h3
             EdgeBasis(35,2) = ((-1.0d0 + u)*(-3.0d0 - 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0)) * h3
             CurlBasis(35,1) = -EdgeBasis(35,2)/h3 * dh3
             CurlBasis(35,2) = EdgeBasis(35,1)/h3 * dh3
             CurlBasis(35,3) = (-Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u + 3.0d0*v)/12.0d0 * h3

             EdgeBasis(36,1) = (v*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0)) * h3
             EdgeBasis(36,2) = -((1.0d0 + u)*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0)) * h3
             CurlBasis(36,1) = -EdgeBasis(36,2)/h3 * dh3
             CurlBasis(36,2) = EdgeBasis(36,1)/h3 * dh3
             CurlBasis(36,3) = (Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u - 3.0d0*v)/12.0d0 * h3

           ELSE
             !--------------------------------------------------------------
             ! The lowest-order element from the optimal family. The optimal
             ! accuracy is obtained also for non-affine meshes.
             ! -------------------------------------------------------------
             ! First nine basis functions associated with the edges
             ! -------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(1,1) = -((-3.0d0 + Sqrt(3.0d0)*v)*(-1.0d0 + w)*w)/12.0d0
             EdgeBasis(1,2) = (u*(-1.0d0 + w)*w)/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(1,3) = 0.0d0
             CurlBasis(1,1) = (u*(1.0d0 - 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(1,2) = -((-3.0d0 + Sqrt(3.0d0)*v)*(-1.0d0 + 2*w))/12.0d0
             CurlBasis(1,3) = ((-1.0d0 + w)*w)/(2.0d0*Sqrt(3.0d0))
             IF (nj<ni) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,:) = -CurlBasis(1,:)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(2,1) = -(v*(-1.0d0 + w)*w)/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(2,2) = ((1.0d0 + u)*(-1.0d0 + w)*w)/(4.0d0*Sqrt(3.0d0)) 
             EdgeBasis(2,3) = 0.0d0
             CurlBasis(2,1) = ((1.0d0 + u)*(1.0d0 - 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(2,2) = (v*(1.0d0 - 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(2,3) = ((-1.0d0 + w)*w)/(2.0d0*Sqrt(3.0d0))
             IF (nj<ni) THEN
               EdgeBasis(2,:) = -EdgeBasis(2,:)
               CurlBasis(2,:) = -CurlBasis(2,:)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(3,1) = -(v*(-1.0d0 + w)*w)/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(3,2) = ((-1.0d0 + u)*(-1.0d0 + w)*w)/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(3,3) = 0.0d0
             CurlBasis(3,1) = ((-1.0d0 + u)*(1.0d0 - 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(3,2) = (v*(1.0d0 - 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(3,3) = ((-1.0d0 + w)*w)/(2.0d0*Sqrt(3.0d0))
             IF (nj<ni) THEN
               EdgeBasis(3,:) = -EdgeBasis(3,:)
               CurlBasis(3,:) = -CurlBasis(3,:)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(4,1) = -((-3.0d0 + Sqrt(3.0d0)*v)*w*(1.0d0 + w))/12.0d0
             EdgeBasis(4,2) = (u*w*(1.0d0 + w))/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(4,3) = 0.0d0
             CurlBasis(4,1) = -(u*(1.0d0 + 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(4,2) = -((-3.0d0 + Sqrt(3.0d0)*v)*(1.0d0 + 2.0d0*w))/12.0d0
             CurlBasis(4,3) = (w*(1.0d0 + w))/(2.0d0*Sqrt(3.0d0))
             IF (nj<ni) THEN
               EdgeBasis(4,:) = -EdgeBasis(4,:)
               CurlBasis(4,:) = -CurlBasis(4,:)
             END IF

             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(5,1) = -(v*w*(1.0d0 + w))/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(5,2) = ((1.0d0 + u)*w*(1.0d0 + w))/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(5,3) = 0.0d0
             CurlBasis(5,1) = -((1.0d0 + u)*(1.0d0 + 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(5,2) = -(v*(1.0d0 + 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(5,3) = (w*(1.0d0 + w))/(2.0d0*Sqrt(3.0d0))
             IF (nj<ni) THEN
               EdgeBasis(5,:) = -EdgeBasis(5,:)
               CurlBasis(5,:) = -CurlBasis(5,:)
             END IF

             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(6,1) = -(v*w*(1.0d0 + w))/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(6,2) = ((-1.0d0 + u)*w*(1.0d0 + w))/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(6,3) = 0.0d0
             CurlBasis(6,1) = -((-1.0d0 + u)*(1.0d0 + 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(6,2) = -(v*(1.0d0 + 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(6,3) = (w*(1.0d0 + w))/(2.0d0*Sqrt(3.0d0))
             IF (nj<ni) THEN
               EdgeBasis(6,:) = -EdgeBasis(6,:)
               CurlBasis(6,:) = -CurlBasis(6,:)
             END IF

             i = EdgeMap(7,1)
             j = EdgeMap(7,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(7,1) = 0.0d0
             EdgeBasis(7,2) = 0.0d0
             EdgeBasis(7,3) = (3*u**2 + v*(-Sqrt(3.0d0) + v) + u*(-3.0d0 + 2*Sqrt(3.0d0)*v))/12.0d0
             CurlBasis(7,1) = (-Sqrt(3.0d0) + 2*Sqrt(3.0d0)*u + 2*v)/12.0d0
             CurlBasis(7,2) = (3.0d0 - 6*u - 2*Sqrt(3.0d0)*v)/12.0d0
             CurlBasis(7,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(7,:) = -EdgeBasis(7,:)
               CurlBasis(7,:) = -CurlBasis(7,:)
             END IF

             i = EdgeMap(8,1)
             j = EdgeMap(8,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(8,1) = 0.0d0
             EdgeBasis(8,2) = 0.0d0
             EdgeBasis(8,3) = (3*u**2 + v*(-Sqrt(3.0d0) + v) + u*(3.0d0 - 2*Sqrt(3.0d0)*v))/12.0d0
             CurlBasis(8,1) = (-Sqrt(3.0d0) - 2*Sqrt(3.0d0)*u + 2*v)/12.0d0
             CurlBasis(8,2) = (-3.0d0 - 6*u + 2*Sqrt(3.0d0)*v)/12.0d0 
             CurlBasis(8,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(8,:) = -EdgeBasis(8,:)
               CurlBasis(8,:) = -CurlBasis(8,:)
             END IF

             i = EdgeMap(9,1)
             j = EdgeMap(9,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(9,1) = 0.0d0
             EdgeBasis(9,2) = 0.0d0
             EdgeBasis(9,3) = (v*(-Sqrt(3.0d0) + 2*v))/6.0d0
             CurlBasis(9,1) = (-Sqrt(3.0d0) + 4*v)/6.0d0
             CurlBasis(9,2) = 0.0d0
             CurlBasis(9,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(9,:) = -EdgeBasis(9,:)
               CurlBasis(9,:) = -CurlBasis(9,:)
             END IF

             ! ---------------------------------------------------------------------
             ! Additional six basis functions on the square faces (two per face).
             ! ---------------------------------------------------------------------         
             PrismSquareFaceMap(1,:) = (/ 1,2,5,4 /)
             PrismSquareFaceMap(2,:) = (/ 2,3,6,5 /)
             PrismSquareFaceMap(3,:) = (/ 3,1,4,6 /)

             ! The first square face:
             WorkBasis(1,1) = ((-3.0d0 + Sqrt(3.0d0)*v)*(-1.0d0 + w**2))/6.0d0
             WorkBasis(1,2) = -(u*(-1.0d0 + w**2))/(2.0d0*Sqrt(3.0d0))
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = (u*w)/Sqrt(3.0d0)
             WorkCurlBasis(1,2) = (-1.0d0 + v/Sqrt(3.0d0))*w
             WorkCurlBasis(1,3) = -((-1.0d0 + w**2)/Sqrt(3.0d0)) 

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = (3.0d0 - 3*u**2 - 2*Sqrt(3.0d0)*v + v**2)/6.0d0
             WorkCurlBasis(2,1) = (-Sqrt(3.0d0) + v)/3.0d0
             WorkCurlBasis(2,2) = u
             WorkCurlBasis(2,3) = 0.0d0

             DO j=1,4
               FaceIndeces(j) = Ind(PrismSquareFaceMap(1,j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(10,:) = D1 * WorkBasis(I1,:)
             CurlBasis(10,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(11,:) = D2 * WorkBasis(I2,:)
             CurlBasis(11,:) = D2 * WorkCurlBasis(I2,:) 

             ! The second square face:
             WorkBasis(1,1) = (v*(-1.0d0 + w**2))/(2.0d0*Sqrt(3.0d0))
             WorkBasis(1,2) = -((1.0d0 + u)*(-1.0d0 + w**2))/(2.0d0*Sqrt(3.0d0))
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = ((1.0d0 + u)*w)/Sqrt(3.0d0)
             WorkCurlBasis(1,2) = (v*w)/Sqrt(3.0d0)
             WorkCurlBasis(1,3) = -((-1.0d0 + w**2)/Sqrt(3.0d0))

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - v)*v)/3.0d0
             WorkCurlBasis(2,1) = (Sqrt(3.0d0) + Sqrt(3.0d0)*u - 2*v)/3.0d0
             WorkCurlBasis(2,2) = -(v/Sqrt(3.0d0))
             WorkCurlBasis(2,3) = 0.0d0 

             DO j=1,4
               FaceIndeces(j) = Ind(PrismSquareFaceMap(2,j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(12,:) = D1 * WorkBasis(I1,:)
             CurlBasis(12,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(13,:) = D2 * WorkBasis(I2,:)
             CurlBasis(13,:) = D2 * WorkCurlBasis(I2,:) 

             ! The third square face:
             WorkBasis(1,1) = (v*(-1.0d0 + w**2))/(2.0d0*SQRT(3.0d0))
             WorkBasis(1,2) = -((-1.0d0 + u)*(-1.0d0 + w**2))/(2.0d0*SQRT(3.0d0))
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = ((-1.0d0 + u)*w)/SQRT(3.0d0)
             WorkCurlBasis(1,2) = (v*w)/SQRT(3.0d0)
             WorkCurlBasis(1,3) = -(-1.0d0 + w**2)/SQRT(3.0d0)

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = -(v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v))/3.0d0
             WorkCurlBasis(2,1) = (Sqrt(3.0d0) - Sqrt(3.0d0)*u - 2*v)/3.0d0
             WorkCurlBasis(2,2) = v/Sqrt(3.0d0)
             WorkCurlBasis(2,3) = 0.0d0

             DO j=1,4
               FaceIndeces(j) = Ind(PrismSquareFaceMap(3,j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(14,:) = D1 * WorkBasis(I1,:)
             CurlBasis(14,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(15,:) = D2 * WorkBasis(I2,:)
             CurlBasis(15,:) = D2 * WorkCurlBasis(I2,:) 
           END IF

         CASE(8)
           !--------------------------------------------------------------
           ! This branch is for handling brick elements
           !--------------------------------------------------------------           
           EdgeMap => LGetEdgeMap(8)
           Ind => Element % Nodeindexes
           
           IF (SecondOrder) THEN
             !---------------------------------------------------------------
             ! The second-order element from the Nedelec's first family 
             ! (note that the lowest-order brick element is from a different 
             ! family). This element may not be optimally accurate if 
             ! the physical element is not affine.
             !--------------------------------------------------------------             
  
             ! Edges 12 and 43 ...
             DO q=1,2
               k = 2*q-1 ! Edge number k: 1 ~ 12 and 3 ~ 43 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               ni = Element % NodeIndexes(i)
               IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
               nj = Element % NodeIndexes(j)
               IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
 
               EdgeBasis(2*(k-1)+1,1) = 0.5d0 * LineNodalPBasis(1,w) * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,2) = 0.5d0 * (-0.5d0) * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,3) = -0.5d0 * LineNodalPBasis(1,w) * dLineNodalPBasis(q,v)
               EdgeBasis(2*(k-1)+2,1) = 1.5d0 * LineNodalPBasis(1,w) * u * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,2) = 1.5d0 * (-0.5d0) * u * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,3) = -1.5d0 * LineNodalPBasis(1,w) * u * dLineNodalPBasis(q,v)
               IF (nj<ni) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO

             ! Edges 56 and 87 ...
             DO q=1,2
               k = 4 + 2*q-1 ! Edge number k: 5 ~ 56 and 7 ~ 87 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               ni = Element % NodeIndexes(i)
               IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
               nj = Element % NodeIndexes(j)
               IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
 
               EdgeBasis(2*(k-1)+1,1) = 0.5d0 * LineNodalPBasis(2,w) * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,2) = 0.5d0 * 0.5d0 * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,3) = -0.5d0 * LineNodalPBasis(2,w) * dLineNodalPBasis(q,v)
               EdgeBasis(2*(k-1)+2,1) = 1.5d0 * LineNodalPBasis(2,w) * u * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,2) = 1.5d0 * 0.5d0 * u * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,3) = -1.5d0 * LineNodalPBasis(2,w) * u * dLineNodalPBasis(q,v)
               IF (nj<ni) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO

             ! Edges 23 and 14 ...
             DO q=1,2
               k = 2*q ! Edge number k: 2 ~ 23 and 4 ~ 14 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               ni = Element % NodeIndexes(i)
               IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
               nj = Element % NodeIndexes(j)
               IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
 
               EdgeBasis(2*(k-1)+1,2) = 0.5d0 * LineNodalPBasis(1,w) * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+1,1) = -0.5d0 * (-0.5d0) * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+1,3) = 0.5d0 * LineNodalPBasis(1,w) * dLineNodalPBasis(3-q,u)
               EdgeBasis(2*(k-1)+2,2) = 1.5d0 * LineNodalPBasis(1,w) * v * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+2,1) = -1.5d0 * (-0.5d0) * v * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+2,3) = 1.5d0 * LineNodalPBasis(1,w) * v * dLineNodalPBasis(3-q,u)
               IF (nj<ni) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO            

             ! Edges 67 and 58 ...
             DO q=1,2
               k = 4+2*q ! Edge number k: 6 ~ 67 and 8 ~ 58 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               ni = Element % NodeIndexes(i)
               IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
               nj = Element % NodeIndexes(j)
               IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
 
               EdgeBasis(2*(k-1)+1,2) = 0.5d0 * LineNodalPBasis(2,w) * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+1,1) = -0.5d0 * 0.5d0 * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+1,3) = 0.5d0 * LineNodalPBasis(2,w) * dLineNodalPBasis(3-q,u)
               EdgeBasis(2*(k-1)+2,2) = 1.5d0 * LineNodalPBasis(2,w) * v * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+2,1) = -1.5d0 * 0.5d0 * v * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+2,3) = 1.5d0 * LineNodalPBasis(2,w) * v * dLineNodalPBasis(3-q,u)
               IF (nj<ni) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO          

             ! Edges 15 and 48 ...
             DO q=1,2
               k = 8+3*(q-1)+1 ! Edge number k: 9 ~ 15 and 12 ~ 48 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               ni = Element % NodeIndexes(i)
               IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
               nj = Element % NodeIndexes(j)
               IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
 
               EdgeBasis(2*(k-1)+1,3) = 0.5d0 * LineNodalPBasis(1,u) * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,1) = 0.5d0 * LineNodalPBasis(1,u) * dLineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,2) = -0.5d0 * dLineNodalPBasis(1,u) * LineNodalPBasis(q,v)
               EdgeBasis(2*(k-1)+2,3) = 1.5d0 * LineNodalPBasis(1,u) * w * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,1) = 1.5d0 * LineNodalPBasis(1,u) * w * dLineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,2) = -1.5d0 * dLineNodalPBasis(1,u) * w * LineNodalPBasis(q,v)
               IF (nj<ni) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO         

             ! Edges 26 and 37 ...
             DO q=1,2
               k = 9+q ! Edge number k: 10 ~ 26 and 11 ~ 37 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               ni = Element % NodeIndexes(i)
               IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
               nj = Element % NodeIndexes(j)
               IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
 
               EdgeBasis(2*(k-1)+1,3) = 0.5d0 * LineNodalPBasis(2,u) * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,1) = 0.5d0 * LineNodalPBasis(2,u) * dLineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,2) = -0.5d0 * dLineNodalPBasis(2,u) * LineNodalPBasis(q,v)
               EdgeBasis(2*(k-1)+2,3) = 1.5d0 * LineNodalPBasis(2,u) * w * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,1) = 1.5d0 * LineNodalPBasis(2,u) * w * dLineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,2) = -1.5d0 * dLineNodalPBasis(2,u) * w * LineNodalPBasis(q,v)
               IF (nj<ni) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO     

             ! ---------------------------------------------------------------------
             ! Additional basis functions on the square faces (four per face).
             ! ---------------------------------------------------------------------         

             ! Faces 1234 and 5678:
             DO q=1,2
               SELECT CASE(q)
               CASE(1)
                 SquareFaceMap(:) = (/ 1,2,3,4 /)
               CASE(2)
                 SquareFaceMap(:) = (/ 5,6,7,8 /)
               END SELECT

               WorkBasis = 0.0d0
               WorkCurlBasis = 0.0d0

               WorkBasis(1,1) = 2.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * LineNodalPBasis(q,w)
               WorkCurlBasis(1,2) = 2.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * dLineNodalPBasis(q,w)
               WorkCurlBasis(1,3) = v * LineNodalPBasis(q,w)

               WorkBasis(2,1) = 12.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * u * LineNodalPBasis(q,w)
               WorkCurlBasis(2,2) = 12.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * u * dLineNodalPBasis(q,w)
               WorkCurlBasis(2,3) = -12.0d0 * (-0.5d0 * v) * u * dLineNodalPBasis(q,w) 

               WorkBasis(3,2) = 2.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * LineNodalPBasis(q,w)
               WorkCurlBasis(3,1) = -2.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * dLineNodalPBasis(q,w)
               WorkCurlBasis(3,3) = -u * LineNodalPBasis(q,w)
               
               WorkBasis(4,2) = 12.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * v * LineNodalPBasis(q,w)
               WorkCurlBasis(4,1) = -12.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * v * dLineNodalPBasis(q,w)
               WorkCurlBasis(4,3) = 12.0d0 * (-0.5d0 * u) * v * LineNodalPBasis(q,w)
               
               DO j=1,4
                 FaceIndeces(j) = Ind(SquareFaceMap(j))
               END DO
               IF (Parallel) THEN
                 DO j=1,4
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                 END DO
               END IF
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

               k = 24
               EdgeBasis(k+4*(q-1)+1,:) = D1 * WorkBasis(2*(I1-1)+1,:)
               CurlBasis(k+4*(q-1)+1,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
               EdgeBasis(k+4*(q-1)+2,:) = WorkBasis(2*(I1-1)+2,:)
               CurlBasis(k+4*(q-1)+2,:) = WorkCurlBasis(2*(I1-1)+2,:)
               EdgeBasis(k+4*(q-1)+3,:) = D2 * WorkBasis(2*(I2-1)+1,:)
               CurlBasis(k+4*(q-1)+3,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
               EdgeBasis(k+4*(q-1)+4,:) = WorkBasis(2*(I2-1)+2,:)
               CurlBasis(k+4*(q-1)+4,:) = WorkCurlBasis(2*(I2-1)+2,:)
             END DO

             ! Faces 1265 and 4378:
             DO q=1,2
               SELECT CASE(q)
               CASE(1)
                 SquareFaceMap(:) = (/ 1,2,6,5 /)
                 k = 32
               CASE(2)
                 SquareFaceMap(:) = (/ 4,3,7,8 /)
                 k = 40
               END SELECT

               WorkBasis = 0.0d0
               WorkCurlBasis = 0.0d0

               WorkBasis(1,1) = 2.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * LineNodalPBasis(q,v)
               WorkCurlBasis(1,2) = 2.0d0 * (-0.5d0 * w) * LineNodalPBasis(q,v)
               WorkCurlBasis(1,3) = -2.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * dLineNodalPBasis(q,v)

               WorkBasis(2,1) = 12.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * u * LineNodalPBasis(q,v)
               WorkCurlBasis(2,2) = 12.0d0 * (-0.5d0 * w) * u * LineNodalPBasis(q,v)
               WorkCurlBasis(2,3) = -12.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * u * dLineNodalPBasis(q,v)

               WorkBasis(3,3) = 2.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * LineNodalPBasis(q,v)
               WorkCurlBasis(3,1) = 2.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * dLineNodalPBasis(q,v)
               WorkCurlBasis(3,2) = u * LineNodalPBasis(q,v)
               
               WorkBasis(4,3) = 12.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * w * LineNodalPBasis(q,v)
               WorkCurlBasis(4,1) = 12.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * w * dLineNodalPBasis(q,v)
               WorkCurlBasis(4,2) = -12.0d0 * (-0.5d0 * u) * w * LineNodalPBasis(q,v)
               
               DO j=1,4
                 FaceIndeces(j) = Ind(SquareFaceMap(j))
               END DO
               IF (Parallel) THEN
                 DO j=1,4
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                 END DO
               END IF
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

               EdgeBasis(k+1,:) = D1 * WorkBasis(2*(I1-1)+1,:)
               CurlBasis(k+1,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
               EdgeBasis(k+2,:) = WorkBasis(2*(I1-1)+2,:)
               CurlBasis(k+2,:) = WorkCurlBasis(2*(I1-1)+2,:)
               EdgeBasis(k+3,:) = D2 * WorkBasis(2*(I2-1)+1,:)
               CurlBasis(k+3,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
               EdgeBasis(k+4,:) = WorkBasis(2*(I2-1)+2,:)
               CurlBasis(k+4,:) = WorkCurlBasis(2*(I2-1)+2,:)
             END DO
             
             ! Faces 2376 and 1485:
             DO q=1,2
               SELECT CASE(q)
               CASE(1)
                 SquareFaceMap(:) = (/ 1,4,8,5 /)
                 k = 44
               CASE(2)
                 SquareFaceMap(:) = (/ 2,3,7,6 /)
                 k = 36
               END SELECT

               WorkBasis = 0.0d0
               WorkCurlBasis = 0.0d0

               WorkBasis(1,2) = 2.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * LineNodalPBasis(q,u)
               WorkCurlBasis(1,1) = -2.0d0 * (-0.5d0 * w) * LineNodalPBasis(q,u)
               WorkCurlBasis(1,3) = 2.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * dLineNodalPBasis(q,u)

               WorkBasis(2,2) = 12.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * v * LineNodalPBasis(q,u)
               WorkCurlBasis(2,1) = -12.0d0 * (-0.5d0 * w) * v * LineNodalPBasis(q,u)
               WorkCurlBasis(2,3) = 12.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * v * dLineNodalPBasis(q,u)

               WorkBasis(3,3) = 2.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * LineNodalPBasis(q,u)
               WorkCurlBasis(3,1) = 2.0d0 * (-0.5d0 * v) * LineNodalPBasis(q,u)
               WorkCurlBasis(3,2) = -2.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * dLineNodalPBasis(q,u)
               
               WorkBasis(4,3) = 12.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * w * LineNodalPBasis(q,u)
               WorkCurlBasis(4,1) = 12.0d0 * (-0.5d0 * v) * w * LineNodalPBasis(q,u)
               WorkCurlBasis(4,2) = -12.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * w * dLineNodalPBasis(q,u)
               
               DO j=1,4
                 FaceIndeces(j) = Ind(SquareFaceMap(j))
               END DO
               IF (Parallel) THEN
                 DO j=1,4
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                 END DO
               END IF
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

               EdgeBasis(k+1,:) = D1 * WorkBasis(2*(I1-1)+1,:)
               CurlBasis(k+1,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
               EdgeBasis(k+2,:) = WorkBasis(2*(I1-1)+2,:)
               CurlBasis(k+2,:) = WorkCurlBasis(2*(I1-1)+2,:)
               EdgeBasis(k+3,:) = D2 * WorkBasis(2*(I2-1)+1,:)
               CurlBasis(k+3,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
               EdgeBasis(k+4,:) = WorkBasis(2*(I2-1)+2,:)
               CurlBasis(k+4,:) = WorkCurlBasis(2*(I2-1)+2,:)
             END DO

             ! Interior basis functions, two per coordinate direction:

             EdgeBasis(49,1) = 8.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * &
                 LineNodalPBasis(1,v) * LineNodalPBasis(2,v)
             CurlBasis(49,2) = 8.0d0 * (-0.5d0 * w) * LineNodalPBasis(1,v) * LineNodalPBasis(2,v)
             CurlBasis(49,3) = -8.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * (-0.5d0 * v)

             EdgeBasis(50,1) = 24.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * u * &
                 LineNodalPBasis(1,v) * LineNodalPBasis(2,v)
             CurlBasis(50,2) = 24.0d0 * (-0.5d0 * w) * u * LineNodalPBasis(1,v) * LineNodalPBasis(2,v)
             CurlBasis(50,3) = -24.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * u *  (-0.5d0 * v)

 
             EdgeBasis(51,2) = 8.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * &
                 LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(51,1) = -8.0d0 * (-0.5d0 * w) * LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(51,3) = 8.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * (-0.5d0 * u)

             EdgeBasis(52,2) = 24.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * v * &
                 LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(52,1) = -24.0d0 * (-0.5d0 * w) * v * LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(52,3) = 24.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * v * (-0.5d0 * u)
            
             EdgeBasis(53,3) = 8.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * &
                 LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(53,1) = 8.0d0 * (-0.5d0 * v) * LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(53,2) = -8.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * (-0.5d0 * u)

             EdgeBasis(54,3) = 24.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * w * &
                 LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(54,1) = 24.0d0 * (-0.5d0 * v) * w * LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(54,2) = -24.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * w * (-0.5d0 * u)

           ELSE
             !--------------------------------------------------------------
             ! The lowest-order element from the optimal family. The optimal
             ! accuracy is obtained also for non-affine meshes.
             ! -------------------------------------------------------------
             ! First twelwe basis functions associated with the edges
             ! -------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(1,1) = ((-1.0d0 + v)*v*(-1.0d0 + w)*w)/8.0d0
             EdgeBasis(1,2) = 0.0d0
             EdgeBasis(1,3) = 0.0d0
             CurlBasis(1,1) = 0.0d0
             CurlBasis(1,2) = ((-1.0d0 + v)*v*(-1.0d0 + 2*w))/8.0d0
             CurlBasis(1,3) = -((-1.0d0 + 2*v)*(-1.0d0 + w)*w)/8.0d0
             IF (nj<ni) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,:) = -CurlBasis(1,:)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(2,1) = 0.0d0
             EdgeBasis(2,2) = (u*(1.0d0 + u)*(-1.0d0 + w)*w)/8.0d0
             EdgeBasis(2,3) = 0.0d0
             CurlBasis(2,1) = -(u*(1.0d0 + u)*(-1.0d0 + 2*w))/8.0d0
             CurlBasis(2,2) = 0.0d0
             CurlBasis(2,3) = ((1.0d0 + 2*u)*(-1.0d0 + w)*w)/8.0d0
             IF (nj<ni) THEN
               EdgeBasis(2,:) = -EdgeBasis(2,:)
               CurlBasis(2,:) = -CurlBasis(2,:)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(3,1) = (v*(1.0d0 + v)*(-1.0d0 + w)*w)/8.0d0
             EdgeBasis(3,2) = 0.0d0
             EdgeBasis(3,3) = 0.0d0
             CurlBasis(3,1) = 0.0d0
             CurlBasis(3,2) = (v*(1.0d0 + v)*(-1.0d0 + 2*w))/8.0d0
             CurlBasis(3,3) = -((1.0d0 + 2*v)*(-1.0d0 + w)*w)/8.0d0
             IF (nj<ni) THEN
               EdgeBasis(3,:) = -EdgeBasis(3,:)
               CurlBasis(3,:) = -CurlBasis(3,:)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(4,1) = 0.0d0
             EdgeBasis(4,2) = ((-1.0d0 + u)*u*(-1.0d0 + w)*w)/8.0d0
             EdgeBasis(4,3) = 0.0d0
             CurlBasis(4,1) = -((-1.0d0 + u)*u*(-1.0d0 + 2*w))/8.0d0
             CurlBasis(4,2) = 0.0d0
             CurlBasis(4,3) = ((-1.0d0 + 2*u)*(-1.0d0 + w)*w)/8.0d0
             IF (nj<ni) THEN
               EdgeBasis(4,:) = -EdgeBasis(4,:)
               CurlBasis(4,:) = -CurlBasis(4,:)
             END IF

             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(5,1) = ((-1.0d0 + v)*v*w*(1.0d0 + w))/8.0d0
             EdgeBasis(5,2) = 0.0d0
             EdgeBasis(5,3) = 0.0d0
             CurlBasis(5,1) = 0.0d0
             CurlBasis(5,2) = ((-1.0d0 + v)*v*(1.0d0 + 2*w))/8.0d0 
             CurlBasis(5,3) = -((-1.0d0 + 2*v)*w*(1.0d0 + w))/8.0d0
             IF (nj<ni) THEN
               EdgeBasis(5,:) = -EdgeBasis(5,:)
               CurlBasis(5,:) = -CurlBasis(5,:)
             END IF

             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(6,1) = 0.0d0
             EdgeBasis(6,2) = (u*(1.0d0 + u)*w*(1.0d0 + w))/8.0d0
             EdgeBasis(6,3) = 0.0d0
             CurlBasis(6,1) = -(u*(1.0d0 + u)*(1.0d0 + 2*w))/8.0d0
             CurlBasis(6,2) = 0.0d0
             CurlBasis(6,3) = ((1.0d0 + 2*u)*w*(1.0d0 + w))/8.0d0
             IF (nj<ni) THEN
               EdgeBasis(6,:) = -EdgeBasis(6,:)
               CurlBasis(6,:) = -CurlBasis(6,:)
             END IF

             i = EdgeMap(7,1)
             j = EdgeMap(7,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(7,1) = (v*(1.0d0 + v)*w*(1.0d0 + w))/8.0d0
             EdgeBasis(7,2) = 0.0d0
             EdgeBasis(7,3) = 0.0d0
             CurlBasis(7,1) = 0.0d0
             CurlBasis(7,2) = (v*(1.0d0 + v)*(1.0d0 + 2*w))/8.0d0
             CurlBasis(7,3) = -((1.0d0 + 2*v)*w*(1.0d0 + w))/8.0d0
             IF (nj<ni) THEN
               EdgeBasis(7,:) = -EdgeBasis(7,:)
               CurlBasis(7,:) = -CurlBasis(7,:)
             END IF

             i = EdgeMap(8,1)
             j = EdgeMap(8,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(8,1) = 0.0d0
             EdgeBasis(8,2) = ((-1.0d0 + u)*u*w*(1.0d0 + w))/8.0d0
             EdgeBasis(8,3) = 0.0d0
             CurlBasis(8,1) = -((-1.0d0 + u)*u*(1.0d0 + 2*w))/8.0d0
             CurlBasis(8,2) = 0.0d0
             CurlBasis(8,3) = ((-1.0d0 + 2*u)*w*(1.0d0 + w))/8.0d0
             IF (nj<ni) THEN
               EdgeBasis(8,:) = -EdgeBasis(8,:)
               CurlBasis(8,:) = -CurlBasis(8,:)
             END IF

             i = EdgeMap(9,1)
             j = EdgeMap(9,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(9,1) = 0.0d0
             EdgeBasis(9,2) = 0.0d0
             EdgeBasis(9,3) = ((-1.0d0 + u)*u*(-1.0d0 + v)*v)/8.0d0
             CurlBasis(9,1) = ((-1.0d0 + u)*u*(-1.0d0 + 2*v))/8.0d0
             CurlBasis(9,2) = -((-1.0d0 + 2*u)*(-1.0d0 + v)*v)/8.0d0
             CurlBasis(9,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(9,:) = -EdgeBasis(9,:)
               CurlBasis(9,:) = -CurlBasis(9,:)
             END IF

             i = EdgeMap(10,1)
             j = EdgeMap(10,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(10,1) = 0.0d0
             EdgeBasis(10,2) = 0.0d0
             EdgeBasis(10,3) = (u*(1.0d0 + u)*(-1.0d0 + v)*v)/8.0d0
             CurlBasis(10,1) = (u*(1.0d0 + u)*(-1.0d0 + 2*v))/8.0d0
             CurlBasis(10,2) = -((1.0d0 + 2*u)*(-1.0d0 + v)*v)/8.0d0
             CurlBasis(10,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(10,:) = -EdgeBasis(10,:)
               CurlBasis(10,:) = -CurlBasis(10,:)
             END IF

             i = EdgeMap(11,1)
             j = EdgeMap(11,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(11,1) = 0.0d0
             EdgeBasis(11,2) = 0.0d0
             EdgeBasis(11,3) = (u*(1.0d0 + u)*v*(1.0d0 + v))/8.0d0
             CurlBasis(11,1) = (u*(1.0d0 + u)*(1.0d0 + 2*v))/8.0d0
             CurlBasis(11,2) = -((1.0d0 + 2*u)*v*(1.0d0 + v))/8.0d0
             CurlBasis(11,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(11,:) = -EdgeBasis(11,:)
               CurlBasis(11,:) = -CurlBasis(11,:)
             END IF

             i = EdgeMap(12,1)
             j = EdgeMap(12,2)
             ni = Element % NodeIndexes(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Element % NodeIndexes(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             EdgeBasis(12,1) = 0.0d0
             EdgeBasis(12,2) = 0.0d0
             EdgeBasis(12,3) = ((-1.0d0 + u)*u*v*(1.0d0 + v))/8.0d0
             CurlBasis(12,1) = ((-1.0d0 + u)*u*(1.0d0 + 2*v))/8.0d0
             CurlBasis(12,2) = -((-1.0d0 + 2*u)*v*(1.0d0 + v))/8.0d0
             CurlBasis(12,3) = 0.0d0
             IF (nj<ni) THEN
               EdgeBasis(12,:) = -EdgeBasis(12,:)
               CurlBasis(12,:) = -CurlBasis(12,:)
             END IF

             ! ---------------------------------------------------------------------
             ! Additional twelwe basis functions on the square faces (two per face).
             ! ---------------------------------------------------------------------         
             BrickFaceMap(1,:) = (/ 1,2,3,4 /)          
             BrickFaceMap(2,:) = (/ 5,6,7,8 /)
             BrickFaceMap(3,:) = (/ 1,2,6,5 /)
             BrickFaceMap(4,:) = (/ 2,3,7,6 /)
             BrickFaceMap(5,:) = (/ 4,3,7,8 /)
             BrickFaceMap(6,:) = (/ 1,4,8,5 /)

             ! The first face:
             WorkBasis(1,1) = -((-1.0d0 + v**2)*(-1.0d0 + w)*w)/4.0d0
             WorkBasis(1,2) = 0.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = 0.0d0
             WorkCurlBasis(1,2) = -((-1.0d0 + v**2)*(-1.0d0 + 2*w))/4.0d0
             WorkCurlBasis(1,3) = (v*(-1.0d0 + w)*w)/2.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = -((-1.0d0 + u**2)*(-1.0d0 + w)*w)/4.0d0
             WorkBasis(2,3) = 0.0d0
             WorkCurlBasis(2,1) = ((-1.0d0 + u**2)*(-1.0d0 + 2*w))/4.0d0
             WorkCurlBasis(2,2) = 0.0d0
             WorkCurlBasis(2,3) = -(u*(-1.0d0 + w)*w)/2.0d0

             DO j=1,4
               FaceIndeces(j) = Ind(BrickFaceMap(1,j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(13,:) = D1 * WorkBasis(I1,:)
             CurlBasis(13,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(14,:) = D2 * WorkBasis(I2,:)
             CurlBasis(14,:) = D2 * WorkCurlBasis(I2,:) 

             ! The second face:
             WorkBasis(1,1) = -((-1.0d0 + v**2)*w*(1.0d0 + w))/4.0d0
             WorkBasis(1,2) = 0.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = 0.0d0
             WorkCurlBasis(1,2) = -((-1.0d0 + v**2)*(1.0d0 + 2*w))/4.0d0
             WorkCurlBasis(1,3) = (v*w*(1.0d0 + w))/2.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = -((-1.0d0 + u**2)*w*(1.0d0 + w))/4.0d0
             WorkBasis(2,3) = 0.0d0
             WorkCurlBasis(2,1) = ((-1.0d0 + u**2)*(1.0d0 + 2*w))/4.0d0
             WorkCurlBasis(2,2) = 0.0d0
             WorkCurlBasis(2,3) = -(u*w*(1.0d0 + w))/2.0d0

             DO j=1,4
               FaceIndeces(j) = Ind(BrickFaceMap(2,j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(15,:) = D1 * WorkBasis(I1,:)
             CurlBasis(15,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(16,:) = D2 * WorkBasis(I2,:)
             CurlBasis(16,:) = D2 * WorkCurlBasis(I2,:) 

             ! The third face:
             WorkBasis(1,1) = -((-1.0d0 + v)*v*(-1.0d0 + w**2))/4.0d0
             WorkBasis(1,2) = 0.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = 0.0d0
             WorkCurlBasis(1,2) = -((-1.0d0 + v)*v*w)/2.0d0
             WorkCurlBasis(1,3) = ((-1.0d0 + 2*v)*(-1.0d0 + w**2))/4.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = -((-1.0d0 + u**2)*(-1.0d0 + v)*v)/4.0d0
             WorkCurlBasis(2,1) = -((-1.0d0 + u**2)*(-1.0d0 + 2*v))/4.0d0
             WorkCurlBasis(2,2) = (u*(-1.0d0 + v)*v)/2.0d0
             WorkCurlBasis(2,3) = 0.0d0

             DO j=1,4
               FaceIndeces(j) = Ind(BrickFaceMap(3,j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(17,:) = D1 * WorkBasis(I1,:)
             CurlBasis(17,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(18,:) = D2 * WorkBasis(I2,:)
             CurlBasis(18,:) = D2 * WorkCurlBasis(I2,:) 

             ! The fourth face:
             WorkBasis(1,1) = 0.0d0
             WorkBasis(1,2) = -(u*(1.0d0 + u)*(-1.0d0 + w**2))/4.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = (u*(1.0d0 + u)*w)/2.0d0
             WorkCurlBasis(1,2) = 0.0d0
             WorkCurlBasis(1,3) = -((1.0d0 + 2*u)*(-1.0d0 + w**2))/4.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = -(u*(1.0d0 + u)*(-1 + v**2))/4.0d0
             WorkCurlBasis(2,1) = -(u*(1.0d0 + u)*v)/2.0d0
             WorkCurlBasis(2,2) = ((1.0d0 + 2*u)*(-1.0d0 + v**2))/4.0d0
             WorkCurlBasis(2,3) = 0.0d0

             DO j=1,4
               FaceIndeces(j) = Ind(BrickFaceMap(4,j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(19,:) = D1 * WorkBasis(I1,:)
             CurlBasis(19,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(20,:) = D2 * WorkBasis(I2,:)
             CurlBasis(20,:) = D2 * WorkCurlBasis(I2,:) 

             ! The fifth face:
             WorkBasis(1,1) = -(v*(1.0d0 + v)*(-1.0d0 + w**2))/4.0d0
             WorkBasis(1,2) = 0.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = 0.0d0
             WorkCurlBasis(1,2) = -(v*(1.0d0 + v)*w)/2.0d0
             WorkCurlBasis(1,3) = ((1.0d0 + 2*v)*(-1.0d0 + w**2))/4.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = -((-1.0d0 + u**2)*v*(1.0d0 + v))/4.0d0
             WorkCurlBasis(2,1) = -((-1.0d0 + u**2)*(1.0d0 + 2*v))/4.0d0
             WorkCurlBasis(2,2) = (u*v*(1.0d0 + v))/2.0d0
             WorkCurlBasis(2,3) = 0.0d0

             DO j=1,4
               FaceIndeces(j) = Ind(BrickFaceMap(5,j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(21,:) = D1 * WorkBasis(I1,:)
             CurlBasis(21,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(22,:) = D2 * WorkBasis(I2,:)
             CurlBasis(22,:) = D2 * WorkCurlBasis(I2,:) 

             ! The sixth face:
             WorkBasis(1,1) = 0.0d0
             WorkBasis(1,2) = -((-1.0d0 + u)*u*(-1.0d0 + w**2))/4.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = ((-1.0d0 + u)*u*w)/2.0d0
             WorkCurlBasis(1,2) = 0.0d0
             WorkCurlBasis(1,3) = -((-1.0d0 + 2*u)*(-1.0d0 + w**2))/4.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = -((-1.0d0 + u)*u*(-1.0d0 + v**2))/4.0d0
             WorkCurlBasis(2,1) = -((-1.0d0 + u)*u*v)/2.0d0
             WorkCurlBasis(2,2) = ((-1.0d0 + 2*u)*(-1.0d0 + v**2))/4.0d0
             WorkCurlBasis(2,3) = 0.0d0

             DO j=1,4
               FaceIndeces(j) = Ind(BrickFaceMap(6,j))
             END DO
             IF (Parallel) THEN
               DO j=1,4
                 FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
               END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)

             EdgeBasis(23,:) = D1 * WorkBasis(I1,:)
             CurlBasis(23,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(24,:) = D2 * WorkBasis(I2,:)
             CurlBasis(24,:) = D2 * WorkCurlBasis(I2,:) 

             ! ------------------------------------------------------------------------
             ! Additional basis functions on the element interior (three per element)
             ! -----------------------------------------------------------------------
             EdgeBasis(25,1) = ((-1.0d0 + v**2)*(-1.0d0 + w**2))/2.0d0
             EdgeBasis(25,2) = 0.0d0
             EdgeBasis(25,3) = 0.0d0
             CurlBasis(25,1) = 0.0d0
             CurlBasis(25,2) = (-1.0d0 + v**2)*w
             CurlBasis(25,3) = v - v*w**2

             EdgeBasis(26,1) = 0.0d0
             EdgeBasis(26,2) = ((-1.0d0 + u**2)*(-1.0d0 + w**2))/2.0d0
             EdgeBasis(26,3) = 0.0d0
             CurlBasis(26,1) = w - u**2*w
             CurlBasis(26,2) = 0.0d0
             CurlBasis(26,3) = u*(-1 + w**2)

             EdgeBasis(27,1) = 0.0d0
             EdgeBasis(27,2) = 0.0d0
             EdgeBasis(27,3) = ((-1.0d0 + u**2)*(-1.0d0 + v**2))/2.0d0
             CurlBasis(27,1) = (-1.0d0 + u**2)*v
             CurlBasis(27,2) = u - u*v**2
             CurlBasis(27,3) = 0.0d0
           END IF

         CASE DEFAULT
           CALL Fatal('ElementDescription::EdgeElementInfo','Unsupported element type')
         END SELECT
       END IF

       IF (cdim == dim) THEN
          !--------------------------------------------------------------------------------
          ! To optimize computation, this branch avoids calling the ElementMetric function
          ! since all necessary data has already been found via PiolaTransformationData.
          !-------------------------------------------------------------------------------
          IF (PerformPiolaTransform) THEN
             DO j=1,DOFs
                DO k=1,dim
                   B(k) = SUM( LG(k,1:dim) * EdgeBasis(j,1:dim) )
                END DO
                EdgeBasis(j,1:dim) = B(1:dim)

                IF (dim == 2) THEN
                   CurlBasis(j,3) = 1.0d0/DetF * CurlBasis(j,3)
                ELSE
                   DO k=1,dim
                      B(k) = 1.0d0/DetF * SUM( LF(k,1:dim) * CurlBasis(j,1:dim) )
                   END DO
                   CurlBasis(j,1:dim) = B(1:dim)
                END IF
             END DO
             ! Make the returned value DetF to act as a metric term for integration
             ! over the volume of the element: 
             DetF = ABS(DetF)
          END IF

          ! ----------------------------------------------------------------------
          ! Get global first derivatives of the nodal basis functions if wanted:
          ! ----------------------------------------------------------------------
          IF ( PRESENT(dBasisdx) ) THEN
             dBasisdx = 0.0d0
             DO i=1,n
                DO j=1,dim
                   DO k=1,dim
                      dBasisdx(i,j) = dBasisdx(i,j) + dLBasisdx(i,k)*LG(j,k)
                   END DO
                END DO
             END DO
          END IF
       ELSE

          IF (PerformPiolaTransform .OR. PRESENT(dBasisdx) .OR. ApplyTraceMapping) THEN
             IF ( .NOT. ElementMetric( n, Element, Nodes, &
                  ElmMetric, detJ, dLBasisdx, LG ) ) THEN
                stat = .FALSE.
                RETURN
             END IF
          END IF

          IF (ApplyTraceMapping .AND. (dim==2) ) THEN
            ! Perform operation b -> b x n. The resulting field transforms under the usual 
            ! Piola transform (like div-conforming field). For a general surface element
            ! embedded in 3D we return B(f(p))=1/sqrt(a) F(b x n) where a is the determinant of
            ! the metric tensor, F=[a1 a2] with a1 and a2 surface basis vectors and (b x n) is
            ! considered to be 2-vector (the trivial component ignored). Note that asking simultaneously 
            ! for the curl of the basis is not an expected combination.
            DO j=1,DOFs
              WorkBasis(1,1:2) = EdgeBasis(j,1:2)
              EdgeBasis(j,1) = WorkBasis(1,2)
              EdgeBasis(j,2) = -WorkBasis(1,1)
            END DO
            IF (PerformPiolaTransform) THEN
              DO j=1,DOFs 
                DO k=1,cdim
                  B(k) = SUM( LF(k,1:dim) * EdgeBasis(j,1:dim) ) / DetJ
                END DO
                EdgeBasis(j,1:cdim) = B(1:cdim)                
              END DO
            END IF
          ELSE
            IF (PerformPiolaTransform) THEN
              DO j=1,DOFs
                DO k=1,cdim
                  B(k) = SUM( LG(k,1:dim) * EdgeBasis(j,1:dim) )
                END DO
                EdgeBasis(j,1:cdim) = B(1:cdim)
                ! The returned spatial curl in the case cdim=3 and dim=2 handled here
                ! has limited usability. This handles only a transformation of
                ! the type x_3 = p_3:
                CurlBasis(j,3) = 1.0d0/DetJ * CurlBasis(j,3)
              END DO
            END IF
          END IF

          ! Make the returned value DetF to act as a metric term for integration
          ! over the volume of the element: 
          DetF = DetJ

          ! ----------------------------------------------------------------------
          ! Get global first derivatives of the nodal basis functions if wanted:
          ! ----------------------------------------------------------------------
          IF ( PRESENT(dBasisdx) ) THEN
             dBasisdx = 0.0d0
             DO i=1,n
                DO j=1,cdim
                   DO k=1,dim
                      dBasisdx(i,j) = dBasisdx(i,j) + dLBasisdx(i,k)*LG(j,k)
                   END DO
                END DO
             END DO
          END IF

       END IF

       IF(PRESENT(F)) F = LF
       IF(PRESENT(G)) G = LG
       IF(PRESENT(RotBasis)) RotBasis(1:DOFs,:) = CurlBasis(1:DOFs,:)
!-----------------------------------------------------------------------------
     END FUNCTION EdgeElementInfo
!------------------------------------------------------------------------------



!----------------------------------------------------------------------------
     SUBROUTINE TriangleFaceDofsOrdering(I1,I2,D1,D2,Ind)       
!-----------------------------------------------------------------------------
! This is used for selecting what additional basis functions are associated
! with a triangular face in the case of second-order approximation.
! ----------------------------------------------------------------------------
       INTEGER ::  I1, I2, Ind(4)
       REAL(KIND=dp) :: D1, D2
!---------------------------------------------------------------------------
       INTEGER ::  k, A
! --------------------------------------------------------------------------
       D1 = 1.0d0
       D2 = 1.0d0
       IF ( Ind(1) < Ind(2) ) THEN
          k = 1
       ELSE
          k = 2
       END IF
       IF ( Ind(k) > Ind(3) ) THEN
          k = 3
       END IF
       A = k

       SELECT CASE(A)
       CASE(1)
          IF (Ind(3) > Ind(2)) THEN
             ! C = 3
             I1 = 1
             I2 = 2
          ELSE
             ! C = 2
             I1 = 2
             I2 = 1             
          END IF
       CASE(2)
         IF (Ind(3) > Ind(1)) THEN
             ! C = 3
             I1 = 1
             I2 = 3
             D1 = -1.0d0
          ELSE
             ! C = 1
             I1 = 3
             I2 = 1
             D2 = -1.0d0             
          END IF
       CASE(3)
          IF (Ind(2) > Ind(1)) THEN
             ! C = 2
             I1 = 2
             I2 = 3
          ELSE
             ! C = 1
             I1 = 3
             I2 = 2
          END IF
          D1 = -1.0d0
          D2 = -1.0d0          
       CASE DEFAULT
          CALL Fatal('ElementDescription::TriangleFaceDofsOrdering','Erratic square face indeces')
       END SELECT
!---------------------------------------------------------
     END SUBROUTINE TriangleFaceDofsOrdering
!-----------------------------------------------------------


!-------------------------------------------------------------
     SUBROUTINE TriangleFaceDofsOrdering2(t,s,Ind)       
!-------------------------------------------------------------------------------
! Returns two unit vectors t and s for spanning constant vector fields
! defined on a triangular face. As a rule for orientation, the vector t is defined 
! as t = Grad L_B - Grad L_A where L_A and L_B are the Lagrange basis functions
! associated with the nodes that has the smallest global indices A and B (A<B).
! Then s = Sqrt(3)* grad L_C, with C corresponding to the largest global index.
!-------------------------------------------------------------------------------
       INTEGER ::  Ind(4)
       REAL(KIND=dp) :: t(3), s(3)
!----------------------------------------------------------
       INTEGER ::  k, A
! -------------------------------------------------------------------
       t = 0.0d0
       s = 0.0d0

       IF ( Ind(1) < Ind(2) ) THEN
          k = 1
       ELSE
          k = 2
       END IF
       IF ( Ind(k) > Ind(3) ) THEN
          k = 3
       END IF
       A = k

       SELECT CASE(A)
       CASE(1)
          IF ( Ind(2) < Ind(3) ) THEN ! B=2, tangent = AB = 12
             t(1) = 1.0d0
             t(2) = 0.0
             s(1) = 0.0d0
             s(2) = 1.0d0
          ELSE ! B=3, tangent = AB = 13
             t(1) = 0.5d0
             t(2) = Sqrt(3.0d0)/2.0d0
             s(1) = Sqrt(3.0d0)/2.0d0
             s(2) = -0.5d0
          END IF
       CASE(2)     
          IF ( Ind(1) < Ind(3) ) THEN ! B=1, tangent = AB = 21
             t(1) = -1.0d0
             t(2) = 0.0
             s(1) = 0.0d0
             s(2) = 1.0d0
          ELSE ! B=3, tangent = AB = 23
             t(1) = -0.5d0
             t(2) = Sqrt(3.0d0)/2.0d0
             s(1) = -Sqrt(3.0d0)/2.0d0
             s(2) = -0.5d0
          END IF
       CASE(3)
          IF ( Ind(1) < Ind(2) ) THEN ! B=1, tangent = AB = 31
             t(1) = -0.5d0
             t(2) = -Sqrt(3.0d0)/2.0d0
             s(1) = Sqrt(3.0d0)/2.0d0
             s(2) = -0.5d0          
          ELSE ! B=2, tangent = AB = 32
             t(1) = 0.5d0
             t(2) = -Sqrt(3.0d0)/2.0d0            
             s(1) = -Sqrt(3.0d0)/2.0d0
             s(2) = -0.5d0       
          END IF
       CASE DEFAULT
          CALL Fatal('ElementDescription::TriangleFaceDofsOrdering','Erratic square face indeces')
       END SELECT
!---------------------------------------------------------
     END SUBROUTINE TriangleFaceDofsOrdering2
!-----------------------------------------------------------


!---------------------------------------------------------
     SUBROUTINE SquareFaceDofsOrdering(I1,I2,D1,D2,Ind)       
!-----------------------------------------------------------
       INTEGER ::  I1, I2, Ind(4)
       REAL(KIND=dp) :: D1, D2
!----------------------------------------------------------
       INTEGER ::  i, j, k, l, A
! -------------------------------------------------------------------
!  Find input for applying an order change and sign reversions to two
!  basis functions associated with a square face. To this end, 
!  find nodes A, B, C such that A has the minimal global index,
!  AB and AC are edges, with C having the largest global index. 
!  Then AB gives the positive direction for the first face DOF and
!  AC gives the positive direction for the second face DOF.
!  REMARK: This convention must be followed when creating basis
!  functions for other element types which are intended to be compatible
!  with the element type to which this rule is applied.
! -------------------------------------------------------------------
       i = 1
       j = 2
       IF ( Ind(i) < Ind(j) ) THEN
          k = i
       ELSE
          k = j
       END IF
       i = 4
       j = 3 
       IF ( Ind(i) < Ind(j) ) THEN
          l = i
       ELSE
          l = j
       END IF
       IF ( Ind(k) > Ind(l) ) THEN
          k = l
       END IF
       A = k

       SELECT CASE(A)
       CASE(1)
          IF ( Ind(2) < Ind(4) ) THEN
             I1 = 1
             I2 = 2
             D1 = 1.0d0
             D2 = 1.0d0
          ELSE
             I1 = 2
             I2 = 1
             D1 = 1.0d0
             D2 = 1.0d0 
          END IF
       CASE(2)
          IF ( Ind(3) < Ind(1) ) THEN
             I1 = 2
             I2 = 1
             D1 = 1.0d0
             D2 = -1.0d0
          ELSE
             I1 = 1
             I2 = 2
             D1 = -1.0d0
             D2 = 1.0d0
          END IF
       CASE(3)
          IF ( Ind(4) < Ind(2) ) THEN
             I1 = 1
             I2 = 2
             D1 = -1.0d0
             D2 = -1.0d0
          ELSE
             I1 = 2
             I2 = 1
             D1 = -1.0d0
             D2 = -1.0d0
          END IF
       CASE(4)
          IF ( Ind(1) < Ind(3) ) THEN
             I1 = 2
             I2 = 1
             D1 = -1.0d0
             D2 = 1.0d0
          ELSE
             I1 = 1
             I2 = 2
             D1 = 1.0d0
             D2 = -1.0d0
          END IF
       CASE DEFAULT
          CALL Fatal('ElementDescription::SquareFaceDofsOrdering','Erratic square face indeces')
       END SELECT
!----------------------------------------------------------
     END SUBROUTINE SquareFaceDofsOrdering
!----------------------------------------------------------

!----------------------------------------------------------------------------------
!>  Returns data for rearranging H(curl)-conforming basis functions so that 
!>  compatibility with the convention for defining global DOFs is attained.
!>  If n basis function value have already been tabulated in the default order
!>  as BasisArray(1:n,:), then SignVec(1:n) * BasisArray(PermVec(1:n),:) gives
!>  the basis vector values corresponding to the global DOFs.
!>  TO DO: support for second-order basis functions, triangles and quads missing
!------------------------------------------------------------------------------------
     SUBROUTINE ReorderingAndSignReversionsData(Element,Nodes,PermVec,SignVec)
!-------------------------------------------------------------------------------------
       IMPLICIT NONE

       TYPE(Element_t), TARGET :: Element        !< Element structure
       TYPE(Nodes_t) :: Nodes                    !< Data corresponding to the classic element nodes
       INTEGER :: PermVec(:)                     !< At exit the permution vector for performing reordering
       REAL(KIND=dp) :: SignVec(:)               !< At exit the vector for performing sign changes
!---------------------------------------------------------------------------------------------------
       TYPE(Mesh_t), POINTER :: Mesh       
       INTEGER, POINTER :: EdgeMap(:,:), Ind(:)
       INTEGER :: SquareFaceMap(4), BrickFaceMap(6,4), PrismSquareFaceMap(3,4), DOFs, i, j, k
       INTEGER :: FaceIndeces(4), I1, I2, ni, nj
       REAL(KIND=dp) :: D1, D2
       LOGICAL :: Parallel
!---------------------------------------------------------------------------------------------------
       Mesh => CurrentModel % Solver % Mesh
       !Parallel = ParEnv % PEs>1       
       Parallel = ASSOCIATED(Mesh % ParallelInfo % Interface)

       SignVec = 1.0d0
       Ind => Element % Nodeindexes

       SELECT CASE( Element % TYPE % ElementCode / 100 )
       !CASE(3) needs to be done

       !CASE(4) needs to be done

       CASE(5)
          ! NOTE: The Nedelec second family is not yet supported
          EdgeMap => LGetEdgeMap(5)
          DO k=1,6
             i = EdgeMap(k,1)
             j = EdgeMap(k,2)
             ni = Ind(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Ind(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) SignVec(k) = -1.0d0
             PermVec(k) = k
          END DO

       CASE(6)
          EdgeMap => LGetEdgeMap(6)
          DO k=1,8
             i = EdgeMap(k,1)
             j = EdgeMap(k,2)
             ni = Ind(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Ind(j) 
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) SignVec(k) = -1.0d0
             PermVec(k) = k
          END DO
          ! -----------------------------------------------------
          ! Additional two basis functions on the square face
          ! -----------------------------------------------------
          SquareFaceMap(:) = (/ 1,2,3,4 /)
          DO j=1,4
             FaceIndeces(j) = Ind(SquareFaceMap(j))
          END DO
          IF (Parallel) THEN
             DO j=1,4
                FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
             END DO
          END IF

          CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)
          i = 8
          PermVec(i+1) = i+I1 
          PermVec(i+2) = i+I2
          SignVec(i+1) = D1
          SignVec(i+2) = D2
 
       CASE(7)
          EdgeMap => LGetEdgeMap(7)
          DO k=1,9
             i = EdgeMap(k,1)
             j = EdgeMap(k,2)
             ni = Ind(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Ind(j)
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) SignVec(k) = -1.0d0
             PermVec(k) = k
          END DO
          ! ---------------------------------------------------------------------
          ! Additional six basis functions on the square faces (two per face).
          ! ---------------------------------------------------------------------         
          PrismSquareFaceMap(1,:) = (/ 1,2,5,4 /)
          PrismSquareFaceMap(2,:) = (/ 2,3,6,5 /)
          PrismSquareFaceMap(3,:) = (/ 3,1,4,6 /)
          DO k=1,3
             DO j=1,4
                FaceIndeces(j) = Ind(PrismSquareFaceMap(k,j))
             END DO
             IF (Parallel) THEN
                DO j=1,4
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)
             i = 9+(k-1)*2
             PermVec(i+1) = i+I1 
             PermVec(i+2) = i+I2
             SignVec(i+1) = D1
             SignVec(i+2) = D2 
          END DO

       CASE(8)
          EdgeMap => LGetEdgeMap(8)
          DO k=1,12
             i = EdgeMap(k,1)
             j = EdgeMap(k,2)
             ni = Ind(i)
             IF (Parallel) ni=Mesh % ParallelInfo % GlobalDOFs(ni)
             nj = Ind(j) 
             IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
             IF (nj<ni) SignVec(k) = -1.0d0
             PermVec(k) = k
          END DO
          ! ---------------------------------------------------------------------
          ! Additional twelwe basis functions on the square faces (two per face).
          ! ---------------------------------------------------------------------         
          BrickFaceMap(1,:) = (/ 1,2,3,4 /)          
          BrickFaceMap(2,:) = (/ 5,6,7,8 /)
          BrickFaceMap(3,:) = (/ 1,2,6,5 /)
          BrickFaceMap(4,:) = (/ 2,3,7,6 /)
          BrickFaceMap(5,:) = (/ 4,3,7,8 /)
          BrickFaceMap(6,:) = (/ 1,4,8,5 /)
          DO k=1,6
             DO j=1,4
                FaceIndeces(j) = Ind(BrickFaceMap(k,j))
             END DO
             IF (Parallel) THEN
                DO j=1,4
                   FaceIndeces(j) = Mesh % ParallelInfo % GlobalDOFs(FaceIndeces(j))
                END DO
             END IF
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndeces)
             i = 12+(k-1)*2
             PermVec(i+1) = i+I1 
             PermVec(i+2) = i+I2
             SignVec(i+1) = D1
             SignVec(i+2) = D2 
          END DO
          PermVec(25) = 25
          PermVec(26) = 26         
          PermVec(27) = 27
           
       CASE DEFAULT
          CALL Fatal('ElementDescription::ReorderingAndSignReversionsData','Unsupported element type')
       END SELECT
!----------------------------------------------------------
     END SUBROUTINE ReorderingAndSignReversionsData
!----------------------------------------------------------


! --------------------------------------------------------------------------------------
!> This subroutine contains an older design for providing edge element basis functions
!> of the lowest-degree. Obtaining optimal accuracy with these elements may require that 
!> the element map is affine, while the edge basis functions given by the newer design 
!> (the function EdgeElementInfo) should also work on general meshes. 
!------------------------------------------------------------------------
   SUBROUTINE GetEdgeBasis( Element, WBasis, RotWBasis, Basis, dBasisdx )
!------------------------------------------------------------------------
     TYPE(Element_t),TARGET :: Element
     REAL(KIND=dp) :: WBasis(:,:), RotWBasis(:,:), Basis(:), dBasisdx(:,:)
!------------------------------------------------------------------------
     TYPE(Element_t),POINTER :: Edge
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Nodes_t), SAVE :: Nodes
     REAL(KIND=dp) :: u,v,w,dudx(3,3),du(3),Base,dBase(3),tBase(3), &
                rBase(3),triBase(3),dtriBase(3,3), G(3,3), F(3,3), detF, detG, &
                EdgeBasis(8,3), CurlBasis(8,3)
     LOGICAL :: Parallel,stat
     INTEGER :: i,j,k,n,nj,nk,i1,i2
     INTEGER, POINTER :: EdgeMap(:,:)
!------------------------------------------------------------------------
     Mesh => CurrentModel % Solver % Mesh
     Parallel = ASSOCIATED(Mesh % ParallelInfo % Interface)

     IF (Element % TYPE % BasisFunctionDegree>1) THEN
       CALL Fatal('GetEdgeBasis',"Can't handle but linear elements, sorry.") 
     END IF

     SELECT CASE(Element % TYPE % ElementCode / 100)
     CASE(4,7,8)
       n = Element % TYPE % NumberOfNodes
       u = SUM(Basis(1:n)*Element % TYPE % NodeU(1:n))
       v = SUM(Basis(1:n)*Element % TYPE % NodeV(1:n))
       w = SUM(Basis(1:n)*Element % TYPE % NodeW(1:n))

       dudx(1,:) = MATMUL(Element % TYPE % NodeU(1:n),dBasisdx(1:n,:))
       dudx(2,:) = MATMUL(Element % TYPE % NodeV(1:n),dBasisdx(1:n,:))
       dudx(3,:) = MATMUL(Element % TYPE % NodeW(1:n),dBasisdx(1:n,:))

       triBase(1) = 1-u-v
       triBase(2) = u
       triBase(3) = v

       dtriBase(1,:) = -dudx(1,:)-dudx(2,:) 
       dtriBase(2,:) =  dudx(1,:)
       dtriBase(3,:) =  dudx(2,:)
     CASE(6)
       n = Element % TYPE % NumberOfNodes
       u = SUM(Basis(1:n)*Element % TYPE % NodeU(1:n))
       v = SUM(Basis(1:n)*Element % TYPE % NodeV(1:n))
       w = SUM(Basis(1:n)*Element % TYPE % NodeW(1:n))

       G(1,:) = MATMUL(Element % TYPE % NodeU(1:n),dBasisdx(1:n,:))
       G(2,:) = MATMUL(Element % TYPE % NodeV(1:n),dBasisdx(1:n,:))
       G(3,:) = MATMUL(Element % TYPE % NodeW(1:n),dBasisdx(1:n,:))            

       detG =  G(1,1) * ( G(2,2)*G(3,3) - G(2,3)*G(3,2) ) + &
                  G(1,2) * ( G(2,3)*G(3,1) - G(2,1)*G(3,3) ) + &
                  G(1,3) * ( G(2,1)*G(3,2) - G(2,2)*G(3,1) )
       detF = 1.0d0/detG
       CALL InvertMatrix3x3(G,F,detG)
       
       !------------------------------------------------------------
       ! The basis functions spanning the reference element space and
       ! their Curl with respect to the local coordinates
       ! ------------------------------------------------------------
       EdgeBasis(1,1) = (1.0d0 - v - w)/4.0d0
       EdgeBasis(1,2) = 0.0d0
       EdgeBasis(1,3) = (u*(-1.0d0 + v + w))/(4.0d0*(-1.0d0 + w))
       CurlBasis(1,1) = u/(4.0d0*(-1.0d0 + w))
       CurlBasis(1,2) = -(-2.0d0 + v + 2.0d0*w)/(4.0d0*(-1.0d0 + w))
       CurlBasis(1,3) = 0.25d0

       EdgeBasis(2,1) = 0.0d0
       EdgeBasis(2,2) = (1.0d0 + u - w)/4.0d0
       EdgeBasis(2,3) = (v*(1.0d0 + u - w))/(4.0d0 - 4.0d0*w)
       CurlBasis(2,1) = (2.0d0 + u - 2.0d0*w)/(4.0d0 - 4.0d0*w)
       CurlBasis(2,2) = v/(4.0d0*(-1.0d0 + w))
       CurlBasis(2,3) = 0.25d0       

       EdgeBasis(3,1) = (1.0d0 + v - w)/4.0d0
       EdgeBasis(3,2) = 0.0d0
       EdgeBasis(3,3) = (u*(1.0d0 + v - w))/(4.0d0 - 4.0d0*w)
       CurlBasis(3,1) = u/(4.0d0 - 4.0d0*w)
       CurlBasis(3,2) = (2.0d0 + v - 2.0d0*w)/(4.0d0*(-1.0d0 + w))
       CurlBasis(3,3) = -0.25d0

       EdgeBasis(4,1) = 0.0d0
       EdgeBasis(4,2) = (1.0d0 - u - w)/4.0d0
       EdgeBasis(4,3) = (v*(-1.0d0 + u + w))/(4.0d0*(-1.0d0 + w))
       CurlBasis(4,1) = (-2.0d0 + u + 2.0d0*w)/(4.0d0*(-1.0d0 + w))
       CurlBasis(4,2) = v/(4.0d0 - 4.0d0*w)
       CurlBasis(4,3) = -0.25d0

       EdgeBasis(5,1) = (w*(-1.0d0 + v + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(5,2) = (w*(-1.0d0 + u + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(5,3) = (-((-1.0d0 + v)*(-1.0d0 + w)**2) + u*(v - (-1.0d0 + w)**2 - 2.0d0*v*w))/&
            (4.0d0*(-1.0d0 + w)**2)
       CurlBasis(5,1) = -(-1.0d0 + u + w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(5,2) = (-1.0d0 + v + w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(5,3) = 0.0d0

       EdgeBasis(6,1) = -(w*(-1.0d0 + v + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(6,2) = (w*(-1.0d0 - u + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(6,3) = (-((-1.0d0 + v)*(-1.0d0 + w)**2) + u*((-1.0d0 + w)**2 + v*(-1.0d0 + 2.0d0*w)))/&
            (4.0d0*(-1.0d0 + w)**2)
       CurlBasis(6,1) = (1.0d0 + u - w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(6,2) = -(-1.0d0 + v + w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(6,3) = 0.0d0    

       EdgeBasis(7,1) = ((1.0d0 + v - w)*w)/(4.0d0*(-1.0d0 + w))
       EdgeBasis(7,2) = ((1.0d0 + u - w)*w)/(4.0d0*(-1.0d0 + w))
       EdgeBasis(7,3) = ((1.0d0 + v)*(-1.0d0 + w)**2 + u*(v + (-1.0d0 + w)**2 - 2.0d0*v*w))/&
            (4.0d0*(-1.0d0 + w)**2)
       CurlBasis(7,1) = (1.0d0 + u - w)/(2.0d0 - 2.0d0*w)
       CurlBasis(7,2) = (1.0d0 + v - w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(7,3) = 0.0d0

       EdgeBasis(8,1) = (w*(-1.0d0 - v + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(8,2) = -(w*(-1.0d0 + u + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(8,3) = ((1.0d0 + v)*(-1.0d0 + w)**2 - u*(v + (-1.0d0 + w)**2 - 2.0d0*v*w))/&
            (4.0d0*(-1.0d0 + w)**2)
       CurlBasis(8,1) = (-1.0d0 + u + w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(8,2) = (1.0d0 + v - w)/(2.0d0 - 2.0d0*w)
       CurlBasis(8,3) = 0.0d0

     END SELECT

     EdgeMap => LGetEdgeMap(Element % TYPE % ElementCode / 100)
     DO i=1,SIZE(Edgemap,1)
       j = EdgeMap(i,1); k = EdgeMap(i,2)

       nj = Element % Nodeindexes(j)
       IF (Parallel) nj=Mesh % ParallelInfo % GlobalDOFs(nj)
       nk = Element % Nodeindexes(k)
       IF (Parallel) nk=Mesh % ParallelInfo % GlobalDOFs(nk)

       SELECT CASE(Element % TYPE % ElementCode / 100)
       CASE(3,5)
         WBasis(i,:) = Basis(j)*dBasisdx(k,:) - Basis(k)*dBasisdx(j,:)

         RotWBasis(i,1) = 2.0_dp * ( dBasisdx(j,2) * dBasisdx(k,3) - &
                       dBasisdx(j,3) * dBasisdx(k,2) )
         RotWBasis(i,2) = 2.0_dp * ( dBasisdx(j,3) * dBasisdx(k,1) - &
                       dBasisdx(j,1) * dBasisdx(k,3) )
         RotWBasis(i,3) = 2.0_dp * ( dBasisdx(j,1) * dBasisdx(k,2) - &
                       dBasisdx(j,2) * dBasisdx(k,1) )

       CASE(6)
          !-----------------------------------------------------------------------
          ! Create the referential description of basis functions and their 
          ! spatial curl on the physical element via applying the Piola transform:
          !-----------------------------------------------------------------------
          DO k=1,3
             WBasis(i,k) = SUM( G(1:3,k) * EdgeBasis(i,1:3) )
          END DO
          DO k=1,3
             RotWBasis(i,k) = 1.0d0/DetF * SUM( F(k,1:3) * CurlBasis(i,1:3) )
          END DO

       CASE(7)
         SELECT CASE(i)
          CASE(1)
            j=1;k=2; Base=(1-w)/2; dBase=-dudx(3,:)/2
          CASE(2)
            j=2;k=3; Base=(1-w)/2; dBase=-dudx(3,:)/2
          CASE(3)
            j=3;k=1; Base=(1-w)/2; dBase=-dudx(3,:)/2
          CASE(4)
            j=1;k=2; Base=(1+w)/2; dBase= dudx(3,:)/2
          CASE(5)
            j=2;k=3; Base=(1+w)/2; dBase= dudx(3,:)/2
          CASE(6)
            j=3;k=1; Base=(1+w)/2; dBase= dudx(3,:)/2
          CASE(7)
            Base=triBase(1); dBase=dtriBase(1,:); du=dudx(3,:)/2
          CASE(8)
            Base=triBase(2); dBase=dtriBase(2,:); du=dudx(3,:)/2
          CASE(9)
            Base=triBase(3); dBase=dtriBase(3,:); du=dudx(3,:)/2
         END SELECT

         IF(i<=6) THEN
            tBase = (triBase(j)*dtriBase(k,:)-triBase(k)*dtriBase(j,:))
            rBase(1) = 2*Base*(dtriBase(j,2)*dtriBase(k,3)-dtriBase(k,2)*dtriBase(j,3)) + &
                              dBase(2)*tBase(3) - dBase(3)*tBase(2)

            rBase(2) = 2*Base*(dtriBase(j,3)*dtriBase(k,1)-dtriBase(k,3)*dtriBase(j,1)) + &
                              dBase(3)*tBase(1) - dBase(1)*tBase(3)

            rBase(3) = 2*Base*(dtriBase(j,1)*dtriBase(k,2)-dtriBase(k,1)*dtriBase(j,2)) + &
                              dBase(1)*tBase(2) - dBase(2)*tBase(1)

            RotWBasis(i,:)=rBase
            WBasis(i,:)=tBase*Base
         ELSE
            WBasis(i,:)=Base*du
            RotWBasis(i,1)=(dBase(2)*du(3) - dBase(3)*du(2))
            RotWBasis(i,2)=(dBase(3)*du(1) - dBase(1)*du(3))
            RotWBasis(i,3)=(dBase(1)*du(2) - dBase(2)*du(1))
         END IF
       CASE(4)
         SELECT CASE(i)
          CASE(1)
             du=dudx(1,:); Base=(1-v)*(1-w)
             dBase(:)=-dudx(2,:)*(1-w)-(1-v)*dudx(3,:)
          CASE(2)
             du=dudx(2,:); Base=(1+u)*(1-w)
             dBase(:)= dudx(1,:)*(1-w)-(1+u)*dudx(3,:)
          CASE(3)
             du=-dudx(1,:); Base=(1+v)*(1-w)
             dBase(:)= dudx(2,:)*(1-w)-(1+v)*dudx(3,:)
          CASE(4)
             du=-dudx(2,:); Base=(1-u)*(1-w)
             dBase(:)=-dudx(1,:)*(1-w)-(1-u)*dudx(3,:)
         END SELECT

         wBasis(i,:) = Base*du/n
         RotWBasis(i,1)=(dBase(2)*du(3) - dBase(3)*du(2))/n
         RotWBasis(i,2)=(dBase(3)*du(1) - dBase(1)*du(3))/n
         RotWBasis(i,3) = (dBase(1)*du(2) - dBase(2)*du(1))/n
       CASE(8)
         SELECT CASE(i)
          CASE(1)
             du=dudx(1,:); Base=(1-v)*(1-w)
             dBase(:)=-dudx(2,:)*(1-w)-(1-v)*dudx(3,:)
          CASE(2)
             du=dudx(2,:); Base=(1+u)*(1-w)
             dBase(:)= dudx(1,:)*(1-w)-(1+u)*dudx(3,:)
          CASE(3)
             du=dudx(1,:); Base=(1+v)*(1-w)
             dBase(:)= dudx(2,:)*(1-w)-(1+v)*dudx(3,:)
          CASE(4)
             du=dudx(2,:); Base=(1-u)*(1-w)
             dBase(:)=-dudx(1,:)*(1-w)-(1-u)*dudx(3,:)
          CASE(5)
             du=dudx(1,:); Base=(1-v)*(1+w)
             dBase(:)=-dudx(2,:)*(1+w)+(1-v)*dudx(3,:)
          CASE(6)
             du=dudx(2,:); Base=(1+u)*(1+w)
             dBase(:)= dudx(1,:)*(1+w)+(1+u)*dudx(3,:)
          CASE(7)
             du=dudx(1,:); Base=(1+v)*(1+w)
             dBase(:)= dudx(2,:)*(1+w)+(1+v)*dudx(3,:)
          CASE(8)
             du=dudx(2,:); Base=(1-u)*(1+w)
             dBase(:)=-dudx(1,:)*(1+w)+(1-u)*dudx(3,:)
          CASE(9)
             du=dudx(3,:); Base=(1-u)*(1-v)
             dBase(:)=-dudx(1,:)*(1-v)-(1-u)*dudx(2,:)
          CASE(10)
             du=dudx(3,:); Base=(1+u)*(1-v)
             dBase(:)= dudx(1,:)*(1-v)-(1+u)*dudx(2,:)
          CASE(11)
             du=dudx(3,:); Base=(1+u)*(1+v)
             dBase(:)= dudx(1,:)*(1+v)+(1+u)*dudx(2,:)
          CASE(12)
             du=dudx(3,:); Base=(1-u)*(1+v)
             dBase(:)=-dudx(1,:)*(1+v)+(1-u)*dudx(2,:)
         END SELECT

         wBasis(i,:)=Base*du/n
         RotWBasis(i,1)=(dBase(2)*du(3) - dBase(3)*du(2))/n
         RotWBasis(i,2)=(dBase(3)*du(1) - dBase(1)*du(3))/n
         RotWBasis(i,3)=(dBase(1)*du(2) - dBase(2)*du(1))/n
       CASE DEFAULT
         CALL Fatal( 'Edge Basis', 'Not implemented for this element type.')
       END SELECT

       IF( nk < nj ) THEN
         WBasis(i,:) = -WBasis(i,:); RotWBasis(i,:) = -RotWBasis(i,:)
       END IF
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE GetEdgeBasis
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Compute contravariant metric tensor (=J^TJ)^-1 of element coordinate
!>    system, and square root of determinant of covariant metric tensor
!>    (=sqrt(det(J^TJ)))
!------------------------------------------------------------------------------
   FUNCTION ElementMetric(nDOFs,Elm,Nodes,Metric,DetG,dLBasisdx,LtoGMap) RESULT(Success)
!------------------------------------------------------------------------------
     INTEGER :: nDOFs                !< Number of active nodes in element
     TYPE(Element_t)  :: Elm         !< Element structure
     TYPE(Nodes_t)    :: Nodes       !< Element nodal coordinates
     REAL(KIND=dp) :: Metric(:,:)    !< Contravariant metric tensor
     REAL(KIND=dp) :: dLBasisdx(:,:) !< Derivatives of element basis function with respect to local coordinates
     REAL(KIND=dp) :: DetG           !< SQRT of determinant of metric tensor
     REAL(KIND=dp) :: LtoGMap(3,3)   !< Transformation to obtain the referencial description of the spatial gradient
     LOGICAL :: Success              !< Returns .FALSE. if element is degenerate
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: dx(3,3),G(3,3),GI(3,3),s
     REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
     INTEGER :: GeomId
     
     INTEGER :: cdim,dim,i,j,k,n
!------------------------------------------------------------------------------
     success = .TRUE.

     x => Nodes % x
     y => Nodes % y
     z => Nodes % z

     cdim = CoordinateSystemDimension()
     n = MIN( SIZE(x), nDOFs )
     dim  = elm % TYPE % DIMENSION

!------------------------------------------------------------------------------
!    Partial derivatives of global coordinates with respect to local coordinates
!------------------------------------------------------------------------------
     DO i=1,dim
       dx(1,i) = SUM( x(1:n) * dLBasisdx(1:n,i) )
       dx(2,i) = SUM( y(1:n) * dLBasisdx(1:n,i) )
       dx(3,i) = SUM( z(1:n) * dLBasisdx(1:n,i) )
     END DO
!------------------------------------------------------------------------------
!    Compute the covariant metric tensor of the element coordinate system
!------------------------------------------------------------------------------
     DO i=1,dim
        DO j=1,dim
           s = 0.0d0
           DO k=1,cdim
             s = s + dx(k,i)*dx(k,j)
           END DO
           G(i,j) = s
        END DO
     END DO
!------------------------------------------------------------------------------
!    Convert the metric to contravariant base, and compute the SQRT(DetG)
!------------------------------------------------------------------------------
     SELECT CASE( dim )
!------------------------------------------------------------------------------
!      Line elements
!------------------------------------------------------------------------------
       CASE (1)
         DetG  = G(1,1)

         IF ( DetG <= TINY( DetG ) ) GOTO 100

         Metric(1,1) = 1.0d0 / DetG
         DetG  = SQRT( DetG )

!------------------------------------------------------------------------------
!      Surface elements
!------------------------------------------------------------------------------
       CASE (2)
         DetG = ( G(1,1)*G(2,2) - G(1,2)*G(2,1) )

         IF ( DetG <= TINY( DetG ) ) GOTO 100

         Metric(1,1) =  G(2,2) / DetG
         Metric(1,2) = -G(1,2) / DetG
         Metric(2,1) = -G(2,1) / DetG
         Metric(2,2) =  G(1,1) / DetG
         DetG = SQRT(DetG)

!------------------------------------------------------------------------------
!      Volume elements
!------------------------------------------------------------------------------
       CASE (3)
         DetG = G(1,1) * ( G(2,2)*G(3,3) - G(2,3)*G(3,2) ) + &
                G(1,2) * ( G(2,3)*G(3,1) - G(2,1)*G(3,3) ) + &
                G(1,3) * ( G(2,1)*G(3,2) - G(2,2)*G(3,1) )

         IF ( DetG <= TINY( DetG ) ) GOTO 100

         CALL InvertMatrix3x3( G,GI,detG )
         Metric = GI
         DetG = SQRT(DetG)
     END SELECT

!--------------------------------------------------------------------------------------
!    Construct a transformation X = LtoGMap such that (grad B)(f(p)) = X(p) Grad b(p),
!    with Grad the gradient with respect to the reference element coordinates p and 
!    the referencial description of the spatial field B(x) satisfying B(f(p)) = b(p).
!    If cdim > dim (e.g. a surface embedded in the 3-dimensional space), X is
!    the pseudo-inverse of (Grad f)^{T}.
!-------------------------------------------------------------------------------
     DO i=1,cdim
       DO j=1,dim
         s = 0.0d0
         DO k=1,dim
           s = s + dx(i,k) * Metric(k,j)
         END DO
         LtoGMap(i,j) = s
       END DO
     END DO

! Return here also implies success = .TRUE.
     RETURN
  

100  Success = .FALSE.
     WRITE( Message,'(A,I0,A,I0)') 'Degenerate ',dim,'D element: ',Elm % ElementIndex
     CALL Error( 'ElementMetric', Message )
     
     IF( ASSOCIATED( Elm % BoundaryInfo ) ) THEN
       WRITE( Message,'(A,I0,A,ES12.3)') 'Boundary Id: ',Elm % BoundaryInfo % Constraint,' DetG:',DetG
     ELSE
       WRITE( Message,'(A,I0,A,ES12.3)') 'Body Id: ',Elm % BodyId,' DetG:',DetG
     END IF
     CALL Info( 'ElementMetric', Message, Level=3 )

     DO i=1,n
       WRITE( Message,'(A,I0,A,3ES12.3)') 'Node: ',i,' Coord:',x(i),y(i),z(i)       
       CALL Info( 'ElementMetric', Message, Level=3 )
     END DO
     DO i=2,n
       WRITE( Message,'(A,I0,A,3ES12.3)') 'Node: ',i,' dCoord:',&
           x(i)-x(1),y(i)-y(1),z(i)-z(1)       
       CALL Info( 'ElementMetric', Message, Level=3 )
     END DO
     IF ( cdim < dim ) THEN
       WRITE( Message,'(A,I0,A,I0)') 'Element dim larger than meshdim: ',dim,' vs. ',cdim
       CALL Info( 'ElementMetric', Message, Level=3 )
     END IF
     
!------------------------------------------------------------------------------
   END FUNCTION ElementMetric
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION ElementMetricVec( Elm, Nodes, nc, ndof, DetJ, nbmax, dLBasisdx, LtoGMap) RESULT(AllSuccess)
!------------------------------------------------------------------------------
     TYPE(Element_t)  :: Elm                                 !< Element structure
     TYPE(Nodes_t)    :: Nodes                               !< element nodal coordinates
     INTEGER, INTENT(IN) :: nc                               !< Number of points to map
     INTEGER :: ndof                                         !< Number of active nodes in element
     REAL(KIND=dp) :: DetJ(VECTOR_BLOCK_LENGTH)              !< SQRT of determinant of element coordinate metric at each point
     INTEGER, INTENT(IN) :: nbmax                            !< Maximum total number of basis functions in local basis
     REAL(KIND=dp) :: dLBasisdx(VECTOR_BLOCK_LENGTH,nbmax,3) !< Derivatives of element basis function with 
                                                             !<  respect to local coordinates at each point
     REAL(KIND=dp) :: LtoGMap(VECTOR_BLOCK_LENGTH,3,3)       !< Mapping between local and global coordinates
     LOGICAL :: AllSuccess                  !< Returns .FALSE. if some point in element is degenerate
!------------------------------------------------------------------------------
!       Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: dx(VECTOR_BLOCK_LENGTH,3,3)
     REAL(KIND=dp) :: Metric(VECTOR_BLOCK_LENGTH,6), &
             G(VECTOR_BLOCK_LENGTH,6)       ! Symmetric Metric(nc,3,3) and G(nc,3,3)

     REAL(KIND=dp) :: s
     INTEGER :: cdim,dim,i,j,k,l,n,ip, jj, kk
     INTEGER :: ldbasis, ldxyz, utind
!DIR$ ATTRIBUTES ALIGN:64::Metric
!DIR$ ATTRIBUTES ALIGN:64::dx
!DIR$ ATTRIBUTES ALIGN:64::G
!DIR$ ASSUME_ALIGNED dLBasisdx:64, LtoGMap:64, DetJ:64
     !------------------------------------------------------------------------------
     AllSuccess = .TRUE.

     ! Coordinates (single array)
     n = MIN( SIZE(Nodes % x, 1), ndof )

     ! Dimensions (coordinate system and element)
     cdim = CoordinateSystemDimension()
     dim  = elm % TYPE % DIMENSION

     ! Leading dimensions for local basis and coordinate arrays
     ldbasis = SIZE(dLBasisdx, 1)
     ldxyz = SIZE(Nodes % xyz, 1)

     ! For linear, extruded and otherwise regular elements mapping has to be computed
     ! only once, the problem is to identify these cases...
     !------------------------------------------------------------------------------
     !       Partial derivatives of global coordinates with respect to local coordinates
     !------------------------------------------------------------------------------
     ! Avoid DGEMM calls for nc small
     IF (nc < VECTOR_SMALL_THRESH) THEN
       DO l=1,dim
         DO j=1,3
           dx(1:nc,j,l)=REAL(0,dp)
           DO k=1,n
!DIR$ UNROLL
             DO i=1,nc
               dx(i,j,l)=dx(i,j,l)+dLBasisdx(i,k,l)*Nodes % xyz(k,j)
             END DO
           END DO
         END DO
       END DO
     ELSE
       DO i=1,dim
         CALL DGEMM('N','N',nc, 3, n, &
                 REAL(1,dp), dLbasisdx(1,1,i), ldbasis, &
                 Nodes % xyz, ldxyz, REAL(0, dp), dx(1,1,i), VECTOR_BLOCK_LENGTH)
       END DO
     END IF
     !------------------------------------------------------------------------------
     !       Compute the covariant metric tensor of the element coordinate system (symmetric)
     !------------------------------------------------------------------------------
     ! Linearized upper triangular indices for accesses to G
     ! | (1,1) (1,2) (1,3) | = | 1 2 4 |
     ! |       (2,2) (2,3) |   |   3 5 |
     ! |             (3,3) |   |     6 |
     ! G is symmetric, compute only the upper triangular part of G=dx^Tdx
!DIR$ LOOP COUNT MAX=3
     DO j=1,dim
!DIR$ LOOP COUNT MAX=3
       DO i=1,j
!DIR$ INLINE
         utind = GetSymmetricIndex(i,j)
         SELECT CASE (cdim)
         CASE(1)
           !_ELMER_OMP_SIMD
           DO l=1,nc
             G(l,utind)=dx(l,1,i)*dx(l,1,j)
           END DO
         CASE(2)
           !_ELMER_OMP_SIMD
           DO l=1,nc
             G(l,utind)=dx(l,1,i)*dx(l,1,j)+dx(l,2,i)*dx(l,2,j)
           END DO
         CASE(3)
           !_ELMER_OMP_SIMD
           DO l=1,nc
             G(l,utind)=dx(l,1,i)*dx(l,1,j)+dx(l,2,i)*dx(l,2,j)+dx(l,3,i)*dx(l,3,j)
           END DO
         END SELECT
       END DO
     END DO

     !------------------------------------------------------------------------------
     !       Convert the metric to contravariant base, and compute the SQRT(DetG)
     !------------------------------------------------------------------------------
     SELECT CASE( dim )
       !------------------------------------------------------------------------------
       !       Line elements
       !------------------------------------------------------------------------------
     CASE (1)
       ! Determinants
       ! DetJ(1:nc)  = G(1:nc,1,1)
       DetJ(1:nc)  = G(1:nc,1)

       DO i=1,nc
         IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
           AllSuccess = .FALSE.
           EXIT
         END IF
       END DO

       IF (AllSuccess) THEN
         !_ELMER_OMP_SIMD
         DO i=1,nc
           ! Metric(i,1,1) = REAL(1,dp)/DetJ(i)
           Metric(i,1) = REAL(1,dp)/DetJ(i)
         END DO
         !_ELMER_OMP_SIMD
         DO i=1,nc
           DetJ(i) = SQRT( DetJ(i))
         END DO
       END IF


       !------------------------------------------------------------------------------
       !       Surface elements
       !------------------------------------------------------------------------------
     CASE (2)
       ! Determinants
       !_ELMER_OMP_SIMD
       DO i=1,nc
         ! DetJ(i) = ( G(i,1,1)*G(i,2,2) - G(i,1,2)*G(i,2,1) )
         ! G is symmetric
         DetJ(i) = G(i,1)*G(i,3)-G(i,2)*G(i,2)
       END DO

       DO i=1,nc
         IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
           AllSuccess = .FALSE.
           EXIT
         END IF
       END DO

       IF (AllSuccess) THEN
         ! Since G=G^T, it holds G^{-1}=(G^T)^{-1}
         !_ELMER_OMP_SIMD
         DO i=1,nc
           s = REAL(1,dp)/DetJ(i)
           ! G is symmetric
           ! All in one go, with redundancies eliminated
           Metric(i,1) =  s*G(i,3)
           Metric(i,2) = -s*G(i,2)
           Metric(i,3) =  s*G(i,1)
         END DO
         !_ELMER_OMP_SIMD
         DO i=1,nc
           DetJ(i) = SQRT(DetJ(i))
         END DO

       END IF
       !------------------------------------------------------------------------------
       !       Volume elements
       !------------------------------------------------------------------------------
     CASE (3)
       ! Determinants
       !_ELMER_OMP_SIMD
       DO i=1,nc
         ! DetJ(i) = G(i,1,1) * ( G(i,2,2)*G(i,3,3) - G(i,2,3)*G(i,3,2) ) + &
         !           G(i,1,2) * ( G(i,2,3)*G(i,3,1) - G(i,2,1)*G(i,3,3) ) + &
         !           G(i,1,3) * ( G(i,2,1)*G(i,3,2) - G(i,2,2)*G(i,3,1) )
         ! G is symmetric
         DetJ(i) = G(i,1)*(G(i,3)*G(i,6)-G(i,5)*G(i,5)) + &
                 G(i,2)*(G(i,5)*G(i,4)-G(i,2)*G(i,6)) + &
                 G(i,4)*(G(i,2)*G(i,5)-G(i,3)*G(i,4))
       END DO

       DO i=1,nc
         IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
           AllSuccess = .FALSE.
           EXIT
         END IF
       END DO

       IF (AllSuccess) THEN
         ! Since G=G^T, it holds G^{-1}=(G^T)^{-1}
         !_ELMER_OMP_SIMD
         DO i=1,nc
           s = REAL(1,dp) / DetJ(i)
           ! Metric(i,1,1) =  s * (G(i,2,2)*G(i,3,3) - G(i,3,2)*G(i,2,3))
           ! Metric(i,2,1) = -s * (G(i,2,1)*G(i,3,3) - G(i,3,1)*G(i,2,3))
           ! Metric(i,3,1) =  s * (G(i,2,1)*G(i,3,2) - G(i,3,1)*G(i,2,2))
           ! G is symmetric

           ! All in one go, with redundancies eliminated
           Metric(i,1)= s*(G(i,3)*G(i,6)-G(i,5)*G(i,5))
           Metric(i,2)=-s*(G(i,2)*G(i,6)-G(i,4)*G(i,5))
           Metric(i,3)= s*(G(i,1)*G(i,6)-G(i,4)*G(i,4))
           Metric(i,4)= s*(G(i,2)*G(i,5)-G(i,3)*G(i,4))
           Metric(i,5)=-s*(G(i,1)*G(i,5)-G(i,2)*G(i,4))
           Metric(i,6)= s*(G(i,1)*G(i,3)-G(i,2)*G(i,2))
         END DO

         !_ELMER_OMP_SIMD
         DO i=1,nc
           DetJ(i) = SQRT(DetJ(i))
         END DO

       END IF
     END SELECT

     IF (AllSuccess) THEN
       SELECT CASE(dim)
       CASE(1)
!DIR$ LOOP COUNT MAX=3
         DO i=1,cdim
           !_ELMER_OMP_SIMD
           DO l=1,nc
             LtoGMap(l,i,1) = dx(l,i,1)*Metric(l,1)
           END DO
         END DO
       CASE(2)
!DIR$ LOOP COUNT MAX=3
         DO i=1,cdim
           !_ELMER_OMP_SIMD
           DO l=1,nc
             LtoGMap(l,i,1) = dx(l,i,1)*Metric(l,1) + dx(l,i,2)*Metric(l,2)
             LtoGMap(l,i,2) = dx(l,i,1)*Metric(l,2) + dx(l,i,2)*Metric(l,3)
           END DO
         END DO
       CASE(3)
!DIR$ LOOP COUNT MAX=3
         DO i=1,cdim
           !_ELMER_OMP_SIMD
           DO l=1,nc
             LtoGMap(l,i,1) = dx(l,i,1)*Metric(l,1) + dx(l,i,2)*Metric(l,2) + dx(l,i,3)*Metric(l,4)
             LtoGMap(l,i,2) = dx(l,i,1)*Metric(l,2) + dx(l,i,2)*Metric(l,3) + dx(l,i,3)*Metric(l,5)
             LtoGMap(l,i,3) = dx(l,i,1)*Metric(l,4) + dx(l,i,2)*Metric(l,5) + dx(l,i,3)*Metric(l,6)
           END DO
         END DO
       END SELECT
     ELSE

       ! Degenerate element!
       WRITE( Message,'(A,I0,A,I0,A,I0)') 'Degenerate ',dim,'D element: ',Elm % ElementIndex, ', pt=', i
       CALL Error( 'ElementMetricVec', Message )
       WRITE( Message,'(A,G10.3)') 'DetG:',DetJ(i)
       CALL Info( 'ElementMetricVec', Message, Level=3 )
       DO i=1,cdim
         WRITE( Message,'(A,I0,A,3G10.3)') 'Dir: ',i,' Coord:',Nodes % xyz(i,1),&
                 Nodes % xyz(i,2), Nodes % xyz(i,3)
         CALL Info( 'ElementMetricVec', Message, Level=3 )
       END DO
       IF (cdim < dim) THEN
         WRITE( Message,'(A,I0,A,I0)') 'Element dim larger than meshdim: ',dim,' vs. ',cdim
         CALL Info( 'ElementMetricVec', Message, Level=3 )
       END IF
     END IF

   CONTAINS

     FUNCTION GetSymmetricIndex(i,j) RESULT(utind)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: i, j
       INTEGER :: utind

       IF (i>j) THEN
         utind = i*(i-1)/2+j
       ELSE
         utind = j*(j-1)/2+i
       END IF
     END FUNCTION GetSymmetricIndex
!------------------------------------------------------------------------------
   END FUNCTION ElementMetricVec
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Given element structure return value of the first partial derivatives with
!>    respect to global coordinates of a quantity x given at element nodes at
!>    local coordinate point u,v,w inside the element. Element basis functions
!>    are used to compute the value. This is internal version,and shoudnt
!>    usually be called directly by the user, but trough the wrapper routine
!>    GlobalFirstDerivatives.
!------------------------------------------------------------------------------
   SUBROUTINE GlobalFirstDerivativesInternal( elm,nodes,df,gx,gy,gz, &
                       Metric,dLBasisdx )
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: element
!      INPUT: element structure
!
!    Type(Nodes_t) :: nodes
!      INPUT: element nodal coordinate arrays
!     
!     REAL(KIND=dp) :: f(:)
!      INPUT: Nodal values of the quantity whose partial derivative we want to know
!
!     REAL(KIND=dp) :: gx = @f(u,v)/@x, gy = @f(u,v)/@y, gz = @f(u,v)/@z
!      OUTPUT: Values of the partial derivatives
!
!     REAL(KIND=dp) :: Metric(:,:)
!      INPUT: Contravariant metric tensor of the element coordinate system
!
!     REAL(KIND=dp), OPTIONAL :: dLBasisdx(:,:)
!      INPUT: Values of partial derivatives with respect to local coordinates
!
!   FUNCTION VALUE:
!      .TRUE. if element is ok, .FALSE. if degenerated
!
!------------------------------------------------------------------------------
   !
   ! Return value of first derivatives of a quantity f in global
   ! coordinates at point (u,v) in gx,gy and gz.
   !
     TYPE(Element_t) :: elm
     TYPE(Nodes_t) :: nodes
 
     REAL(KIND=dp) :: df(:),Metric(:,:)
     REAL(KIND=dp) :: gx,gy,gz
     REAL(KIND=dp) :: dLBasisdx(:,:)

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
     REAL(KIND=dp) :: dx(3,3),dfc(3),s

     INTEGER :: cdim,dim,i,j,n,NB
!------------------------------------------------------------------------------

     n    = elm % TYPE % NumberOfNodes
     dim  = elm % TYPE % DIMENSION
     cdim = CoordinateSystemDimension()

     x => nodes % x
     y => nodes % y
     z => nodes % z
!------------------------------------------------------------------------------
!    Partial derivatives of global coordinates with respect to local, and
!    partial derivatives of the quantity given, also with respect to local
!    coordinates
!------------------------------------------------------------------------------
     SELECT CASE(cdim)
       CASE(1)
         DO i=1,dim
            dx(1,i) = SUM( x(1:n)*dLBasisdx(1:n,i) )
         END DO

       CASE(2)
         DO i=1,dim
            dx(1,i) = SUM( x(1:n)*dLBasisdx(1:n,i) )
            dx(2,i) = SUM( y(1:n)*dLBasisdx(1:n,i) )
         END DO

       CASE(3)
         DO i=1,dim
            dx(1,i) = SUM( x(1:n)*dLBasisdx(1:n,i) )
            dx(2,i) = SUM( y(1:n)*dLBasisdx(1:n,i) )
            dx(3,i) = SUM( z(1:n)*dLBasisdx(1:n,i) )
         END DO
     END SELECT
!------------------------------------------------------------------------------
!    Contravariant components of partials in element coordinates
!------------------------------------------------------------------------------
     DO i=1,dim
       s = 0.0d0
       DO j=1,dim
         s = s + Metric(i,j) * df(j)
       END DO
       dfc(i) = s
     END DO
!------------------------------------------------------------------------------
!    Transform partials to space coordinates
!------------------------------------------------------------------------------
     gx = 0.0d0
     gy = 0.0d0
     gz = 0.0d0
     SELECT CASE(cdim)
       CASE(1)
         gx = SUM( dx(1,1:dim) * dfc(1:dim) )

       CASE(2)
         gx = SUM( dx(1,1:dim) * dfc(1:dim) )
         gy = SUM( dx(2,1:dim) * dfc(1:dim) )

       CASE(3)
         gx = SUM( dx(1,1:dim) * dfc(1:dim) )
         gy = SUM( dx(2,1:dim) * dfc(1:dim) )
         gz = SUM( dx(3,1:dim) * dfc(1:dim) )
     END SELECT

   END SUBROUTINE GlobalFirstDerivativesInternal
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of the first partial derivative with
!>   respect to global coordinates of a quantity f given at element nodes at
!>   local coordinate point u,v,w inside the element. Element basis functions
!>   are used to compute the value.
!------------------------------------------------------------------------------
   SUBROUTINE GlobalFirstDerivatives( Elm, Nodes, df, gx, gy, gz, &
                    Metric, dLBasisdx )
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!
!   Type(Nodes_t) :: nodes
!     INPUT: element nodal coordinate arrays
!     
!   REAL(KIND=dp) :: f(:)
!     INPUT: Nodal values of the quantity whose partial derivatives we want
!            to know
!
!   REAL(KIND=dp) :: gx=@f(u,v,w)/@x, gy=@f(u,v,w)/@y, gz=@f(u,v,w)/@z
!     OUTPUT: Values of the partial derivatives
!
!   REAL(KIND=dp) :: u,v,w
!     INPUT: Point at which to evaluate the partial derivative
!
!   REAL(KIND=dp)L :: dLBasisdx(:,:)
!     INPUT: Values of partial derivatives of basis functions with respect to
!            local coordinates
!
!   REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)
!     INPUT: Values of partial derivatives of basis functions with respect to
!            global coordinates can be given here, if known, otherwise they
!            will be computed from the element basis functions.
!
!------------------------------------------------------------------------------

     TYPE(Element_t) :: elm
     TYPE(Nodes_t) :: nodes

     REAL(KIND=dp) :: gx,gy,gz
     REAL(KIND=dp) :: dLBasisdx(:,:),Metric(:,:),df(:)

!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: n
!------------------------------------------------------------------------------

    CALL GlobalFirstDerivativesInternal( Elm, Nodes, df, &
              gx, gy, gz, Metric, dLBasisdx )

   END SUBROUTINE GlobalFirstDerivatives
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of a quantity x given at element nodes
!>   at local coordinate point u inside the element. Element basis functions are
!>   used to compute the value. This is just a wrapper routine and will call the
!>   real function according to element dimension.   
!------------------------------------------------------------------------------
   FUNCTION InterpolateInElement( elm,f,u,v,w,Basis ) RESULT(VALUE)
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
     REAL(KIND=dp) :: VALUE
     INTEGER :: n

     IF ( PRESENT( Basis ) ) THEN
!------------------------------------------------------------------------------
!      Basis function values given, just sum the result ...
!------------------------------------------------------------------------------
       n = elm % TYPE % NumberOfNodes
       VALUE = SUM( f(1:n)*Basis(1:n) )
     ELSE
!------------------------------------------------------------------------------
!      ... otherwise compute from the definition.
!------------------------------------------------------------------------------
       SELECT CASE (elm % TYPE % DIMENSION)
         CASE (0)
           VALUE = f(1)
         CASE (1)
           VALUE = InterpolateInElement1D( elm,f,u )
         CASE (2)
           VALUE = InterpolateInElement2D( elm,f,u,v )
         CASE (3)
           VALUE = InterpolateInElement3D( elm,f,u,v,w )
       END SELECT
     END IF
  
   END FUNCTION InterpolateInElement
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>          Compute elementwise matrix of second partial derivatives
!>          at given point u,v,w in global coordinates.
!------------------------------------------------------------------------------
   SUBROUTINE GlobalSecondDerivatives(elm,nodes,f,values,u,v,w,Metric,dBasisdx)
!------------------------------------------------------------------------------
!  
!       Parameters:
!  
!           Input:   (Element_t) structure describing the element
!                    (Nodes_t)   element nodal coordinates
!                    (double precision) F nodal values of the quantity
!                    (double precision) u,v point at which to evaluate
!  
!           Output:   3x3 matrix (values) of partial derivatives
!  
!------------------------------------------------------------------------------

     TYPE(Nodes_t)   :: nodes
     TYPE(Element_t) :: elm
 
     REAL(KIND=dp) :: u,v,w
     REAL(KIND=dp) ::  f(:),Metric(:,:)
     REAL(KIND=dp) ::  values(:,:)
     REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,l,dim,cdim

     REAL(KIND=dp), DIMENSION(3,3,3) :: C1,C2,ddx
     REAL(KIND=dp), DIMENSION(3)     :: df
     REAL(KIND=dp), DIMENSION(3,3)   :: cddf,ddf,dx

     REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
     REAL(KIND=dp) :: s

     INTEGER :: n
!------------------------------------------------------------------------------
#if 1
!
! This is actually not quite correct...
!
     IF ( elm % TYPE % BasisFunctionDegree <= 1 ) RETURN
#else
!
! this is ...
!
     IF ( elm % TYPE % ElementCode <= 202 .OR. &
          elm % TYPE % ElementCode == 303 .OR. &
          elm % TYPE % ElementCode == 504 ) RETURN
#endif

     n  = elm % TYPE % NumberOfNodes
     x => nodes % x
     y => nodes % y
     z => nodes % z

     dim  = elm % TYPE % DIMENSION
     cdim = CoordinateSystemDimension()

!------------------------------------------------------------------------------
!    Partial derivatives of the basis functions are given, just
!    sum for the first partial derivatives...
!------------------------------------------------------------------------------
     dx = 0.0d0
     df = 0.0d0
     SELECT CASE( cdim )
       CASE(1)
         DO i=1,dim
           dx(1,i) = SUM( x(1:n)*dBasisdx(1:n,i) )
           df(i)   = SUM( f(1:n)*dBasisdx(1:n,i) )
         END DO

       CASE(2)
         DO i=1,dim
           dx(1,i) = SUM( x(1:n)*dBasisdx(1:n,i) )
           dx(2,i) = SUM( y(1:n)*dBasisdx(1:n,i) )
           df(i)   = SUM( f(1:n)*dBasisdx(1:n,i) )
         END DO

       CASE(3)
         DO i=1,dim
           dx(1,i) = SUM( x(1:n)*dBasisdx(1:n,i) )
           dx(2,i) = SUM( y(1:n)*dBasisdx(1:n,i) )
           dx(3,i) = SUM( z(1:n)*dBasisdx(1:n,i) )
           df(i)   = SUM( f(1:n)*dBasisdx(1:n,i) )
         END DO
     END SELECT
!------------------------------------------------------------------------------
!     Get second partial derivatives with respect to local coordinates
!------------------------------------------------------------------------------
     SELECT CASE( dim )
       CASE(1)
!------------------------------------------------------------------------------
!        Line elements
!------------------------------------------------------------------------------
         ddx(1,1,1) = SecondDerivatives1D( elm,x,u )
         ddx(2,1,1) = SecondDerivatives1D( elm,y,u )
         ddx(3,1,1) = SecondDerivatives1D( elm,z,u )

       CASE(2)
!------------------------------------------------------------------------------
!        Surface elements
!------------------------------------------------------------------------------
         ddx(1,1:2,1:2) = SecondDerivatives2D( elm,x,u,v )
         ddx(2,1:2,1:2) = SecondDerivatives2D( elm,y,u,v )
         ddx(3,1:2,1:2) = SecondDerivatives2D( elm,z,u,v )

       CASE(3)
!------------------------------------------------------------------------------
!        Volume elements
!------------------------------------------------------------------------------
         ddx(1,1:3,1:3) = SecondDerivatives3D( elm,x,u,v,w )
         ddx(2,1:3,1:3) = SecondDerivatives3D( elm,y,u,v,w )
         ddx(3,1:3,1:3) = SecondDerivatives3D( elm,z,u,v,w )
      END SELECT
!
!------------------------------------------------------------------------------
!    Christoffel symbols of the second kind of the element coordinate system
!------------------------------------------------------------------------------
      DO i=1,dim
        DO j=1,dim
          DO k=1,dim
            s = 0.0d0
            DO l=1,cdim
              s = s + ddx(l,i,j)*dx(l,k)
            END DO
            C2(i,j,k) = s
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
!    Christoffel symbols of the first kind
!------------------------------------------------------------------------------
      DO i=1,dim
        DO j=1,dim
          DO k=1,dim
            s = 0.0d0
            DO l=1,dim
              s = s + Metric(k,l)*C2(i,j,l)
            END DO
            C1(i,j,k) = s
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
!     First add ordinary partials (change of the quantity with coordinates)...
!------------------------------------------------------------------------------
      SELECT CASE(dim)
        CASE(1)
          ddf(1,1) = SecondDerivatives1D( elm,f,u )

        CASE(2)
          ddf(1:2,1:2) = SecondDerivatives2D( elm,f,u,v )

        CASE(3)
          ddf(1:3,1:3) = SecondDerivatives3D( elm,f,u,v,w )
      END SELECT
!------------------------------------------------------------------------------
!     ... then add change of coordinates
!------------------------------------------------------------------------------
      DO i=1,dim
        DO j=1,dim
          s = 0.0d0
          DO k=1,dim
            s = s - C1(i,j,k)*df(k)
          END DO
          ddf(i,j) = ddf(i,j) + s
        END DO
      END DO
!------------------------------------------------------------------------------
!     Convert to contravariant base
!------------------------------------------------------------------------------
      DO i=1,dim
        DO j=1,dim
          s = 0.0d0
          DO k=1,dim
            DO l=1,dim
              s = s + Metric(i,k)*Metric(j,l)*ddf(k,l)
            END DO
          END DO
          cddf(i,j) = s
        END DO
      END DO
!------------------------------------------------------------------------------
!    And finally transform to global coordinates 
!------------------------------------------------------------------------------
      Values = 0.0d0
      DO i=1,cdim
        DO j=1,cdim
          s = 0.0d0
          DO k=1,dim
            DO l=1,dim
              s = s + dx(i,k)*dx(j,l)*cddf(k,l)    
            END DO
          END DO
          Values(i,j) = s
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE GlobalSecondDerivatives
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
 FUNCTION LGetEdgeMap( ElementFamily ) RESULT(EdgeMap)
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
      CALL Fatal( 'LGetEdgeMap', Message )
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
  END FUNCTION LGetEdgeMap
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Figure out element diameter parameter for stablization.
!------------------------------------------------------------------------------
   FUNCTION ElementDiameter( elm, nodes, UseLongEdge ) RESULT(hK)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!     
!    Type(Nodes_t) :: nodes
!     INPUT: Nodal coordinate arrays of the element
!
!  FUNCTION VALUE:
!     REAL(KIND=dp) :: hK
!    
!------------------------------------------------------------------------------
     TYPE(Element_t) :: elm
     TYPE(Nodes_t) :: nodes
     LOGICAL, OPTIONAL :: UseLongEdge
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp), DIMENSION(:), POINTER :: X,Y,Z
     INTEGER :: i,j,k,Family
     INTEGER, POINTER :: EdgeMap(:,:)
     REAL(KIND=dp) :: x0,y0,z0,hK,A,S,CX,CY,CZ
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
         EdgeMap => LGetEdgeMap(Family)

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
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>     Figure out if given point x,y,z is inside a triangle, whose node
!>     coordinates are given in nx,ny,nz. Method: Invert the basis
!>     functions....
!------------------------------------------------------------------------------
  FUNCTION TriangleInside( nx,ny,nz,x,y,z ) RESULT(inside)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    REAL(KIND=dp) :: nx(:),ny(:),nz(:)
!      INPUT:  Node coordinate arrays
!
!    REAL(KIND=dp) :: x,y,z
!      INPUT: point which to consider
!
!  FUNCTION VALUE:
!    LOGICAL :: inside
!       result of the in/out test
!    
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: nx(:),ny(:),nz(:),x,y,z

!------------------------------------------------------------------------------
!   Local variables
!------------------------------------------------------------------------------
    LOGICAL :: inside

    REAL(KIND=dp) :: a00,a01,a10,a11,b00,b01,b10,b11,detA,px,py,u,v
!------------------------------------------------------------------------------

    inside = .FALSE.

    IF ( MAXVAL(nx) < x .OR. MAXVAL(ny) < y ) RETURN
    IF ( MINVAL(nx) > x .OR. MINVAL(ny) > y ) RETURN

    A00 = nx(2) - nx(1)
    A01 = nx(3) - nx(1)
    A10 = ny(2) - ny(1)
    A11 = ny(3) - ny(1)

    detA = A00*A11 - A01*A10
    IF ( ABS(detA) < AEPS ) RETURN

    detA = 1 / detA

    B00 =  A11*detA
    B01 = -A01*detA
    B10 = -A10*detA
    B11 =  A00*detA

    px = x - nx(1)
    py = y - ny(1)
    u = 0.0d0
    v = 0.0d0

    u = B00*px + B01*py
    IF ( u < 0.0d0 .OR. u > 1.0d0 ) RETURN

    v = B10*px + B11*py
    IF ( v < 0.0d0 .OR. v > 1.0d0 ) RETURN

    inside = (u + v <=  1.0d0)
!------------------------------------------------------------------------------
   END FUNCTION TriangleInside
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>     Figure out if given point x,y,z is inside a quadrilateral, whose
!>     node coordinates are given in nx,ny,nz. Method: Invert the
!>     basis functions....
!------------------------------------------------------------------------------
   FUNCTION QuadInside( nx,ny,nz,x,y,z ) RESULT(inside)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    REAL(KIND=dp) :: nx(:),ny(:),nz(:)
!      INPUT:  Node coordinate arrays
!
!    REAL(KIND=dp) :: x,y,z
!      INPUT: point which to consider
!
!  FUNCTION VALUE:
!    LOGICAL :: inside
!       result of the in/out test
!    
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: nx(:),ny(:),nz(:),x,y,z
!------------------------------------------------------------------------------
!   Local variables
!------------------------------------------------------------------------------
    LOGICAL :: inside

    REAL(KIND=dp) :: r,a,b,c,d,ax,bx,cx,dx,ay,by,cy,dy,px,py,u,v
!------------------------------------------------------------------------------
    inside = .FALSE.

    IF ( MAXVAL(nx) < x .OR. MAXVAL(ny) < y ) RETURN
    IF ( MINVAL(nx) > x .OR. MINVAL(ny) > y ) RETURN

    ax = 0.25*(  nx(1) + nx(2) + nx(3) + nx(4) )
    bx = 0.25*( -nx(1) + nx(2) + nx(3) - nx(4) )
    cx = 0.25*( -nx(1) - nx(2) + nx(3) + nx(4) )
    dx = 0.25*(  nx(1) - nx(2) + nx(3) - nx(4) )

    ay = 0.25*(  ny(1) + ny(2) + ny(3) + ny(4) )
    by = 0.25*( -ny(1) + ny(2) + ny(3) - ny(4) )
    cy = 0.25*( -ny(1) - ny(2) + ny(3) + ny(4) )
    dy = 0.25*(  ny(1) - ny(2) + ny(3) - ny(4) )

    px = x - ax
    py = y - ay

    a = cy*dx - cx*dy
    b = bx*cy - by*cx + dy*px - dx*py
    c = by*px - bx*py

    u = 0.0d0
    v = 0.0d0

    IF ( ABS(a) < AEPS ) THEN
      r = -c / b
      IF ( r < -1.0d0 .OR. r > 1.0d0 ) RETURN

      v = r
      u = (px - cx*r)/(bx + dx*r)
      inside = (u >= -1.0d0 .AND. u <= 1.0d0)
      RETURN
    END IF

    d = b*b - 4*a*c
    IF ( d < 0.0d0 ) RETURN

    d = SQRT(d)
    IF ( b>0 ) THEN
      r = -2*c/(b+d)
    ELSE
      r = (-b+d)/(2*a)
    END IF
    IF ( r >= -1.0d0 .AND. r <= 1.0d0 ) THEN
      v = r
      u = (px - cx*r)/(bx + dx*r)
        
      IF ( u >= -1.0d0 .AND. u <= 1.0d0 ) THEN
        inside = .TRUE.
        RETURN
      END IF
    END IF

    IF ( b>0 ) THEN
      r = -(b+d)/(2*a)
    ELSE
      r = 2*c/(-b+d)
    END IF
    IF ( r >= -1.0d0 .AND. r <= 1.0d0 ) THEN
      v = r
      u = (px - cx*r)/(bx + dx*r)
      inside = u >= -1.0d0 .AND. u <= 1.0d0
      RETURN
    END IF
!------------------------------------------------------------------------------
  END FUNCTION QuadInside
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>     Figure out if given point x,y,z is inside a tetrahedron, whose
!>     node coordinates are given in nx,ny,nz. Method: Invert the
!>     basis functions....
!------------------------------------------------------------------------------
  FUNCTION TetraInside( nx,ny,nz,x,y,z ) RESULT(inside)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    REAL(KIND=dp) :: nx(:),ny(:),nz(:)
!      INPUT:  Node coordinate arrays
!
!    REAL(KIND=dp) :: x,y,z
!      INPUT: point which to consider
!
!  FUNCTION VALUE:
!    LOGICAL :: inside
!       result of the in/out test
!    
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: nx(:),ny(:),nz(:),x,y,z

!------------------------------------------------------------------------------
!   Local variables
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: A00,A01,A02,A10,A11,A12,A20,A21,A22,detA
    REAL(KIND=dp) :: B00,B01,B02,B10,B11,B12,B20,B21,B22

    LOGICAL :: inside

    REAL(KIND=dp) :: px,py,pz,u,v,w
!------------------------------------------------------------------------------
    inside = .FALSE.

    IF ( MAXVAL(nx) < x .OR. MAXVAL(ny) < y .OR. MAXVAL(nz) < z ) RETURN
    IF ( MINVAL(nx) > x .OR. MINVAL(ny) > y .OR. MINVAL(nz) > z ) RETURN

    A00 = nx(2) - nx(1)
    A01 = nx(3) - nx(1)
    A02 = nx(4) - nx(1)

    A10 = ny(2) - ny(1)
    A11 = ny(3) - ny(1)
    A12 = ny(4) - ny(1)

    A20 = nz(2) - nz(1)
    A21 = nz(3) - nz(1)
    A22 = nz(4) - nz(1)

    detA =        A00*(A11*A22 - A12*A21)
    detA = detA + A01*(A12*A20 - A10*A22)
    detA = detA + A02*(A10*A21 - A11*A20)
    IF ( ABS(detA) < AEPS ) RETURN

    detA = 1 / detA

    px = x - nx(1)
    py = y - ny(1)
    pz = z - nz(1)

    B00 = (A11*A22 - A12*A21)*detA
    B01 = (A21*A02 - A01*A22)*detA
    B02 = (A01*A12 - A11*A02)*detA

    u = B00*px + B01*py + B02*pz
    IF ( u < 0.0d0 .OR. u > 1.0d0 ) RETURN


    B10 = (A12*A20 - A10*A22)*detA
    B11 = (A00*A22 - A20*A02)*detA
    B12 = (A10*A02 - A00*A12)*detA

    v = B10*px + B11*py + B12*pz
    IF ( v < 0.0d0 .OR. v > 1.0d0 ) RETURN


    B20 = (A10*A21 - A11*A20)*detA
    B21 = (A01*A20 - A00*A21)*detA
    B22 = (A00*A11 - A10*A01)*detA

    w = B20*px + B21*py + B22*pz
    IF ( w < 0.0d0 .OR. w > 1.0d0 ) RETURN

    inside = (u + v + w) <= 1.0d0
!------------------------------------------------------------------------------
  END FUNCTION TetraInside
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>     Figure out if given point x,y,z is inside a brick, whose node coordinates
!>     are given in nx,ny,nz. Method: Divide to tetrahedrons.
!------------------------------------------------------------------------------
  FUNCTION BrickInside( nx,ny,nz,x,y,z ) RESULT(inside)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    REAL(KIND=dp) :: nx(:),ny(:),nz(:)
!      INPUT:  Node coordinate arrays
!
!    REAL(KIND=dp) :: x,y,z
!      INPUT: point which to consider
!
!  FUNCTION VALUE:
!    LOGICAL :: inside
!       result of the in/out test
!    
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: nx(:),ny(:),nz(:),x,y,z

!------------------------------------------------------------------------------
!   Local variables
!------------------------------------------------------------------------------
    LOGICAL :: inside

    INTEGER :: i,j
    REAL(KIND=dp) :: px(4),py(4),pz(4),r,s,t,maxx,minx,maxy,miny,maxz,minz

    INTEGER :: map(3,12)
!------------------------------------------------------------------------------
    map = RESHAPE( [ 0,1,2,   0,2,3,   4,5,6,   4,6,7,   3,2,6,   3,6,7,  &
     1,5,6,   1,6,2,   0,4,7,   0,7,3,   0,1,5,   0,5,4 ], [ 3,12 ] ) + 1
    
    inside = .FALSE.

    IF ( MAXVAL(nx) < x .OR. MAXVAL(ny) < y .OR. MAXVAL(nz) < z ) RETURN
    IF ( MINVAL(nx) > x .OR. MINVAL(ny) > y .OR. MINVAL(nz) > z ) RETURN

    px(1) = 0.125d0 * SUM(nx)
    py(1) = 0.125d0 * SUM(ny)
    pz(1) = 0.125d0 * SUM(nz)

    DO i=1,12
      px(2:4) = nx(map(1:3,i))
      py(2:4) = ny(map(1:3,i))
      pz(2:4) = nz(map(1:3,i))

      IF ( TetraInside( px,py,pz,x,y,z ) ) THEN
        inside = .TRUE.
        RETURN
      END IF
    END DO
!------------------------------------------------------------------------------
  END FUNCTION BrickInside
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
    CHARACTER(LEN=MAX_NAME_LEN) :: PassName, PrevPassName
    LOGICAL :: NoPassiveElements = .FALSE.

    SAVE Passive, PrevPassName, NoPassiveElements
    !$OMP THREADPRIVATE(Passive, PrevPassName, NoPassiveElements)
    !------------------------------------------------------------------------------
    IsPassive = .FALSE.


    nlen = CurrentModel % Solver % Variable % NameLen
    PassName = GetVarName(CurrentModel % Solver % Variable) // ' Passive'     

    IF( PassName(1:nlen) == PrevPassName(1:nlen) ) THEN
      IF( NoPassiveElements ) RETURN
    ELSE
      NoPassiveElements = .NOT. ListCheckPresentAnyBodyForce( CurrentModel, PassName )
      PrevPassName = PassName
      IF( NoPassiveElements ) RETURN       
    END IF

    IF (PRESENT(UElement)) THEN
      Element => UElement
    ELSE
#ifdef _OPENMP
      IF (omp_in_parallel()) THEN
        CALL Fatal('CheckPassiveElement', &
             'Need an element to update inside a threaded region')
      END IF
#endif
      Element => CurrentModel % CurrentElement
    END IF

    body_id = Element % BodyId 
    IF ( body_id <= 0 )  RETURN   ! body_id == 0 for boundary elements

    bf_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
         'Body Force', Found, minv=1,maxv=CurrentModel % NumberOfBodyForces )
    IF ( .NOT. Found )  RETURN

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
!>   Normal will point from more dense material to less dense
!>   or outwards, if no elements on the other side.
!------------------------------------------------------------------------------
  SUBROUTINE CheckNormalDirection( Boundary,Normal,x,y,z,turn )
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: Boundary
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: Normal(3),x,y,z
    LOGICAL, OPTIONAL :: turn
!------------------------------------------------------------------------------

    TYPE (Element_t), POINTER :: Element,LeftElement,RightElement

    INTEGER :: LMat,RMat,n,k

    REAL(KIND=dp) :: x1,y1,z1
    REAL(KIND=dp), ALLOCATABLE :: nx(:),ny(:),nz(:)
    LOGICAL :: LPassive
!------------------------------------------------------------------------------
    
    IF(.NOT. ASSOCIATED( Boundary % BoundaryInfo ) )  RETURN
    
    k = Boundary % BoundaryInfo % OutBody

    LeftElement => Boundary % BoundaryInfo % Left

    Element => Null()
    IF ( ASSOCIATED(LeftELement) ) THEN
       RightElement => Boundary % BoundaryInfo % Right
       IF ( ASSOCIATED( RightElement ) ) THEN ! we have a body-body boundary        
         IF ( k > 0 ) THEN ! declared outbody 
           IF ( LeftElement % BodyId == k ) THEN
             Element => RightElement
           ELSE
             Element => LeftElement
           END IF
         ELSE IF (LeftElement % BodyId > RightElement % BodyId) THEN ! normal pointing into body with lower body ID
             Element => LeftElement
         ELSE IF (LeftElement % BodyId < RightElement % BodyId) THEN! normal pointing into body with lower body ID
           Element => RightElement
         ELSE ! active/passive boundary
           LPassive = CheckPassiveElement( LeftElement )
           IF (LPassive .NEQV. CheckPassiveElement( RightElement )) THEN 
             IF(LPassive) THEN
               Element => RightElement
             ELSE
               Element => LeftElement
             END IF
           END IF
         END IF
       ELSE ! body-vacuum boundary from left->right
         Element => LeftElement
       END IF
    ELSE! body-vacuum boundary from right->left
       Element => Boundary % BoundaryInfo % Right
    END IF

    IF ( .NOT. ASSOCIATED(Element) ) RETURN

    n = Element % TYPE % NumberOfNodes

    ALLOCATE( nx(n), ny(n), nz(n) )

    nx(1:n) = CurrentModel % Nodes % x(Element % NodeIndexes)
    ny(1:n) = CurrentModel % Nodes % y(Element % NodeIndexes)
    nz(1:n) = CurrentModel % Nodes % z(Element % NodeIndexes)

    SELECT CASE( Element % TYPE % ElementCode / 100 )

    CASE(2,4,8)
       x1 = InterpolateInElement( Element, nx, 0.0d0, 0.0d0, 0.0d0 )
       y1 = InterpolateInElement( Element, ny, 0.0d0, 0.0d0, 0.0d0 )
       z1 = InterpolateInElement( Element, nz, 0.0d0, 0.0d0, 0.0d0 )
    CASE(3)
       x1 = InterpolateInElement( Element, nx, 1.0d0/3, 1.0d0/3, 0.0d0 )
       y1 = InterpolateInElement( Element, ny, 1.0d0/3, 1.0d0/3, 0.0d0 )
       z1 = InterpolateInElement( Element, nz, 1.0d0/3, 1.0d0/3, 0.0d0 )
    CASE(5)
       x1 = InterpolateInElement( Element, nx, 1.0d0/4, 1.0d0/4, 1.0d0/4 )
       y1 = InterpolateInElement( Element, ny, 1.0d0/4, 1.0d0/4, 1.0d0/4 )
       z1 = InterpolateInElement( Element, nz, 1.0d0/4, 1.0d0/4, 1.0d0/4 )
    CASE(6)
       x1 = InterpolateInElement( Element, nx, 0.0d0, 0.0d0, 1.0d0/3 )
       y1 = InterpolateInElement( Element, ny, 0.0d0, 0.0d0, 1.0d0/3 )
       z1 = InterpolateInElement( Element, nz, 0.0d0, 0.0d0, 1.0d0/3 )
    CASE(7)
       x1 = InterpolateInElement( Element, nx, 1.0d0/3, 1.0d0/3, 0.0d0 )
       y1 = InterpolateInElement( Element, ny, 1.0d0/3, 1.0d0/3, 0.0d0 )
       z1 = InterpolateInElement( Element, nz, 1.0d0/3, 1.0d0/3, 0.0d0 )
    CASE DEFAULT
       CALL Fatal('CheckNormalDirection','Invalid elementcode for parent element!')   

    END SELECT
    x1 = x1 - x
    y1 = y1 - y
    z1 = z1 - z

    IF ( PRESENT(turn) ) turn = .FALSE.
    IF ( x1*Normal(1) + y1*Normal(2) + z1*Normal(3) > 0 ) THEN
       IF ( Element % BodyId /= k ) THEN
          Normal = -Normal
          IF ( PRESENT(turn) ) turn = .TRUE.
       END IF
    ELSE IF (  Element % BodyId == k ) THEN
       Normal = -Normal
       IF ( PRESENT(turn) ) turn = .TRUE.
    END IF
    DEALLOCATE( nx,ny,nz )
!------------------------------------------------------------------------------
  END SUBROUTINE CheckNormalDirection
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Gives the normal vector of a boundary element.
!> For noncurved elements the normal vector does not depend on the local coordinate
!> while otherwise it does. There are different uses of the function where some
!> do not have the luxury of knowing the local coordinates and hence the center
!> point is used as default.
!------------------------------------------------------------------------------
  FUNCTION NormalVector( Boundary,BoundaryNodes,u0,v0,Check ) RESULT(Normal)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Boundary
    TYPE(Nodes_t)   :: BoundaryNodes
    REAL(KIND=dp), OPTIONAL :: u0,v0
    LOGICAL, OPTIONAL :: Check 
    REAL(KIND=dp) :: Normal(3)
!------------------------------------------------------------------------------
    LOGICAL :: DoCheck
    TYPE(ElementType_t),POINTER :: elt
    REAL(KIND=dp) :: u,v,Auu,Auv,Avu,Avv,detA,x,y,z
    REAL(KIND=dp) :: dxdu,dxdv,dydu,dydv,dzdu,dzdv
    REAL(KIND=dp), DIMENSION(:), POINTER :: nx,ny,nz

!------------------------------------------------------------------------------

    nx => BoundaryNodes % x
    ny => BoundaryNodes % y
    nz => BoundaryNodes % z
    
    SELECT CASE ( Boundary % TYPE % DIMENSION )

    CASE ( 0 ) 
      Normal(1) = 1.0_dp
      Normal(2:3) = 0.0_dp

    CASE ( 1 ) 
      IF( PRESENT( u0 ) ) THEN
        u = u0
      ELSE
        u = 0.0_dp
      END IF

      dxdu = FirstDerivative1D( Boundary,nx,u )
      dydu = FirstDerivative1D( Boundary,ny,u )
 
      detA = dxdu*dxdu + dydu*dydu
      IF ( detA <= 0._dp ) THEN
        Normal = 0._dp
        RETURN
      END IF
      detA = 1.0_dp / SQRT(detA)
      Normal(1) = -dydu * detA
      Normal(2) =  dxdu * detA
      Normal(3) =  0.0d0
    
    CASE ( 2 ) 
      IF( PRESENT( u0 ) ) THEN
        u = u0
        v = v0
      ELSE
        IF( Boundary % TYPE % ElementCode / 100 == 3 ) THEN
          u = 1.0_dp/3
          v = 1.0_dp/3
        ELSE
          u = 0.0_dp
          v = 0.0_dp
        END IF
      END IF

      dxdu = FirstDerivativeInU2D( Boundary,nx,u,v )
      dydu = FirstDerivativeInU2D( Boundary,ny,u,v )
      dzdu = FirstDerivativeInU2D( Boundary,nz,u,v )

      dxdv = FirstDerivativeInV2D( Boundary,nx,u,v )
      dydv = FirstDerivativeInV2D( Boundary,ny,u,v )
      dzdv = FirstDerivativeInV2D( Boundary,nz,u,v )

      Auu = dxdu*dxdu + dydu*dydu + dzdu*dzdu
      Auv = dxdu*dxdv + dydu*dydv + dzdu*dzdv
      Avv = dxdv*dxdv + dydv*dydv + dzdv*dzdv

      detA = 1.0d0 / SQRT(Auu*Avv - Auv*Auv)

      Normal(1) = (dydu * dzdv - dydv * dzdu) * detA
      Normal(2) = (dxdv * dzdu - dxdu * dzdv) * detA
      Normal(3) = (dxdu * dydv - dxdv * dydu) * detA
    
    CASE DEFAULT
      CALL Fatal('NormalVector','Invalid dimension for determining normal!')
      
    END SELECT


    DoCheck = .FALSE.
    IF ( PRESENT(Check) ) DoCheck = Check

    IF ( DoCheck ) THEN
      SELECT CASE( Boundary % TYPE % ElementCode / 100 ) 
        
      CASE(1)
        x = nx(1)
        y = nx(1)
        z = nz(1)

      CASE(2,4)
        x = InterpolateInElement( Boundary,nx,0.0d0,0.0d0,0.0d0 )
        y = InterpolateInElement( Boundary,ny,0.0d0,0.0d0,0.0d0 )
        z = InterpolateInElement( Boundary,nz,0.0d0,0.0d0,0.0d0 )

      CASE(3)
        x = InterpolateInElement( Boundary,nx,1.0d0/3,1.0d0/3,0.0d0)
        y = InterpolateInElement( Boundary,ny,1.0d0/3,1.0d0/3,0.0d0)
        z = InterpolateInElement( Boundary,nz,1.0d0/3,1.0d0/3,0.0d0)
      END SELECT

      CALL CheckNormalDirection( Boundary,Normal,x,y,z )

    END IF

!------------------------------------------------------------------------------
  END FUNCTION NormalVector
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Returns a point that is most importantly supposed to be on the surface
!> For noncurved elements this may simply be the mean while otherwise
!> there may be a need to find the surface node using the local coordinates.
!> Hence the optional parameters. Typically the NormalVector and SurfaceVector
!> should be defined at the same position.
!------------------------------------------------------------------------------
  FUNCTION SurfaceVector( Boundary,BoundaryNodes,u,v ) RESULT(Surface)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Boundary
    TYPE(Nodes_t)   :: BoundaryNodes
    REAL(KIND=dp),OPTIONAL :: u,v
    REAL(KIND=dp) :: Surface(3)
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(:), POINTER :: nx,ny,nz
    INTEGER :: i,n
!------------------------------------------------------------------------------

    nx => BoundaryNodes % x
    ny => BoundaryNodes % y
    nz => BoundaryNodes % z
    n = Boundary % TYPE % NumberOfNodes

    IF( .NOT. PRESENT( u ) ) THEN
      Surface(1) = SUM( nx ) / n
      Surface(2) = SUM( ny ) / n
      Surface(3) = SUM( nz ) / n
    ELSE
      IF( Boundary % TYPE % DIMENSION == 1 ) THEN
        Surface(1) = InterpolateInElement( Boundary,nx,u,0.0_dp,0.0_dp)
        Surface(2) = InterpolateInElement( Boundary,ny,u,0.0_dp,0.0_dp)
        Surface(3) = InterpolateInElement( Boundary,nz,u,0.0_dp,0.0_dp)
      ELSE 
        Surface(1) = InterpolateInElement( Boundary,nx,u,v,0.0_dp)
        Surface(2) = InterpolateInElement( Boundary,ny,u,v,0.0_dp)
        Surface(3) = InterpolateInElement( Boundary,nz,u,v,0.0_dp)        
      END IF
    END IF

!------------------------------------------------------------------------------
  END FUNCTION SurfaceVector
!------------------------------------------------------------------------------


!---------------------------------------------------------------------------
!> This subroutine tests where the intersection between the line defined by two 
!> points and a plane (or line) defined by a boundary element meet. There is
!> an intersection if ( 0 < Lambda < 1 ). Of all intersections the first one is 
!> that with the smallest positive lambda. 
!---------------------------------------------------------------------------
  FUNCTION LineFaceIntersection(FaceElement,FaceNodes,&
      Rinit,Rfin,u,v) RESULT ( Lambda )
!---------------------------------------------------------------------------
    TYPE(Nodes_t) :: FaceNodes
    TYPE(Element_t), POINTER   :: FaceElement
    REAL(KIND=dp) :: Rinit(3),Rfin(3)
    REAL(KIND=dp),OPTIONAL :: u,v
    REAL(KIND=dp) :: Lambda

    REAL (KIND=dp) :: Surface(3),t1(3),t2(3),Normal(3),Rproj
    REAL (KIND=dp) :: Lambda0
    INTEGER :: third

    third = 3

100 CONTINUE

    ! For higher order elements this may be a necessity
    IF( PRESENT( u ) .AND. PRESENT(v) ) THEN
      Surface = SurfaceVector( FaceElement, FaceNodes, u, v )
      Normal = NormalVector( FaceElement, FaceNodes, u, v )

    ELSE IF( FaceElement % TYPE % DIMENSION == 2 ) THEN
      ! Any point known to be at the surface, even corner node
      Surface(1) = FaceNodes % x(1)
      Surface(2) = FaceNodes % y(1)
      Surface(3) = FaceNodes % z(1)

      ! Tangent vector, nor normalized to unity!
      t1(1) = FaceNodes % x(2) - Surface(1)
      t1(2) = FaceNodes % y(2) - Surface(2)
      t1(3) = FaceNodes % z(2) - Surface(3)

      t2(1) = FaceNodes % x(third) - Surface(1)
      t2(2) = FaceNodes % y(third) - Surface(2)
      t2(3) = FaceNodes % z(third) - Surface(3)

      ! Normal vector obtained from the cross product of tangent vectoes
      ! This is not normalized to unity as value of lambda does not depend on its magnitude
      Normal(1) = t1(2)*t2(3) - t1(3)*t2(2)
      Normal(2) = t1(3)*t2(1) - t1(1)*t2(3)
      Normal(3) = t1(1)*t2(2) - t1(2)*t2(1)
    ELSE
      Surface(1) = FaceNodes % x(1)
      Surface(2) = FaceNodes % y(1)
      Surface(3) = 0.0_dp

      Normal(1) = Surface(2) - FaceNodes % y(2)
      Normal(2) = FaceNodes % x(2) - Surface(1)
      Normal(3) = 0.0_dp      
    END IF

    ! Project of the line to the face normal
    Rproj = SUM( (Rfin - Rinit) * Normal )
    
    IF( ABS( Rproj ) < TINY( Rproj ) ) THEN
      ! if the intersection cannot be defined make it an impossible one
      Lambda = -HUGE( Lambda ) 
    ELSE
      Lambda = SUM( ( Surface - Rinit ) * Normal ) / Rproj
    END IF

    IF( FaceElement % NDofs == 4 ) THEN
      IF( third == 3 ) THEN
        third = 4
	Lambda0 = Lambda
        GOTO 100
      END IF
      IF( ABS( Lambda0 ) < ABS( Lambda) ) THEN
        Lambda = Lambda0 
      END IF
   END IF


  END FUNCTION LineFaceIntersection
  

!---------------------------------------------------------------------------
!> This subroutine performs a similar test as above using slightly different 
!> strategy.
!---------------------------------------------------------------------------
  FUNCTION LineFaceIntersection2(FaceElement,FaceNodes,Rinit,Rfin,Intersect) RESULT ( Lambda ) 

    TYPE(Nodes_t) :: FaceNodes
    TYPE(Element_t), POINTER   :: FaceElement
    REAL(KIND=dp) :: Rinit(3), Rfin(3),Lambda
    LOGICAL :: Intersect
!----------------------------------------------------------------------------
    REAL (KIND=dp) :: A(3,3),B(3),C(3),Eps,Eps2,Eps3,detA,absA,ds
    INTEGER :: split, i, n, notriangles, triangle, ElemDim

    Eps = EPSILON( Eps )
    Eps2 = SQRT(TINY(Eps2))    
    Eps3 = 1.0d-12
    Lambda = -HUGE( Lambda )
    Intersect = .FALSE.
    ElemDim = FaceElement % TYPE % DIMENSION 

    ! Then solve the exact points of intersection from a 3x3 or 2x2 linear system
    !--------------------------------------------------------------------------
    IF( ElemDim == 2 ) THEN
      n = FaceElement % NDofs
      ! In 3D rectangular faces are treated as two triangles
      IF( n == 4 .OR. n == 8 .OR. n == 9 ) THEN
        notriangles = 2
      ELSE
        notriangles = 1
      END IF

      DO triangle=1,notriangles
          
        A(1:3,1) = Rfin(1:3) - Rinit(1:3)
        
        IF(triangle == 1) THEN
          A(1,2) = FaceNodes % x(1) - FaceNodes % x(2)
          A(2,2) = FaceNodes % y(1) - FaceNodes % y(2)
          A(3,2) = FaceNodes % z(1) - FaceNodes % z(2)
        ELSE 
          A(1,2) = FaceNodes % x(1) - FaceNodes % x(4)
          A(2,2) = FaceNodes % y(1) - FaceNodes % y(4)
          A(3,2) = FaceNodes % z(1) - FaceNodes % z(4)
        END IF

        A(1,3) = FaceNodes % x(1) - FaceNodes % x(3)
        A(2,3) = FaceNodes % y(1) - FaceNodes % y(3)
        A(3,3) = FaceNodes % z(1) - FaceNodes % z(3)
        
        ! Check for linearly dependent vectors
        detA = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
             - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
             + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
        absA = SUM(ABS(A(1,1:3))) * SUM(ABS(A(2,1:3))) * SUM(ABS(A(3,1:3))) 

        IF(ABS(detA) <= eps * absA + Eps2) CYCLE
!        print *,'detA',detA

        B(1) = FaceNodes % x(1) - Rinit(1)
        B(2) = FaceNodes % y(1) - Rinit(2)
        B(3) = FaceNodes % z(1) - Rinit(3)
        
        CALL InvertMatrix( A,3 )
        C(1:3) = MATMUL( A(1:3,1:3),B(1:3) )
        
        IF( ANY(C(2:3) < -Eps3) .OR. ANY(C(2:3) > 1.0_dp + Eps3 ) ) CYCLE
        IF( C(2)+C(3) > 1.0_dp + Eps3 ) CYCLE

        ! Relate the point of intersection to local coordinates
        !IF(corners < 4) THEN
        !  u = C(2)
        !  v = C(3)
        !ELSE IF(corners == 4 .AND. split == 0) THEN
        !  u = 2*(C(2)+C(3))-1
        !  v = 2*C(3)-1
        !ELSE 
        !  ! For the 2nd split of the rectangle the local coordinates switched
        !  v = 2*(C(2)+C(3))-1
        !  u = 2*C(3)-1        
        !END IF
        
        Intersect = .TRUE.
        Lambda = C(1)
        EXIT
 
      END DO
    ELSE
      ! In 2D the intersection is between two lines
      
      A(1:2,1) = Rfin(1:2) - Rinit(1:2)
      A(1,2) = FaceNodes % x(1) - FaceNodes % x(2)
      A(2,2) = FaceNodes % y(1) - FaceNodes % y(2)

      detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
      absA = SUM(ABS(A(1,1:2))) * SUM(ABS(A(2,1:2)))

      ! Lines are almost parallel => no intersection possible
      IF(ABS(detA) <= eps * absA + Eps2) RETURN

      B(1) = FaceNodes % x(1) - Rinit(1)
      B(2) = FaceNodes % y(1) - Rinit(2)

      CALL InvertMatrix( A,2 )
      C(1:2) = MATMUL(A(1:2,1:2),B(1:2))
     
      IF(C(2) < -Eps3 .OR. C(2) > 1.0_dp + Eps3 ) RETURN

      Intersect = .TRUE.
      Lambda = C(1)

!      u = -1.0d0 + 2.0d0 * C(2)

    END IF

!    IF(.NOT. Inside) RETURN

!    stat = ElementInfo( Element, FaceNodes, U, V, W, SqrtElementMetric, &
!        Basis, dBasisdx )
    
!    Weights(1:n) = Basis(1:n)
!    MaxInd = 1
!    DO i=2,n
!      IF(Weights(MaxInd) < Weights(i)) MaxInd = i
!    END DO

  END FUNCTION LineFaceIntersection2
  
 

!---------------------------------------------------------------------------
!> This subroutine computes the signed distance of a point from a surface.
!---------------------------------------------------------------------------
  FUNCTION PointFaceDistance(BoundaryElement,BoundaryNodes,&
      Coord,Normal,u0,v0) RESULT ( Dist )
!---------------------------------------------------------------------------
    TYPE(Nodes_t) :: BoundaryNodes
    TYPE(Element_t), POINTER   :: BoundaryElement
    REAL(KIND=dp) :: Coord(3),Normal(3)
    REAL(KIND=dp),OPTIONAL :: u0,v0
    REAL(KIND=dp) :: Dist

    REAL (KIND=dp) :: Surface(3),t1(3),t2(3),u,v

    ! For higher order elements this may be a necessity
    IF( PRESENT( u0 ) .AND. PRESENT(v0) ) THEN
      u = u0
      v = v0
      Surface = SurfaceVector( BoundaryElement, BoundaryNodes, u, v )
    ELSE
      u = 0.0_dp
      v = 0.0_dp

      ! Any point known to be at the surface, even corner node
      Surface(1) = BoundaryNodes % x(1)
      Surface(2) = BoundaryNodes % y(1)
      Surface(3) = BoundaryNodes % z(1)
    END IF

    Normal = NormalVector( BoundaryElement, BoundaryNodes, u, v, .TRUE. )

    ! Project of the line to the face normal
    Dist = SUM( (Surface - Coord ) * Normal ) 
END FUNCTION PointFaceDistance



!------------------------------------------------------------------------------
!> Convert global coordinates x,y,z inside element to local coordinates
!> u,v,w of the element.
!> @todo Change to support p elements
!------------------------------------------------------------------------------
  SUBROUTINE GlobalToLocal( u,v,w,x,y,z,Element,ElementNodes )
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp) :: x,y,z,u,v,w
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: MaxIter = 50
    INTEGER :: i,n
    REAL(KIND=dp) :: r,s,t,delta(3),prevdelta(3),J(3,3),J1(3,2),det,swap,acc,err
    LOGICAL :: Converged
!------------------------------------------------------------------------------

    u = 0._dp
    v = 0._dp
    w = 0._dp
    IF (Element % TYPE % DIMENSION==0) RETURN

    n = Element % TYPE % NumberOfNodes

    ! @todo Not supported yet
!   IF (ASSOCIATED(Element % PDefs)) THEN
!      CALL Fatal('GlobalToLocal','P elements not supported yet!')
!   END IF
    acc = EPSILON(1.0_dp)
    Converged = .FALSE.

     delta = 0._dp

!------------------------------------------------------------------------------
    DO i=1,Maxiter
!------------------------------------------------------------------------------
      r = InterpolateInElement(Element,ElementNodes % x(1:n),u,v,w) - x
      s = InterpolateInElement(Element,ElementNodes % y(1:n),u,v,w) - y
      t = InterpolateInElement(Element,ElementNodes % z(1:n),u,v,w) - z

      err = r**2 + s**2 + t**2 

      IF ( err < acc ) THEN
        Converged = .TRUE.
        EXIT
      END IF

      prevdelta = delta
      delta = 0.d0

      SELECT CASE( Element % TYPE % DIMENSION )
      CASE(1)

        J(1,1) = FirstDerivative1D( Element, ElementNodes % x, u )
        J(2,1) = FirstDerivative1D( Element, ElementNodes % y, u )
        J(3,1) = FirstDerivative1D( Element, ElementNodes % z, u )

        det = SUM( J(1:3,1)**2 )
        delta(1) = (r*J(1,1)+s*J(2,1)+t*J(3,1))/det

      CASE(2)

         J(1,1) = FirstDerivativeInU2D( Element, ElementNodes % x,u,v )
         J(1,2) = FirstDerivativeInV2D( Element, ElementNodes % x,u,v )
         J(2,1) = FirstDerivativeInU2D( Element, ElementNodes % y,u,v )
         J(2,2) = FirstDerivativeInV2D( Element, ElementNodes % y,u,v )

        SELECT CASE( CoordinateSystemDimension() )
           CASE(3)
              J(3,1) = FirstDerivativeInU2D( Element, ElementNodes % z, u, v )
              J(3,2) = FirstDerivativeInV2D( Element, ElementNodes % z, u, v )

              delta(1) = r
              delta(2) = s
              delta(3) = t
              delta(1:2) = MATMUL( TRANSPOSE(J(1:3,1:2)), delta )
              r = delta(1)
              s = delta(2)

              J(1:2,1:2) = MATMUL( TRANSPOSE(J(1:3,1:2)), J(1:3,1:2) )
              delta(3)   = 0.0d0
         END SELECT

         CALL SolveLinSys2x2( J(1:2,1:2), delta(1:2), [ r, s] )

      CASE(3)
        J(1,1) = FirstDerivativeInU3D( Element, ElementNodes % x, u, v, w )
        J(1,2) = FirstDerivativeInV3D( Element, ElementNodes % x, u, v, w )
        J(1,3) = FirstDerivativeInW3D( Element, ElementNodes % x, u, v, w )

        J(2,1) = FirstDerivativeInU3D( Element, ElementNodes % y, u, v, w )
        J(2,2) = FirstDerivativeInV3D( Element, ElementNodes % y, u, v, w )
        J(2,3) = FirstDerivativeInW3D( Element, ElementNodes % y, u, v, w )

        J(3,1) = FirstDerivativeInU3D( Element, ElementNodes % z, u, v, w )
        J(3,2) = FirstDerivativeInV3D( Element, ElementNodes % z, u, v, w )
        J(3,3) = FirstDerivativeInW3D( Element, ElementNodes % z, u, v, w )

        CALL SolveLinSys3x3( J, delta, [ r, s, t ] )

      END SELECT

      IF( i > 10 ) THEN
        ! If the same values is suggested over and over again, then exit
        ! This may be a sign that the node is off-plane and cannot be 
        ! described within the element.
        IF( SUM( ABS( delta - prevdelta ) ) < acc ) EXIT

        ! Use sloppier criteria when iteration still unsuccesfull
        IF( i > 20 ) THEN
          IF( SUM( ABS( delta - prevdelta ) ) < SQRT( acc ) ) EXIT         
        END IF

        ! If the iteration does not proceed try with some relaxation
        delta = 0.5_dp * delta 
      END IF

      u = u - delta(1)
      v = v - delta(2)
      w = w - delta(3)


!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

    IF ( .NOT. Converged ) THEN        
      IF( err > SQRT( acc ) ) THEN
        IF( i > MaxIter ) THEN	
          CALL Warn( 'GlobalToLocal', 'did not converge.')
          PRINT *,'rst',i,r,s,t
          PRINT *,'err',err,acc,SQRT(acc)
          PRINT *,'delta',delta,prevdelta
          PRINT *,'uvw',u,v,w
          PRINT *,'code',Element % TYPE % ElementCode
          PRINT *,'x:',x,ElementNodes % x(1:n)
          PRINT *,'y:',y,ElementNodes % y(1:n)
          PRINT *,'z:',z,ElementNodes % z(1:n)
        ELSE
!          CALL Warn( 'GlobalToLocal', 'Node may be out of element')
!          PRINT *,'rst',i,r,s,t,acc
        END IF
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GlobalToLocal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE InvertMatrix3x3( G,GI,detG )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: G(3,3),GI(3,3)
    REAL(KIND=dp) :: detG, s
!------------------------------------------------------------------------------
    s = 1.0 / DetG
    
    GI(1,1) =  s * (G(2,2)*G(3,3) - G(3,2)*G(2,3));
    GI(2,1) = -s * (G(2,1)*G(3,3) - G(3,1)*G(2,3));
    GI(3,1) =  s * (G(2,1)*G(3,2) - G(3,1)*G(2,2));
    
    GI(1,2) = -s * (G(1,2)*G(3,3) - G(3,2)*G(1,3));
    GI(2,2) =  s * (G(1,1)*G(3,3) - G(3,1)*G(1,3));
    GI(3,2) = -s * (G(1,1)*G(3,2) - G(3,1)*G(1,2));

    GI(1,3) =  s * (G(1,2)*G(2,3) - G(2,2)*G(1,3));
    GI(2,3) = -s * (G(1,1)*G(2,3) - G(2,1)*G(1,3));
    GI(3,3) =  s * (G(1,1)*G(2,2) - G(2,1)*G(1,2));
!------------------------------------------------------------------------------
  END SUBROUTINE InvertMatrix3x3
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>     Given element and its face map (for some triangular face of element ), 
!>     this routine returns global direction of triangle face so that 
!>     functions are continuous over element boundaries
!------------------------------------------------------------------------------
  FUNCTION getTriangleFaceDirection( Element, FaceMap ) RESULT(globalDir)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to get direction to
!
!    INTEGER :: FaceMap(3)
!      INPUT: Element triangular face map
!
!  FUNCTION VALUE:
!    INTEGER :: globalDir(3)
!       Global direction of triangular face as local node numbers.
!    
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t) :: Element
    INTEGER :: i, FaceMap(3), globalDir(3), nodes(3)

    nodes = 0
    
    ! Put global nodes of face into sorted order
    nodes(1:3) = Element % NodeIndexes( FaceMap )
    CALL sort(3, nodes)
    
    globalDir = 0
    ! Find local numbers of sorted nodes. These local nodes 
    ! span continuous functions over element boundaries
    DO i=1,Element % TYPE % NumberOfNodes
       IF (nodes(1) == Element % NodeIndexes(i)) THEN
          globalDir(1) = i
       ELSE IF (nodes(2) == Element % NodeIndexes(i)) THEN
          globalDir(2) = i
       ELSE IF (nodes(3) == Element % NodeIndexes(i)) THEN
          globalDir(3) = i
       END IF
    END DO
  END FUNCTION getTriangleFaceDirection


!------------------------------------------------------------------------------
!>     Given element and its face map (for some square face of element ), 
!>     this routine returns global direction of square face so that 
!>     functions are continuous over element boundaries
!------------------------------------------------------------------------------
  FUNCTION getSquareFaceDirection( Element, FaceMap ) RESULT(globalDir)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to get direction to
!
!    INTEGER :: FaceMap(4)
!      INPUT: Element square face map
!
!  FUNCTION VALUE:
!    INTEGER :: globalDir(3)
!       Global direction of square face as local node numbers.
!    
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t) :: Element
    INTEGER :: i, A,B,C,D, FaceMap(4), globalDir(4), nodes(4), minGlobal

    ! Get global nodes 
    nodes(1:4) = Element % NodeIndexes( FaceMap )
    ! Find min global node
    minGlobal = nodes(1)
    A = 1
    DO i=2,4
       IF (nodes(i) < minGlobal) THEN
          A = i
          minGlobal = nodes(i)
       END IF
    END DO

    ! Now choose node B as the smallest node NEXT to min node
    B = MOD(A,4)+1
    C = MOD(A+3,4)
    IF (C == 0) C = 4
    D = MOD(A+2,4)
    IF (D == 0) D = 4
    IF (nodes(B) > nodes(C)) THEN
       i = B
       B = C
       C = i
    END IF

    ! Finally find local numbers of nodes A,B and C. They uniquely
    ! define a global face so that basis functions are continuous 
    ! over element boundaries
    globalDir = 0
    DO i=1,Element % TYPE % NumberOfNodes
       IF (nodes(A) == Element % NodeIndexes(i)) THEN
          globalDir(1) = i
       ELSE IF (nodes(B) == Element % NodeIndexes(i)) THEN
          globalDir(2) = i
       ELSE IF (nodes(C) == Element % NodeIndexes(i)) THEN
          globalDir(4) = i
       ELSE IF (nodes(D) == Element % NodeIndexes(i)) THEN
          globalDir(3) = i
       END IF
    END DO
  END FUNCTION getSquareFaceDirection


!------------------------------------------------------------------------------
!>     Function checks if given local numbering of a square face
!>     is legal for wedge element
!------------------------------------------------------------------------------
  FUNCTION wedgeOrdering( ordering ) RESULT(retVal)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!
!    INTEGER :: ordering(4)
!      INPUT: Local ordering of a wedge square face
!
!  FUNCTION VALUE:
!    INTEGER :: retVal
!       .TRUE. if given ordering is legal for wedge square face,
!       .FALSE. otherwise
!    
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, DIMENSION(4), INTENT(IN) :: ordering
    LOGICAL :: retVal

    retVal = .FALSE.
    IF ((ordering(1) >= 1 .AND. ordering(1) <= 3 .AND.&
         ordering(2) >= 1 .AND. ordering(2) <= 3) .OR. &
       (ordering(1) >= 4 .AND. ordering(1) <= 6 .AND.&
       ordering(2) >= 4 .AND. ordering(2) <= 6)) THEN
       retVal = .TRUE.
    END IF
  END FUNCTION wedgeOrdering

  !---------------------------------------------------------
  !> Computes the 3D rotation matrix for a given 
  !> surface normal vector
  !---------------------------------------------------------
  FUNCTION ComputeRotationMatrix(PlaneVector) RESULT ( RotMat )

    REAL(KIND=dp) :: PlaneVector(3), RotMat(3,3), ex(3), ey(3), ez(3)
    INTEGER :: i, MinIndex, MidIndex, MaxIndex

    !Ensure PlaneVector is the unit normal
    PlaneVector = PlaneVector / SQRT( SUM(PlaneVector ** 2) )
    
    !The new z-axis is normal to the defined surface
    ez = PlaneVector

    MaxIndex = MAXLOC(ABS(ez),1)
    MinIndex = MINLOC(ABS(ez),1)

    !Special case when calving front perfectly aligned to either
    ! x or y axis. In this case, make minindex = 3 (ex points upwards)
    IF(ABS(ez(3)) == ABS(ez(2)) .OR. ABS(ez(3)) == ABS(ez(1))) &
         MinIndex = 3

    DO i=1,3
       IF(i == MaxIndex .OR. i == MinIndex) CYCLE
       MidIndex = i
    END DO

    ex(MinIndex) = 1.0
    ex(MidIndex) = 0.0
    
    ex(MaxIndex) = -ez(MinIndex)/ez(MaxIndex)
    ex = ex / SQRT( SUM(ex ** 2) )

    !The new y-axis is orthogonal to new x and z axes
    ey = CrossProduct(ez, ex)
    ey = ey / SQRT( SUM(ey ** 2) ) !just in case...

    RotMat(1,:) = ex
    RotMat(2,:) = ey
    RotMat(3,:) = ez

  END FUNCTION ComputeRotationMatrix

END MODULE ElementDescription


!> \}
